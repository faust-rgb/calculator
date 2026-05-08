// ============================================================================
// 节点驻留与结构键缓存
// ============================================================================
//
// 本文件实现符号表达式节点的驻留机制：
// - 相同结构的表达式节点共享同一内存实例
// - 使用 LRU 缓存管理驻留节点（最大 8192 个）
// - 减少内存分配，加速结构比较
// ============================================================================

#include "symbolic/symbolic_expression_internal.h"

#include "core/format_utils.h"
#include "math/mymath.h"
#include "math/mymath_dual.h"

#include <algorithm>
#include <cctype>
#include <cstring>
#include <iomanip>
#include <list>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <mutex>

namespace symbolic_expression_internal {

namespace {

std::string long_double_bits_key(long double value) {
    unsigned char bytes[sizeof(long double)] = {};
    std::memcpy(bytes, &value, sizeof(value));

    std::ostringstream out;
    out << std::hex << std::setfill('0');
    for (unsigned char byte : bytes) {
        out << std::setw(2) << static_cast<unsigned int>(byte);
    }
    return out.str();
}

std::string scalar_structural_number_key(Scalar value) {
    return long_double_bits_key(value.hi) + ":" + long_double_bits_key(value.lo);
}

} // namespace

int& mutable_display_precision() {
    static int precision = 12;
    return precision;
}

int clamp_display_precision(int precision) {
    return std::clamp(precision, 1, 17);
}


// ============================================================================
// 节点驻留（Interning）
// ============================================================================

/**
 * @brief 驻留表达式节点
 *
 * 节点驻留确保相同结构的表达式共享同一实例：
 * - 减少内存分配：相同子表达式只存储一份
 * - 加速比较：指针比较可判断结构相等
 * - 支持缓存：结构键可作为缓存键
 *
 * 使用 LRU 策略管理驻留池，最大 8192 个节点（可自适应调整）。
 * 当池满时，淘汰最久未使用的节点。
 *
 * @param node 待驻留的节点
 * @return 驻留后的节点指针（可能是已存在的实例）
 */
std::shared_ptr<SymbolicExpression::Node> intern_node(
    std::shared_ptr<SymbolicExpression::Node> node) {
    // 驻留池：存储结构键到节点的弱引用
    using InternEntry =
        std::pair<std::string, std::weak_ptr<SymbolicExpression::Node>>;

    // 全局驻留池，带互斥锁保证线程安全
    static std::list<InternEntry> interned_order;  // LRU 顺序
    static std::unordered_map<std::string,
                                           std::list<InternEntry>::iterator>
        interned_nodes;  // 快速查找
    static std::mutex intern_mutex;
    static std::size_t kMaxInternedNodes = 8192;  // 最大驻留数

    std::lock_guard<std::mutex> lock(intern_mutex);

    // 计算结构键
    const std::string key = node_structural_key(node);

    // 查找是否已驻留
    const auto found = interned_nodes.find(key);
    if (found != interned_nodes.end()) {
        // 尝试提升弱引用为共享指针
        if (std::shared_ptr<SymbolicExpression::Node> existing = found->second->second.lock()) {
            // 移动到 LRU 队列前端
            interned_order.splice(interned_order.begin(), interned_order, found->second);
            return existing;
        }
        // 弱引用已失效，清理条目
        interned_order.erase(found->second);
        interned_nodes.erase(found);
    }

    // 池满时清理过期条目并淘汰 LRU
    if (interned_nodes.size() >= kMaxInternedNodes) {
        // 增量清理
        int checked = 0;
        const int kMaxCheckPerInsert = 32; 
        
        for (auto it = interned_order.rbegin(); it != interned_order.rend() && checked < kMaxCheckPerInsert; ) {
            if (it->second.expired()) {
                auto erase_it = std::next(it).base();
                interned_nodes.erase(erase_it->first);
                interned_order.erase(erase_it);
                it = interned_order.rbegin();
                checked++;
            } else {
                ++it;
            }
        }

        // 自适应调整：如果几乎所有节点都在使用中，扩容
        if (interned_nodes.size() >= kMaxInternedNodes) {
             if (kMaxInternedNodes < 65536) {
                 kMaxInternedNodes *= 2;
             } else {
                // 强制移除末尾最旧的条目
                while (interned_nodes.size() >= kMaxInternedNodes && !interned_order.empty()) {
                    interned_nodes.erase(interned_order.back().first);
                    interned_order.pop_back();
                }
             }
        }
    }

    // 插入新条目
    interned_order.push_front({key, node});
    interned_nodes[key] = interned_order.begin();
    return node;
}

// ============================================================================
// 简化结果缓存
// ============================================================================

SymbolicExpressionLruCache::SymbolicExpressionLruCache(std::size_t capacity)
    : capacity_(capacity) {}

bool SymbolicExpressionLruCache::get(const std::string& key, SymbolicExpression* value) {
    const auto found = index_.find(key);
    if (found == index_.end()) {
        return false;
    }
    entries_.splice(entries_.begin(), entries_, found->second);
    *value = found->second->second;
    return true;
}

void SymbolicExpressionLruCache::put(const std::string& key,
                                     const SymbolicExpression& value) {
    const auto found = index_.find(key);
    if (found != index_.end()) {
        found->second->second = value;
        entries_.splice(entries_.begin(), entries_, found->second);
        return;
    }

    entries_.push_front({key, value});
    index_[key] = entries_.begin();
    while (entries_.size() > capacity_) {
        index_.erase(entries_.back().first);
        entries_.pop_back();
    }
}

// ============================================================================
// 数值格式化
// ============================================================================

/**
 * @brief 格式化数值为字符串
 *
 * 格式化策略：
 * 1. 接近零：输出 "0"
 * 2. 接近整数：输出整数形式
 * 3. 可近似为分数：输出分数形式（如 "1/2"）
 * 4. 否则：输出浮点数（12 位精度）
 *
 * @param value 数值
 * @return 格式化字符串
 */
std::string format_number(Scalar value) {
    if (mymath::is_near_zero(value, kFormatEps)) {
        return "0";
    }
    if (value < Scalar(0.0L)) {
        return "-" + format_number(-value);
    }

    if (mymath::is_integer(value, Scalar(1e-10L))) {
        return std::to_string(static_cast<long long>(static_cast<long double>(value.to_long_double()) + 0.5));
    }

    if (mutable_display_precision() < 12) {
        return format_decimal(value);
    }

    long long n, d;
    long double ld_val = static_cast<long double>(value.to_long_double());
    if (mymath::approximate_fraction(ld_val, &n, &d, 999, 1e-10)) {
        if (d == 1) return std::to_string(n);
        return std::to_string(n) + "/" + std::to_string(d);
    }

    std::ostringstream out;
    out.precision(mutable_display_precision());
    out << value;
    return out.str();
}

// ============================================================================
// 节点构造函数
// ============================================================================

/** @brief 创建数值节点（带驻留） */
std::shared_ptr<SymbolicExpression::Node> make_number(Scalar value) {
    return intern_node(std::make_shared<SymbolicExpression::Node>(value));
}

/** @brief 创建变量节点（带驻留） */
std::shared_ptr<SymbolicExpression::Node> make_variable(const std::string& name) {
    std::shared_ptr<SymbolicExpression::Node> node =
        std::make_shared<SymbolicExpression::Node>();
    if (name == "pi") {
        node->type = NodeType::kPi;
    } else if (name == "e") {
        node->type = NodeType::kE;
    } else {
        node->type = NodeType::kVariable;
        node->text = name;
    }
    return intern_node(node);
}

/** @brief 创建无穷大节点（带驻留） */
std::shared_ptr<SymbolicExpression::Node> make_infinity(bool positive) {
    std::shared_ptr<SymbolicExpression::Node> node =
        std::make_shared<SymbolicExpression::Node>();
    node->type = NodeType::kInfinity;
    node->number_value = positive ? 1.0L : -1.0L;  // 1.0L 表示 +inf, -1.0L 表示 -inf
    return intern_node(node);
}

/**
 * @brief 创建一元节点（取负或函数调用）
 * @param type 节点类型
 * @param operand 操作数
 * @param text 函数名（仅 kFunction 使用）
 */
std::shared_ptr<SymbolicExpression::Node> make_unary(NodeType type,
                                                     std::shared_ptr<SymbolicExpression::Node> operand,
                                                     const std::string& text) {
    std::shared_ptr<SymbolicExpression::Node> node =
        std::make_shared<SymbolicExpression::Node>();
    node->type = type;
    node->left = std::move(operand);
    node->text = text;
    return intern_node(node);
}

/**
 * @brief 创建二元运算节点
 * @param type 运算类型
 * @param left 左操作数
 * @param right 右操作数
 */
std::shared_ptr<SymbolicExpression::Node> make_binary(NodeType type,
                                                      std::shared_ptr<SymbolicExpression::Node> left,
                                                      std::shared_ptr<SymbolicExpression::Node> right) {
    std::shared_ptr<SymbolicExpression::Node> node =
        std::make_shared<SymbolicExpression::Node>();
    node->type = type;
    node->left = std::move(left);
    node->right = std::move(right);
    return intern_node(node);
}

// ============================================================================
// 运算优先级
// ============================================================================

/**
 * @brief 获取节点的运算优先级
 *
 * 优先级用于决定输出时是否需要括号：
 * - 较低优先级的子表达式需要括号
 * - 例如 (a + b) * c 中，a + b 需要括号
 *
 * 优先级表（从低到高）：
 * - 1: 加法、减法
 * - 2: 乘法、除法
 * - 3: 幂运算
 * - 4: 取负
 * - 5: 函数、变量、数值
 */
int precedence(const std::shared_ptr<SymbolicExpression::Node>& node) {
    switch (node->type) {
        case NodeType::kAdd:
        case NodeType::kSubtract:
            return 1;
        case NodeType::kMultiply:
        case NodeType::kDivide:
            return 2;
        case NodeType::kPower:
            return 3;
        case NodeType::kNegate:
            return 4;
        case NodeType::kFunction:
        case NodeType::kNumber:
        case NodeType::kVariable:
        case NodeType::kPi:
        case NodeType::kE:
        case NodeType::kInfinity:
        case NodeType::kVector:
        case NodeType::kTensor:
        case NodeType::kDifferentialOp:
        case NodeType::kRootOf:
            return 5;
    }
    return 5;
}

// ============================================================================
// 字符串转换
// ============================================================================

/**
 * @brief 将表达式节点转换为字符串
 *
 * 递归遍历表达式树，生成可读的字符串表示。
 * 根据优先级智能插入括号，避免冗余。
 *
 * @param node 表达式节点
 * @param parent_precedence 父节点的优先级
 * @return 表达式字符串
 */
std::string to_string_impl(const std::shared_ptr<SymbolicExpression::Node>& node, int parent_precedence) {
    static thread_local int depth = 0;
    static constexpr int kMaxDepth = 1000;
    if (++depth > kMaxDepth) {
        --depth;
        return "...";
    }

    struct DepthGuard {
        int* depth;
        ~DepthGuard() { if (depth) (*depth)--; }
    } guard{&depth};

    std::string text;
    switch (node->type) {
        case NodeType::kNumber:
            text = format_number(node->number_value);
            break;
        case NodeType::kVariable:
            text = node->text;
            break;
        case NodeType::kPi:
            text = "pi";
            break;
        case NodeType::kE:
            text = "e";
            break;
        case NodeType::kInfinity:
            text = (node->number_value > 0) ? "inf" : "-inf";
            break;
        case NodeType::kNegate:
            text = "-" + to_string_impl(node->left, precedence(node));
            break;
        case NodeType::kFunction:
            text = node->text + "(" + to_string_impl(node->left, 0) + ")";
            break;
        case NodeType::kAdd:
            // 特殊处理：将 a + (-b) 显示为 a - b
            if (node->right->type == NodeType::kNegate) {
                text = to_string_impl(node->left, precedence(node)) + " - " +
                       to_string_impl(node->right->left, precedence(node) + 1);
            } else {
                text = to_string_impl(node->left, precedence(node)) + " + " +
                       to_string_impl(node->right, precedence(node));
            }
            break;
        case NodeType::kSubtract:
            if (node->right->type == NodeType::kNumber &&
                node->right->number_value < Scalar(0)) {
                text = to_string_impl(node->left, precedence(node)) + " + " +
                       format_number(-node->right->number_value);
            } else {
                text = to_string_impl(node->left, precedence(node)) + " - " +
                       to_string_impl(node->right, precedence(node) + 1);
            }
            break;
        case NodeType::kMultiply:
            text = to_string_impl(node->left, precedence(node)) + " * " +
                   to_string_impl(node->right, precedence(node));
            break;
        case NodeType::kDivide:
            text = to_string_impl(node->left, precedence(node)) + " / " +
                   to_string_impl(node->right, precedence(node) + 1);
            break;
        case NodeType::kPower:
            text = to_string_impl(node->left, precedence(node)) + " ^ " +
                   to_string_impl(node->right, precedence(node));
            break;
        case NodeType::kVector: {
            text = "[";
            for (std::size_t i = 0; i < node->children.size(); ++i) {
                if (i != 0) text += ", ";
                text += to_string_impl(node->children[i], 0);
            }
            text += "]";
            break;
        }
        case NodeType::kTensor: {
            text = "[";
            for (std::size_t i = 0; i < node->children.size(); ++i) {
                if (i != 0) text += "; ";
                text += to_string_impl(node->children[i], 0);
            }
            text += "]";
            break;
        }
        case NodeType::kDifferentialOp:
            text = node->text + "(" + to_string_impl(node->left, 0) + ")";
            break;
        case NodeType::kRootOf: {
            // RootOf(poly, var, index)
            text = "RootOf(";
            if (!node->children.empty()) {
                text += to_string_impl(node->children[0], 0);
            }
            text += ", " + node->text;
            text += ", " + std::to_string(static_cast<int>(node->number_value));
            text += ")";
            break;
        }
    }

    // 分数作幂运算底数时需要括号：如 (1/2)^x
    if (node->type == NodeType::kNumber &&
        text.find('/') != std::string::npos &&
        parent_precedence >= 3) {
        return "(" + text + ")";
    }

    // 优先级低于父节点时需要括号
    if (precedence(node) < parent_precedence) {
        return "(" + text + ")";
    }
    return text;
}

// ============================================================================
// 结构键计算
// ============================================================================

/**
 * @brief 计算节点的结构键
 *
 * 结构键是唯一标识表达式结构的字符串：
 * - 相同结构的表达式具有相同的结构键
 * - 用于节点驻留、缓存查找和表达式比较
 *
 * 格式示例：
 * - 数值：N(3.14)
 * - 变量：V(x)
 * - 加法：ADD(N(1),V(x))
 * - 函数：F(sin:V(x))
 *
 * 结果缓存在节点的 structural_key_cache 字段中。
 */
std::string node_structural_key(const std::shared_ptr<SymbolicExpression::Node>& node) {
    if (!node->structural_key_cache.empty()) {
        return node->structural_key_cache;
    }

    std::string key;
    switch (node->type) {
        case NodeType::kNumber:
            // 使用原始位模式保证内部键的唯一性，避免受显示格式和浮点文本格式化影响。
            key = "N(" + scalar_structural_number_key(node->number_value) + ")";
            break;
        case NodeType::kVariable:
            key = "V(" + node->text + ")";
            break;
        case NodeType::kPi:
            key = "PI";
            break;
        case NodeType::kE:
            key = "E";
            break;
        case NodeType::kInfinity:
            key = (node->number_value > 0) ? "INF" : "NINF";
            break;
        case NodeType::kNegate:
            key = "NEG(" + node_structural_key(node->left) + ")";
            break;
        case NodeType::kFunction:
            key = "F(" + node->text + ":" + node_structural_key(node->left) + ")";
            break;
        case NodeType::kAdd:
            key = "ADD(" + node_structural_key(node->left) + "," +
                  node_structural_key(node->right) + ")";
            break;
        case NodeType::kSubtract:
            key = "SUB(" + node_structural_key(node->left) + "," +
                  node_structural_key(node->right) + ")";
            break;
        case NodeType::kMultiply:
            key = "MUL(" + node_structural_key(node->left) + "," +
                  node_structural_key(node->right) + ")";
            break;
        case NodeType::kDivide:
            key = "DIV(" + node_structural_key(node->left) + "," +
                  node_structural_key(node->right) + ")";
            break;
        case NodeType::kPower:
            key = "POW(" + node_structural_key(node->left) + "," +
                  node_structural_key(node->right) + ")";
            break;
        case NodeType::kVector: {
            key = "VEC(";
            for (std::size_t i = 0; i < node->children.size(); ++i) {
                if (i > 0) key += ",";
                key += node_structural_key(node->children[i]);
            }
            key += ")";
            break;
        }
        case NodeType::kTensor: {
            key = "TEN(";
            for (std::size_t i = 0; i < node->children.size(); ++i) {
                if (i > 0) key += ";";
                key += node_structural_key(node->children[i]);
            }
            key += ")";
            break;
        }
        case NodeType::kDifferentialOp:
            key = "DOP(" + node->text + ":" + node_structural_key(node->left) + ")";
            break;
        case NodeType::kRootOf:
            key = "ROOTOF(" + node->text + ")";
            break;
    }
    node->structural_key_cache = key;
    return node->structural_key_cache;
}

}  // namespace symbolic_expression_internal
