// ============================================================================
// 符号表达式解析与节点管理模块
// ============================================================================
//
// 本文件实现符号表达式的核心基础设施：
//
// 1. 节点驻留（Interning）
//    - 相同结构的表达式节点共享同一内存实例
//    - 使用 LRU 缓存管理驻留节点（最大 8192 个）
//    - 减少内存分配，加速结构比较
//
// 2. 结构键缓存
//    - 每个节点缓存其结构字符串键
//    - 用于快速比较、缓存查找和规范化排序
//
// 3. 表达式解析器
//    - 递归下降解析器，支持标准数学语法
//    - 运算优先级：加减 < 乘除 < 幂 < 一元 < 函数/原子
//    - 支持变量、常量（pi, e）和函数调用
//
// 4. 字符串格式化
//    - 将表达式树转换为可读字符串
//    - 智能括号插入，避免冗余括号
//    - 数值格式化：整数、分数、浮点数
//
// 5. 简化缓存
//    - LRU 缓存简化结果（最大 4096 条目）
//    - 线程局部存储，避免锁竞争
// ============================================================================

#include "symbolic_expression_internal.h"

#include "mymath.h"

#include <algorithm>
#include <cctype>
#include <list>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace symbolic_expression_internal {

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
 * 使用 LRU 策略管理驻留池，最大 8192 个节点。
 * 当池满时，淘汰最久未使用的节点。
 *
 * @param node 待驻留的节点
 * @return 驻留后的节点指针（可能是已存在的实例）
 */
std::shared_ptr<SymbolicExpression::Node> intern_node(
    std::shared_ptr<SymbolicExpression::Node> node) {
    // 驻留池：存储结构键到节点的弱引用
    // 使用弱引用允许节点在无外部引用时被释放
    using InternEntry =
        std::pair<std::string, std::weak_ptr<SymbolicExpression::Node>>;

    // 线程局部的驻留池，避免锁竞争
    static thread_local std::list<InternEntry> interned_order;  // LRU 顺序
    static thread_local std::unordered_map<std::string,
                                           std::list<InternEntry>::iterator>
        interned_nodes;  // 快速查找
    static constexpr std::size_t kMaxInternedNodes = 8192;  // 最大驻留数

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
        // 增量清理：从末尾开始检查一定数量的条目，释放失效的弱引用或直接淘汰最旧条目
        // 这样避免了 O(N) 的全量扫描，保证了单次插入的性能平稳
        int checked = 0;
        const int kMaxCheckPerInsert = 32; // 每次插入最多检查的过期条目数
        
        for (auto it = interned_order.rbegin(); it != interned_order.rend() && checked < kMaxCheckPerInsert; ) {
            if (it->second.expired()) {
                auto erase_it = std::next(it).base(); // 转换为正向迭代器
                interned_nodes.erase(erase_it->first);
                interned_order.erase(erase_it);
                it = interned_order.rbegin(); // 结构改变，重新开始（或者简单的 --it）
                checked++;
            } else {
                ++it;
            }
        }

        // 如果清理后依然满，则强制移除末尾最旧的条目
        while (interned_nodes.size() >= kMaxInternedNodes && !interned_order.empty()) {
            interned_nodes.erase(interned_order.back().first);
            interned_order.pop_back();
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

/**
 * @class SymbolicExpressionLruCache
 * @brief 表达式 LRU 缓存
 *
 * 用于缓存 simplify() 和 derivative() 的结果。
 * 键为表达式的结构键，值为简化/求导后的表达式。
 *
 * 线程局部存储，每个线程独立缓存，避免锁竞争。
 */
class SymbolicExpressionLruCache {
public:
    explicit SymbolicExpressionLruCache(std::size_t capacity)
        : capacity_(capacity) {}

    /**
     * @brief 查找缓存
     * @param key 结构键
     * @param value 输出参数，找到时写入缓存值
     * @return true 如果命中缓存
     */
    bool get(const std::string& key, SymbolicExpression* value) {
        const auto found = index_.find(key);
        if (found == index_.end()) {
            return false;
        }
        // 移动到前端（LRU 更新）
        entries_.splice(entries_.begin(), entries_, found->second);
        *value = found->second->second;
        return true;
    }

    /**
     * @brief 插入或更新缓存
     * @param key 结构键
     * @param value 缓存值
     */
    void put(const std::string& key, const SymbolicExpression& value) {
        const auto found = index_.find(key);
        if (found != index_.end()) {
            // 更新已存在的条目
            found->second->second = value;
            entries_.splice(entries_.begin(), entries_, found->second);
            return;
        }

        // 插入新条目
        entries_.push_front({key, value});
        index_[key] = entries_.begin();

        // 超容量时淘汰 LRU
        while (entries_.size() > capacity_) {
            index_.erase(entries_.back().first);
            entries_.pop_back();
        }
    }

private:
    std::size_t capacity_ = 0;  // 缓存容量
    std::list<std::pair<std::string, SymbolicExpression>> entries_;  // LRU 队列
    std::unordered_map<std::string,
                       std::list<std::pair<std::string, SymbolicExpression>>::iterator>
        index_;  // 快速查找索引
};

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
std::string format_number(double value) {
    // 接近零归零，避免 -0 输出
    if (mymath::is_near_zero(value, kFormatEps)) {
        return "0";
    }

    // 接近整数按整数打印
    if (mymath::is_integer(value, 1e-10)) {
        long long rounded = static_cast<long long>(value >= 0.0 ? value + 0.5 : value - 0.5);
        return std::to_string(rounded);
    }

    // 尝试分数近似（分母最大 999）
    long long numerator = 0;
    long long denominator = 1;
    if (mymath::approximate_fraction(value,
                                     &numerator,
                                     &denominator,
                                     999,
                                     1e-10)) {
        if (value < 0.0) {
            numerator = -numerator;
        }
        if (denominator == 1) {
            return std::to_string(numerator);
        }
        return std::to_string(numerator) + "/" + std::to_string(denominator);
    }

    // 一般浮点数，12 位精度
    std::ostringstream out;
    out.precision(mutable_display_precision());
    out << value;
    return out.str();
}

// ============================================================================
// 节点构造函数
// ============================================================================

/** @brief 创建数值节点（带驻留） */
std::shared_ptr<SymbolicExpression::Node> make_number(double value) {
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
            text = to_string_impl(node->left, precedence(node)) + " - " +
                   to_string_impl(node->right, precedence(node) + 1);
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
    // 使用缓存避免重复计算
    if (!node->structural_key_cache.empty()) {
        return node->structural_key_cache;
    }

    std::string key;
    switch (node->type) {
        case NodeType::kNumber:
            key = "N(" + format_number(node->number_value) + ")";
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
    }
    node->structural_key_cache = key;
    return node->structural_key_cache;
}

// ============================================================================
// 表达式解析器
// ============================================================================

/**
 * @class Parser
 * @brief 递归下降表达式解析器
 *
 * 解析数学表达式字符串，构建表达式树。
 *
 * 文法规则：
 * ```
 * expression -> term (('+' | '-') term)*
 * term       -> unary (('*' | '/') unary)*
 * unary      -> ('+' | '-') unary | power
 * power      -> primary ('^' unary)?
 * primary    -> '(' expression ')' | identifier '(' expression ')' | identifier | number
 * ```
 *
 * 运算优先级（从低到高）：
 * 1. 加法、减法
 * 2. 乘法、除法
 * 3. 幂运算（右结合）
 * 4. 一元正负
 * 5. 函数调用、原子表达式
 */
class Parser {
public:
    explicit Parser(std::string source) : source_(std::move(source)) {}

    /**
     * @brief 解析表达式
     * @return 解析后的符号表达式
     * @throw std::runtime_error 当语法错误时
     */
    SymbolicExpression parse() {
        SymbolicExpression expression = parse_expression();
        skip_spaces();
        if (pos_ != source_.size()) {
            throw std::runtime_error("unexpected token near: " + source_.substr(pos_, 1));
        }
        return expression;
    }

private:
    std::string source_;
    std::size_t pos_ = 0;
    int depth_ = 0;
    static constexpr int kMaxDepth = 256;

    struct DepthGuard {
        int& depth_;
        explicit DepthGuard(int& depth) : depth_(depth) {
            if (++depth_ > kMaxDepth) {
                throw std::runtime_error("expression too complex: maximum parsing depth exceeded");
            }
        }
        ~DepthGuard() { --depth_; }
    };

    // ========================================================================
    // 解析规则实现
    // ========================================================================

    /**
     * @brief 解析加减法表达式
     *
     * expression -> term (('+' | '-') term)*
     */
    SymbolicExpression parse_expression() {
        DepthGuard guard(depth_);
        SymbolicExpression value = parse_term();
        while (true) {
            skip_spaces();
            if (match('+')) {
                // 加法在解析时简化，减少后续工作量
                value = SymbolicExpression(make_binary(NodeType::kAdd, value.simplify().node_, parse_term().simplify().node_));
            } else if (match('-')) {
                value = SymbolicExpression(make_binary(NodeType::kSubtract, value.simplify().node_, parse_term().simplify().node_));
            } else {
                return value;
            }
        }
    }

    /**
     * @brief 解析乘除法表达式
     *
     * term -> unary (('*' | '/') unary)*
     */
    SymbolicExpression parse_term() {
        SymbolicExpression value = parse_unary();
        while (true) {
            skip_spaces();
            if (match('*')) {
                value = SymbolicExpression(make_binary(NodeType::kMultiply, value.node_, parse_unary().node_));
            } else if (match('/')) {
                value = SymbolicExpression(make_binary(NodeType::kDivide, value.node_, parse_unary().node_));
            } else {
                return value;
            }
        }
    }

    /**
     * @brief 解析幂运算表达式
     *
     * power -> primary ('^' unary)?
     * 注意：幂运算是右结合的
     */
    SymbolicExpression parse_power() {
        SymbolicExpression value = parse_primary();
        skip_spaces();
        if (match('^')) {
            // 幂运算右结合：x^y^z = x^(y^z)
            return SymbolicExpression(make_binary(NodeType::kPower, value.node_, parse_unary().node_));
        }
        return value;
    }

    /**
     * @brief 解析一元表达式
     *
     * unary -> ('+' | '-') unary | power
     */
    SymbolicExpression parse_unary() {
        skip_spaces();
        if (match('+')) {
            // 正号可忽略
            return parse_unary();
        }
        if (match('-')) {
            // 负号创建取负节点
            return SymbolicExpression(make_unary(NodeType::kNegate, parse_unary().node_));
        }
        return parse_power();
    }

    /**
     * @brief 解析基本表达式
     *
     * primary -> '(' expression ')' | identifier '(' expression ')' | identifier | number
     */
    SymbolicExpression parse_primary() {
        skip_spaces();

        // 括号表达式
        if (match('(')) {
            SymbolicExpression value = parse_expression();
            skip_spaces();
            expect(')');
            return value;
        }

        // 标识符（变量或函数）
        if (peek_is_alpha()) {
            std::string identifier = parse_identifier();

            // 特殊常量
            if (identifier == "pi") {
                return SymbolicExpression(make_unary(NodeType::kPi, nullptr));
            }
            if (identifier == "e") {
                return SymbolicExpression(make_unary(NodeType::kE, nullptr));
            }
            if (identifier == "inf" || identifier == "infinity") {
                return SymbolicExpression::number(mymath::infinity());
            }
            if (identifier == "nan") {
                return SymbolicExpression::number(0.0 / 0.0);
            }

            skip_spaces();

            // 函数调用
            if (match('(')) {
                // 别名规范化
                if (identifier == "u" || identifier == "heaviside") {
                    identifier = "step";
                } else if (identifier == "impulse") {
                    identifier = "delta";
                }
                SymbolicExpression argument = parse_expression();
                skip_spaces();
                expect(')');
                return SymbolicExpression(make_unary(NodeType::kFunction, argument.node_, identifier));
            }

            // 变量
            return SymbolicExpression(make_variable(identifier));
        }

        // 数值
        return parse_number();
    }

    /**
     * @brief 解析数值
     *
     * 支持整数和小数，不支持科学计数法。
     */
    SymbolicExpression parse_number() {
        skip_spaces();
        const std::size_t start = pos_;
        bool has_digit = false;
        bool seen_dot = false;

        // 收集数字字符
        while (pos_ < source_.size()) {
            const char ch = source_[pos_];
            if (std::isdigit(static_cast<unsigned char>(ch))) {
                has_digit = true;
                ++pos_;
            } else if (ch == '.' && !seen_dot) {
                seen_dot = true;
                ++pos_;
            } else {
                break;
            }
        }

        if (!has_digit) {
            throw std::runtime_error("expected number");
        }

        // 手动解析数值（避免 strtod 的依赖和区域设置问题）
        double value = 0.0;
        std::size_t index = start;
        while (index < pos_ && source_[index] != '.') {
            value = value * 10.0 + static_cast<double>(source_[index] - '0');
            ++index;
        }
        if (index < pos_ && source_[index] == '.') {
            ++index;
            double place = 0.1;
            while (index < pos_) {
                value += static_cast<double>(source_[index] - '0') * place;
                place *= 0.1;
                ++index;
            }
        }

        // 解析科学计数法
        if (pos_ < source_.size() && (source_[pos_] == 'e' || source_[pos_] == 'E')) {
            std::size_t exp_pos = pos_ + 1;
            bool exp_negative = false;
            if (exp_pos < source_.size() && source_[exp_pos] == '+') {
                ++exp_pos;
            } else if (exp_pos < source_.size() && source_[exp_pos] == '-') {
                exp_negative = true;
                ++exp_pos;
            }
            
            if (exp_pos < source_.size() && std::isdigit(static_cast<unsigned char>(source_[exp_pos]))) {
                double exponent = 0.0;
                while (exp_pos < source_.size() && std::isdigit(static_cast<unsigned char>(source_[exp_pos]))) {
                    exponent = exponent * 10.0 + static_cast<double>(source_[exp_pos] - '0');
                    ++exp_pos;
                }
                pos_ = exp_pos;
                if (exp_negative) {
                    value /= mymath::pow(10.0, exponent);
                } else {
                    value *= mymath::pow(10.0, exponent);
                }
            }
        }

        return SymbolicExpression(make_number(value));
    }

    // ========================================================================
    // 辅助函数
    // ========================================================================

    /** @brief 解析标识符（字母开头，含字母数字下划线） */
    std::string parse_identifier() {
        const std::size_t start = pos_;
        while (pos_ < source_.size()) {
            const char ch = source_[pos_];
            if (std::isalnum(static_cast<unsigned char>(ch)) || ch == '_') {
                ++pos_;
            } else {
                break;
            }
        }
        return source_.substr(start, pos_ - start);
    }

    /** @brief 检查当前位置是否为字母 */
    bool peek_is_alpha() const {
        return pos_ < source_.size() &&
               std::isalpha(static_cast<unsigned char>(source_[pos_]));
    }

    /** @brief 尝试匹配指定字符，成功则前进 */
    bool match(char ch) {
        if (pos_ >= source_.size() || source_[pos_] != ch) {
            return false;
        }
        ++pos_;
        return true;
    }

    /** @brief 期望指定字符，不存在则抛出异常 */
    void expect(char ch) {
        if (!match(ch)) {
            throw std::runtime_error(std::string("expected '") + ch + "'");
        }
    }

    /** @brief 跳过空白字符 */
    void skip_spaces() {
        while (pos_ < source_.size() &&
               std::isspace(static_cast<unsigned char>(source_[pos_]))) {
            ++pos_;
        }
    }
};

bool expr_is_number(const SymbolicExpression& expression, double* value);

SymbolicExpression simplify_impl(const SymbolicExpression& expression);
SymbolicExpression simplify_once(const SymbolicExpression& expression);

SymbolicExpression substitute_impl(const SymbolicExpression& expression,
                                  const std::string& variable_name,
                                  const SymbolicExpression& replacement) {
    const auto& node = expression.node_;
    switch (node->type) {
        case NodeType::kNumber:
        case NodeType::kPi:
        case NodeType::kE:
            return expression;
        case NodeType::kVariable:
            if (node->text == variable_name) {
                return replacement;
            }
            return expression;
        case NodeType::kNegate:
            return SymbolicExpression(
                       make_unary(NodeType::kNegate,
                                  substitute_impl(SymbolicExpression(node->left),
                                                  variable_name,
                                                  replacement)
                                      .node_))
                .simplify();
        case NodeType::kFunction:
            return SymbolicExpression(
                       make_unary(NodeType::kFunction,
                                  substitute_impl(SymbolicExpression(node->left),
                                                  variable_name,
                                                  replacement)
                                      .node_,
                                  node->text))
                .simplify();
        case NodeType::kAdd:
        case NodeType::kSubtract:
        case NodeType::kMultiply:
        case NodeType::kDivide:
        case NodeType::kPower:
            return SymbolicExpression(
                       make_binary(node->type,
                                   substitute_impl(SymbolicExpression(node->left),
                                                   variable_name,
                                                   replacement)
                                       .node_,
                                   substitute_impl(SymbolicExpression(node->right),
                                                   variable_name,
                                                   replacement)
                                       .node_))
                .simplify();
    }
    throw std::runtime_error("unsupported symbolic substitution");
}

bool try_evaluate_numeric_node(const std::shared_ptr<SymbolicExpression::Node>& node,
                               double* value) {
    switch (node->type) {
        case NodeType::kNumber:
            *value = node->number_value;
            return true;
        case NodeType::kPi:
        case NodeType::kE:
            return false;
        case NodeType::kVariable:
            if (node->text == "1 / 2") {
                *value = 0.5;
                return true;
            }
            return false;
        case NodeType::kNegate: {
            double operand = 0.0;
            if (!try_evaluate_numeric_node(node->left, &operand)) {
                return false;
            }
            *value = -operand;
            return true;
        }
        case NodeType::kAdd:
        case NodeType::kSubtract:
        case NodeType::kMultiply:
        case NodeType::kDivide:
        case NodeType::kPower: {
            double left = 0.0;
            double right = 0.0;
            if (!try_evaluate_numeric_node(node->left, &left) ||
                !try_evaluate_numeric_node(node->right, &right)) {
                return false;
            }

            switch (node->type) {
                case NodeType::kAdd:
                    *value = left + right;
                    return true;
                case NodeType::kSubtract:
                    *value = left - right;
                    return true;
                case NodeType::kMultiply:
                    *value = left * right;
                    return true;
                case NodeType::kDivide:
                    *value = left / right;
                    return true;
                case NodeType::kPower:
                    *value = mymath::pow(left, right);
                    return true;
                case NodeType::kNumber:
                case NodeType::kPi:
                case NodeType::kE:
                case NodeType::kVariable:
                case NodeType::kNegate:
                case NodeType::kFunction:
                    break;
            }
            return false;
        }
        case NodeType::kFunction: {
            double argument = 0.0;
            if (!try_evaluate_numeric_node(node->left, &argument)) {
                return false;
            }
            if (node->text == "asin") {
                *value = mymath::asin(argument);
                return true;
            }
            if (node->text == "acos") {
                *value = mymath::acos(argument);
                return true;
            }
            if (node->text == "atan") {
                *value = mymath::atan(argument);
                return true;
            }
            if (node->text == "sin") {
                *value = mymath::sin(argument);
                return true;
            }
            if (node->text == "cos") {
                *value = mymath::cos(argument);
                return true;
            }
            if (node->text == "tan") {
                *value = mymath::tan(argument);
                return true;
            }
            if (node->text == "exp") {
                *value = mymath::exp(argument);
                return true;
            }
            if (node->text == "sinh") {
                *value = mymath::sinh(argument);
                return true;
            }
            if (node->text == "cosh") {
                *value = mymath::cosh(argument);
                return true;
            }
            if (node->text == "tanh") {
                *value = mymath::tanh(argument);
                return true;
            }
            if (node->text == "ln") {
                *value = mymath::ln(argument);
                return true;
            }
            if (node->text == "sqrt") {
                double root = mymath::sqrt(argument);
                // Only evaluate to number if it's a perfect square
                if (mymath::is_near_zero(root * root - argument, 1e-12) && mymath::is_integer(root, 1e-10)) {
                    *value = root;
                    return true;
                }
                return false;
            }
            if (node->text == "erf") {
                *value = mymath::erf(argument);
                return true;
            }
            if (node->text == "erfc") {
                *value = mymath::erfc(argument);
                return true;
            }
            if (node->text == "gamma") {
                *value = mymath::gamma(argument);
                return true;
            }
            if (node->text == "abs") {
                *value = mymath::abs(argument);
                return true;
            }
            if (node->text == "floor") {
                *value = static_cast<double>(
                    static_cast<long long>(argument < 0.0 &&
                                                   static_cast<double>(static_cast<long long>(argument)) != argument
                                               ? argument - 1.0
                                               : argument));
                return true;
            }
            if (node->text == "ceil") {
                long long truncated = static_cast<long long>(argument);
                if (argument > 0.0 && static_cast<double>(truncated) != argument) {
                    ++truncated;
                }
                *value = static_cast<double>(truncated);
                return true;
            }
            if (node->text == "cbrt") {
                *value = mymath::cbrt(argument);
                return true;
            }
            if (node->text == "sign") {
                if (mymath::is_near_zero(argument, kFormatEps)) {
                    *value = 0.0;
                } else {
                    *value = argument > 0.0 ? 1.0 : -1.0;
                }
                return true;
            }
            if (node->text == "step") {
                *value = argument >= 0.0 ? 1.0 : 0.0;
                return true;
            }
            if (node->text == "delta") {
                *value = mymath::is_near_zero(argument, kFormatEps) ? 1.0 : 0.0;
                return true;
            }
            return false;
        }
    }
    return false;
}

bool expr_is_variable(const SymbolicExpression& expression, const std::string& name) {
    if (name == "pi") return expression.node_->type == NodeType::kPi;
    if (name == "e") return expression.node_->type == NodeType::kE;
    return expression.node_->type == NodeType::kVariable && expression.node_->text == name;
}

}  // namespace symbolic_expression_internal

using namespace symbolic_expression_internal;

SymbolicExpression::SymbolicExpression() : node_(make_number(0.0)) {}

SymbolicExpression::SymbolicExpression(std::shared_ptr<Node> node)
    : node_(std::move(node)) {}

SymbolicExpression SymbolicExpression::parse(const std::string& text) {
    Parser parser(text);
    return parser.parse().simplify();
}

SymbolicExpression SymbolicExpression::number(double value) {
    return SymbolicExpression(make_number(value));
}

void SymbolicExpression::set_display_precision(int precision) {
    mutable_display_precision() = clamp_display_precision(precision);
}

SymbolicExpression operator+(const SymbolicExpression& lhs, const SymbolicExpression& rhs) {
    return make_add(lhs, rhs);
}

SymbolicExpression operator-(const SymbolicExpression& lhs, const SymbolicExpression& rhs) {
    return make_subtract(lhs, rhs);
}

SymbolicExpression operator*(const SymbolicExpression& lhs, const SymbolicExpression& rhs) {
    return make_multiply(lhs, rhs);
}

SymbolicExpression operator/(const SymbolicExpression& lhs, const SymbolicExpression& rhs) {
    return make_divide(lhs, rhs);
}

SymbolicExpression operator^(const SymbolicExpression& lhs, const SymbolicExpression& rhs) {
    return make_power(lhs, rhs);
}

SymbolicExpression operator-(const SymbolicExpression& expr) {
    return make_negate(expr);
}

SymbolicExpression SymbolicExpression::variable(const std::string& name) {
    return SymbolicExpression(make_variable(name));
}

std::string SymbolicExpression::to_string() const {
    return to_string_impl(simplify().node_, 0);
}

bool SymbolicExpression::is_constant(const std::string& variable_name) const {
    switch (node_->type) {
        case NodeType::kNumber:
        case NodeType::kPi:
        case NodeType::kE:
            return true;
        case NodeType::kVariable:
            return node_->text != variable_name;
        case NodeType::kNegate:
        case NodeType::kFunction:
            return SymbolicExpression(node_->left).is_constant(variable_name);
        case NodeType::kAdd:
        case NodeType::kSubtract:
        case NodeType::kMultiply:
        case NodeType::kDivide:
        case NodeType::kPower:
            return SymbolicExpression(node_->left).is_constant(variable_name) &&
                   SymbolicExpression(node_->right).is_constant(variable_name);
    }
    return false;
}

bool SymbolicExpression::is_number(double* value) const {
    double numeric = 0.0;
    if (!try_evaluate_numeric_node(node_, &numeric)) {
        return false;
    }
    if (value != nullptr) {
        *value = numeric;
    }
    return true;
}

bool SymbolicExpression::is_variable_named(const std::string& variable_name) const {
    return node_->type == NodeType::kVariable && node_->text == variable_name;
}

bool SymbolicExpression::polynomial_coefficients(
    const std::string& variable_name,
    std::vector<double>* coefficients) const {
    const SymbolicExpression simplified = simplify();
    return polynomial_coefficients_from_simplified(simplified, variable_name, coefficients);
}

std::vector<std::string> SymbolicExpression::identifier_variables() const {
    std::vector<std::string> names;
    collect_identifier_variables(simplify(), &names);
    std::sort(names.begin(), names.end());
    names.erase(std::unique(names.begin(), names.end()), names.end());
    return names;
}

namespace {
void collect_subexpression_counts(const std::shared_ptr<SymbolicExpression::Node>& node,
                                  std::unordered_map<std::string, std::pair<std::shared_ptr<SymbolicExpression::Node>, int>>& counts) {
    if (!node) return;
    
    // 忽略叶子节点（数字、变量、常数），因为它们作为 CSE 提取没有意义
    if (node->type != NodeType::kNumber && node->type != NodeType::kVariable &&
        node->type != NodeType::kPi && node->type != NodeType::kE) {
        const std::string key = node_structural_key(node);
        auto& entry = counts[key];
        if (entry.second == 0) {
            entry.first = node;
        }
        entry.second++;
    }

    collect_subexpression_counts(node->left, counts);
    collect_subexpression_counts(node->right, counts);
}
} // namespace

std::vector<std::pair<SymbolicExpression, int>> SymbolicExpression::common_subexpressions() const {
    std::unordered_map<std::string, std::pair<std::shared_ptr<Node>, int>> counts;
    collect_subexpression_counts(node_, counts);

    std::vector<std::pair<SymbolicExpression, int>> result;
    for (auto& [key, entry] : counts) {
        if (entry.second > 1) {
            result.push_back({SymbolicExpression(entry.first), entry.second});
        }
    }

    // 按出现次数降序排列，次数相同时按长度（通过字符串表示的长度估计）降序
    std::sort(result.begin(), result.end(), [](const auto& a, const auto& b) {
        if (a.second != b.second) return a.second > b.second;
        return a.first.to_string().size() > b.first.to_string().size();
    });

    return result;
}

SymbolicExpression SymbolicExpression::simplify() const {
    static constexpr std::size_t kMaxSimplifyCacheSize = 4096;
    static thread_local SymbolicExpressionLruCache cache(kMaxSimplifyCacheSize);

    const std::string key = node_structural_key(node_);
    SymbolicExpression cached;
    if (cache.get(key, &cached)) {
        return cached;
    }

    SymbolicExpression simplified = simplify_impl(*this);
    cache.put(key, simplified);
    return simplified;
}

SymbolicExpression SymbolicExpression::expand() const {
    return expand_impl(*this);
}

SymbolicExpression SymbolicExpression::substitute(
    const std::string& variable_name,
    const SymbolicExpression& replacement) const {
    if (!is_identifier_variable_name(variable_name) ||
        variable_name == "pi" || variable_name == "e" || variable_name == "i") {
        throw std::runtime_error(
            "symbolic substitution variable must be a non-reserved identifier");
    }
    return substitute_impl(*this, variable_name, replacement).simplify();
}
