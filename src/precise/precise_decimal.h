// ============================================================================
// 精确小数类型
// ============================================================================
//
// 本文件定义了 PreciseDecimal 类型，用于精确表示小数数值。
// 与浮点数不同，PreciseDecimal 使用字符串存储数字，完全避免浮点误差。
//
// 主要特性：
// - 支持任意精度的小数表示和运算
// - 内部使用基数为 10^9 的大整数数组存储有效数字
// - 提供完整的算术运算和数学函数支持
// - 适用于需要精确计算的场景，如货币计算、科学计算等
//
// 使用示例：
//   PreciseDecimal a("0.1");
//   PreciseDecimal b("0.2");
//   PreciseDecimal c = a + b;  // 结果精确为 0.3
// ============================================================================

#ifndef TYPES_PRECISE_DECIMAL_H
#define TYPES_PRECISE_DECIMAL_H

#include <string>
#include <vector>
#include <cstdint>
#include <iostream>
#include "math/types/float128.h"

/**
 * @struct PrecisionContext
 * @brief 全局高精度计算配置上下文
 *
 * 用于管理全局的精度设置，如默认的小数位数。
 * 通过 get_default_scale() 和 set_default_scale() 可以查询和修改默认精度。
 */
struct PrecisionContext {
    /**
     * @brief 获取默认精度（小数点后的位数）
     * @return 当前默认精度值
     */
    static int get_default_scale();

    /**
     * @brief 设置默认精度（小数点后的位数）
     * @param scale 新的默认精度值
     */
    static void set_default_scale(int scale);
};

/**
 * @struct ScopedPrecision
 * @brief 局部高精度上下文（RAII 机制），用于自动管理计算时的保护位
 */
struct ScopedPrecision {
    int old_scale;
    ScopedPrecision(int extra_scale) {
        old_scale = PrecisionContext::get_default_scale();
        PrecisionContext::set_default_scale(old_scale + extra_scale);
    }
    ~ScopedPrecision() {
        PrecisionContext::set_default_scale(old_scale);
    }
};

/**
 * @struct BigIntData
 * @brief 具有小对象优化（SSO）的大整数存储结构
 *
 * 当数据较小时直接存储在栈上，避免堆分配。
 */
struct BigIntData {
    static constexpr size_t SSO_CAP = 16;
    uint32_t* ptr;
    uint32_t size_;
    uint32_t capacity_;
    uint32_t sso[SSO_CAP];

    BigIntData();
    BigIntData(size_t n, uint32_t val = 0);
    BigIntData(const BigIntData& other);
    BigIntData(BigIntData&& other) noexcept;
    BigIntData(std::initializer_list<uint32_t> list);
    BigIntData(const uint32_t* first, const uint32_t* last);
    ~BigIntData();

    BigIntData& operator=(const BigIntData& other);
    BigIntData& operator=(BigIntData&& other) noexcept;

    size_t size() const { return size_; }
    bool empty() const { return size_ == 0; }
    uint32_t& operator[](size_t i) { return ptr[i]; }
    const uint32_t& operator[](size_t i) const { return ptr[i]; }
    uint32_t& back() { return ptr[size_ - 1]; }
    const uint32_t& back() const { return ptr[size_ - 1]; }

    void push_back(uint32_t val);
    void pop_back() { if (size_ > 0) size_--; }
    void clear() { size_ = 0; }
    void reserve(size_t n);
    void resize(size_t n, uint32_t val = 0);

    uint32_t* begin() { return ptr; }
    const uint32_t* begin() const { return ptr; }
    uint32_t* end() { return ptr + size_; }
    const uint32_t* end() const { return ptr + size_; }
    uint32_t* data() { return ptr; }
    const uint32_t* data() const { return ptr; }

    void erase(const uint32_t* first, const uint32_t* last);
    void insert(uint32_t* pos, const uint32_t* first, const uint32_t* last);
    void insert(uint32_t* pos, size_t count, uint32_t val);

    /**
     * @brief 获取当前容量
     */
    size_t capacity() const { return capacity_; }
};

/**
 * @class BigIntArena
 * @brief 简单的内存池，用于 NTT 等大规模运算
 *
 * 预分配一大块连续内存，通过偏移量管理内部的小型数组需求。
 * 避免频繁的 new/delete 调用。
 */
class BigIntArena {
public:
    static constexpr size_t DEFAULT_SIZE = 1024 * 1024;  // 4MB (1M uint32_t)

    explicit BigIntArena(size_t initial_size = DEFAULT_SIZE)
        : buffer_(new uint32_t[initial_size]), capacity_(initial_size), offset_(0) {}

    ~BigIntArena() { delete[] buffer_; }

    // 禁止拷贝
    BigIntArena(const BigIntArena&) = delete;
    BigIntArena& operator=(const BigIntArena&) = delete;

    // 允许移动
    BigIntArena(BigIntArena&& other) noexcept
        : buffer_(other.buffer_), capacity_(other.capacity_), offset_(other.offset_) {
        other.buffer_ = nullptr;
        other.capacity_ = 0;
        other.offset_ = 0;
    }

    BigIntArena& operator=(BigIntArena&& other) noexcept {
        if (this != &other) {
            delete[] buffer_;
            buffer_ = other.buffer_;
            capacity_ = other.capacity_;
            offset_ = other.offset_;
            other.buffer_ = nullptr;
            other.capacity_ = 0;
            other.offset_ = 0;
        }
        return *this;
    }

    /**
     * @brief 分配 n 个 uint32_t
     * @return 指向分配内存的指针
     */
    uint32_t* allocate(size_t n) {
        if (offset_ + n > capacity_) {
            // 扩容：新大小为 max(2*capacity_, offset_ + n)
            size_t new_capacity = std::max(capacity_ * 2, offset_ + n);
            uint32_t* new_buffer = new uint32_t[new_capacity];
            std::copy(buffer_, buffer_ + offset_, new_buffer);
            delete[] buffer_;
            buffer_ = new_buffer;
            capacity_ = new_capacity;
        }
        uint32_t* result = buffer_ + offset_;
        offset_ += n;
        return result;
    }

    /**
     * @brief 分配并初始化为指定值
     */
    uint32_t* allocate(size_t n, uint32_t val) {
        uint32_t* result = allocate(n);
        std::fill(result, result + n, val);
        return result;
    }

    /**
     * @brief 重置内存池（保留内存，只重置偏移量）
     */
    void reset() { offset_ = 0; }

    /**
     * @brief 获取当前使用量
     */
    size_t used() const { return offset_; }

    /**
     * @brief 获取总容量
     */
    size_t capacity() const { return capacity_; }

private:
    uint32_t* buffer_;
    size_t capacity_;
    size_t offset_;
};

/**
 * @struct BinaryBigInt
 * @brief 二进制基数的大整数表示（基数 2^32）
 *
 * 用于内部高效运算。使用 2^32 作为基数时：
 * - 乘法可以通过位移和掩码快速完成
 * - 模运算和除法通过位操作完成
 * - 仅在输入/输出时进行进制转换
 */
struct BinaryBigInt {
    static constexpr uint64_t kBase = 0x100000000ULL;  // 2^32
    static constexpr uint32_t kMask = 0xFFFFFFFF;       // 32位掩码

    BigIntData limbs;  ///< 小端序存储，每个 limb 是 32 位

    BinaryBigInt() = default;
    explicit BinaryBigInt(uint64_t val);
    explicit BinaryBigInt(const BigIntData& decimal_data);  // 从十进制基数转换
    explicit BinaryBigInt(BigIntData&& decimal_data);       // 从十进制基数转换（移动）

    /**
     * @brief 转换回十进制基数表示
     */
    BigIntData to_decimal() const;

    /**
     * @brief 检查是否为零
     */
    bool is_zero() const { return limbs.empty() || (limbs.size() == 1 && limbs[0] == 0); }

    /**
     * @brief 获取有效位数（去除前导零）
     */
    size_t effective_size() const;

    /**
     * @brief 规范化（去除前导零）
     */
    void normalize();

    // ==================== 算术运算 ====================

    /**
     * @brief 二进制大整数加法
     */
    static BinaryBigInt add(const BinaryBigInt& a, const BinaryBigInt& b);

    /**
     * @brief 二进制大整数减法（要求 a >= b）
     */
    static BinaryBigInt sub(const BinaryBigInt& a, const BinaryBigInt& b);

    /**
     * @brief 二进制大整数乘法（使用快速算法）
     */
    static BinaryBigInt mul(const BinaryBigInt& a, const BinaryBigInt& b);

    /**
     * @brief 二进制大整数乘以小整数
     */
    static BinaryBigInt mul_uint32(const BinaryBigInt& a, uint32_t b);

    /**
     * @brief 二进制大整数除法
     */
    static void div(const BinaryBigInt& num, const BinaryBigInt& den,
                    BinaryBigInt* quotient, BinaryBigInt* remainder);

    /**
     * @brief 比较两个二进制大整数
     * @return -1 表示 a < b，0 表示相等，1 表示 a > b
     */
    static int compare(const BinaryBigInt& a, const BinaryBigInt& b);

    /**
     * @brief 左移（乘以 2^k）
     */
    BinaryBigInt operator<<(int bits) const;

    /**
     * @brief 右移（除以 2^k）
     */
    BinaryBigInt operator>>(int bits) const;
};

/**
 * @struct PreciseDecimal
 * @brief 精确小数表示结构体
 *
 * 使用大整数数组存储数字，实现任意精度的十进制小数表示。
 *
 * 内部表示说明：
 * - data: 存储以 10^9 为基数的有效数字数组（小端序）
 * - scale: 表示小数点后的总位数
 * - negative: 表示是否为负数
 *
 * 示例：
 * - 123.45 → data={12345}, scale=2
 * - 0.000001 → data={1}, scale=6
 * - 1000000000000 → data={0, 1}, scale=0 （因为 10^9 为基数）
 */
struct PreciseDecimal {
    BigIntData data;                  ///< 基数为 10^9 的有效数字数组（小端序）
    int scale = 0;                    ///< 小数点后的位数
    bool negative = false;            ///< 是否为负数
    bool is_inf = false;              ///< 是否为无穷大
    bool is_nan = false;              ///< 是否为 NaN

    // ==================== 构造函数 ====================

    /** @brief 默认构造函数，初始化为零 */
    PreciseDecimal() = default;

    /**
     * @brief 从 long long 整数构造
     * @param value 整数值
     */
    PreciseDecimal(long long value);

    /**
     * @brief 从 unsigned long long 整数构造
     * @param value 无符号整数值
     */
    PreciseDecimal(unsigned long long value) : PreciseDecimal(static_cast<long long>(value)) {}

    /**
     * @brief 从 int 整数构造
     * @param value 整数值
     */
    PreciseDecimal(int value) : PreciseDecimal(static_cast<long long>(value)) {}

    /**
     * @brief 从 long double 浮点数构造（可能有精度损失）
     * @param value 浮点数值
     */
    PreciseDecimal(long double value);

    /**
     * @brief 从 double 浮点数构造（可能有精度损失）
     * @param value 浮点数值
     */
    PreciseDecimal(double value) : PreciseDecimal(static_cast<long double>(value)) {}

    /**
     * @brief 从 float 浮点数构造（可能有精度损失）
     * @param value 浮点数值
     */
    PreciseDecimal(float value) : PreciseDecimal(static_cast<long double>(value)) {}

    /**
     * @brief 从字符串构造
     * @param token 数值字符串，如 "123.45", "1e-6", "-0.5"
     */
    explicit PreciseDecimal(const std::string& token);

    // ==================== 类型转换运算符 ====================

    /** @brief 显式转换为 long double */
    explicit operator long double() const { return to_double(); }

    /** @brief 显式转换为 double */
    explicit operator double() const { return static_cast<double>(to_double()); }

    /** @brief 显式转换为 float */
    explicit operator float() const { return static_cast<float>(to_double()); }

    /** @brief 显式转换为 int */
    explicit operator int() const { return static_cast<int>(to_double()); }

    /** @brief 显式转换为 long long */
    explicit operator long long() const { return static_cast<long long>(to_double()); }

    // ==================== 基本操作 ====================

    /**
     * @brief 规范化表示
     *
     * 去除前导零和末尾零，确保表示的唯一性。
     * 例如：data={100}, scale=2 规范化为 data={1}, scale=0
     */
    void normalize();

    /**
     * @brief 检查值是否为零
     * @return 如果值为零返回 true
     */
    bool is_zero() const;

    /**
     * @brief 转换为字符串表示
     * @return 数值的字符串形式，如 "123.45", "-0.001"
     */
    std::string to_string() const;

    /**
     * @brief 转换为 long double（可能有精度损失）
     * @return 近似的浮点数值
     */
    long double to_double() const;

    /**
     * @brief 转换为 long double（to_double 的别名）
     * @return 近似的浮点数值
     */
    long double to_long_double() const { return to_double(); }

    // ==================== 算术运算符（成员函数） ====================

    /** @brief 取负运算 */
    PreciseDecimal operator-() const;

    /** @brief 加法赋值运算 */
    PreciseDecimal& operator+=(const PreciseDecimal& rhs);

    /** @brief 减法赋值运算 */
    PreciseDecimal& operator-=(const PreciseDecimal& rhs);

    /** @brief 乘法赋值运算 */
    PreciseDecimal& operator*=(const PreciseDecimal& rhs);

    /** @brief 除法赋值运算 */
    PreciseDecimal& operator/=(const PreciseDecimal& rhs);

    // ==================== 混合类型运算辅助 ====================

    /** @brief 与 long double 的加法赋值 */
    PreciseDecimal& operator+=(long double rhs) { return *this += PreciseDecimal(rhs); }

    /** @brief 与 long double 的减法赋值 */
    PreciseDecimal& operator-=(long double rhs) { return *this -= PreciseDecimal(rhs); }

    /** @brief 与 long double 的乘法赋值 */
    PreciseDecimal& operator*=(long double rhs) { return *this *= PreciseDecimal(rhs); }

    /** @brief 与 long double 的除法赋值 */
    PreciseDecimal& operator/=(long double rhs) { return *this /= PreciseDecimal(rhs); }

    // ==================== 比较运算符 ====================

    bool operator==(const PreciseDecimal& rhs) const;
    bool operator!=(const PreciseDecimal& rhs) const { return !(*this == rhs); }
    bool operator<(const PreciseDecimal& rhs) const;
    bool operator>(const PreciseDecimal& rhs) const { return rhs < *this; }
    bool operator<=(const PreciseDecimal& rhs) const { return !(*this > rhs); }
    bool operator>=(const PreciseDecimal& rhs) const { return !(*this < rhs); }

    // ==================== 静态工厂方法 ====================

    /**
     * @brief 从原始数字字符串构造
     * @param raw_digits 数字字符串（不含小数点和符号）
     * @param raw_scale 小数点后的位数
     * @param is_negative 是否为负数
     * @return 构造的 PreciseDecimal 对象
     */
    static PreciseDecimal from_digits(std::string raw_digits,
                                      int raw_scale,
                                      bool is_negative);

    /**
     * @brief 从整数字符串构造
     * @param integer_text 整数字符串
     * @param is_negative 是否为负数
     * @return 构造的 PreciseDecimal 对象
     */
    static PreciseDecimal from_integer_string(const std::string& integer_text,
                                              bool is_negative);

    /**
     * @brief 从小数字面量字符串构造
     * @param token 小数字面量，如 "123.45", "1e-6", "-0.5"
     * @return 构造的 PreciseDecimal 对象
     */
    static PreciseDecimal from_decimal_literal(const std::string& token);

    /**
     * @brief 创建正无穷大
     * @return 表示正无穷大的 PreciseDecimal 对象
     */
    static PreciseDecimal infinity() {
        PreciseDecimal result;
        result.is_inf = true;
        result.negative = false;
        return result;
    }

    /**
     * @brief 创建负无穷大
     * @return 表示负无穷大的 PreciseDecimal 对象
     */
    static PreciseDecimal neg_infinity() {
        PreciseDecimal result;
        result.is_inf = true;
        result.negative = true;
        return result;
    }

    /**
     * @brief 创建 NaN
     * @return 表示 NaN 的 PreciseDecimal 对象
     */
    static PreciseDecimal nan() {
        PreciseDecimal result;
        result.is_nan = true;
        return result;
    }

    /**
     * @brief 检查是否为无穷大
     * @return 如果是无穷大返回 true
     */
    bool is_infinity() const { return is_inf; }

    /**
     * @brief 检查是否为 NaN
     * @return 如果是 NaN 返回 true
     */
    bool is_nan_value() const { return is_nan; }
};

// ============================================================================
// 非成员二元运算符（支持混合类型运算）
// ============================================================================

// 加法运算符
PreciseDecimal operator+(PreciseDecimal lhs, const PreciseDecimal& rhs);
PreciseDecimal operator+(PreciseDecimal lhs, long double rhs);
PreciseDecimal operator+(long double lhs, const PreciseDecimal& rhs);
PreciseDecimal operator+(PreciseDecimal lhs, double rhs);
PreciseDecimal operator+(double lhs, const PreciseDecimal& rhs);
PreciseDecimal operator+(PreciseDecimal lhs, int rhs);
PreciseDecimal operator+(int lhs, const PreciseDecimal& rhs);

// 减法运算符
PreciseDecimal operator-(PreciseDecimal lhs, const PreciseDecimal& rhs);
PreciseDecimal operator-(PreciseDecimal lhs, long double rhs);
PreciseDecimal operator-(long double lhs, const PreciseDecimal& rhs);
PreciseDecimal operator-(PreciseDecimal lhs, double rhs);
PreciseDecimal operator-(double lhs, const PreciseDecimal& rhs);
PreciseDecimal operator-(PreciseDecimal lhs, int rhs);
PreciseDecimal operator-(int lhs, const PreciseDecimal& rhs);

// 乘法运算符
PreciseDecimal operator*(PreciseDecimal lhs, const PreciseDecimal& rhs);
PreciseDecimal operator*(PreciseDecimal lhs, long double rhs);
PreciseDecimal operator*(long double lhs, const PreciseDecimal& rhs);
PreciseDecimal operator*(PreciseDecimal lhs, double rhs);
PreciseDecimal operator*(double lhs, const PreciseDecimal& rhs);
PreciseDecimal operator*(PreciseDecimal lhs, int rhs);
PreciseDecimal operator*(int lhs, const PreciseDecimal& rhs);

// 除法运算符
PreciseDecimal operator/(PreciseDecimal lhs, const PreciseDecimal& rhs);
PreciseDecimal operator/(PreciseDecimal lhs, long double rhs);
PreciseDecimal operator/(long double lhs, const PreciseDecimal& rhs);
PreciseDecimal operator/(PreciseDecimal lhs, double rhs);
PreciseDecimal operator/(double lhs, const PreciseDecimal& rhs);
PreciseDecimal operator/(PreciseDecimal lhs, int rhs);
PreciseDecimal operator/(int lhs, const PreciseDecimal& rhs);

// 混合类型比较运算符
bool operator==(const PreciseDecimal& lhs, long double rhs);
bool operator==(long double lhs, const PreciseDecimal& rhs);
bool operator!=(const PreciseDecimal& lhs, long double rhs);
bool operator!=(long double lhs, const PreciseDecimal& rhs);
bool operator<(const PreciseDecimal& lhs, long double rhs);
bool operator<(long double lhs, const PreciseDecimal& rhs);
bool operator>(const PreciseDecimal& lhs, long double rhs);
bool operator>(long double lhs, const PreciseDecimal& rhs);
bool operator<=(const PreciseDecimal& lhs, long double rhs);
bool operator<=(long double lhs, const PreciseDecimal& rhs);
bool operator>=(const PreciseDecimal& lhs, long double rhs);
bool operator>=(long double lhs, const PreciseDecimal& rhs);

bool operator==(const PreciseDecimal& lhs, double rhs);
bool operator==(double lhs, const PreciseDecimal& rhs);
bool operator!=(const PreciseDecimal& lhs, double rhs);
bool operator!=(double lhs, const PreciseDecimal& rhs);
bool operator<(const PreciseDecimal& lhs, double rhs);
bool operator<(double lhs, const PreciseDecimal& rhs);
bool operator>(const PreciseDecimal& lhs, double rhs);
bool operator>(double lhs, const PreciseDecimal& rhs);
bool operator<=(const PreciseDecimal& lhs, double rhs);
bool operator<=(double lhs, const PreciseDecimal& rhs);
bool operator>=(const PreciseDecimal& lhs, double rhs);
bool operator>=(double lhs, const PreciseDecimal& rhs);

bool operator==(const PreciseDecimal& lhs, int rhs);
bool operator==(int lhs, const PreciseDecimal& rhs);
bool operator!=(const PreciseDecimal& lhs, int rhs);
bool operator!=(int lhs, const PreciseDecimal& rhs);
bool operator<(const PreciseDecimal& lhs, int rhs);
bool operator<(int lhs, const PreciseDecimal& rhs);
bool operator>(const PreciseDecimal& lhs, int rhs);
bool operator>(int lhs, const PreciseDecimal& rhs);
bool operator<=(const PreciseDecimal& lhs, int rhs);
bool operator<=(int lhs, const PreciseDecimal& rhs);
bool operator>=(const PreciseDecimal& lhs, int rhs);
bool operator>=(int lhs, const PreciseDecimal& rhs);

// 流输出运算符
std::ostream& operator<<(std::ostream& os, const PreciseDecimal& val);

// 无穷大函数
inline PreciseDecimal infinity() { return PreciseDecimal::infinity(); }
inline PreciseDecimal neg_infinity() { return PreciseDecimal::neg_infinity(); }

// 字符串转换函数
inline PreciseDecimal from_string(const std::string& s) { return PreciseDecimal(s); }

// ============================================================================
// 高精度数学函数命名空间
// ============================================================================

/**
 * @namespace precise
 * @brief 高精度数学函数集合
 *
 * 提供一系列高精度数学函数，包括基本运算、三角函数、
 * 指数对数函数等，所有函数都保证计算精度。
 */
namespace precise {

// ==================== 基本数学函数 ====================

/**
 * @brief 绝对值函数
 * @param val 输入值
 * @return 输入值的绝对值
 */
PreciseDecimal abs(const PreciseDecimal& val);

/**
 * @brief 平方根函数
 * @param val 输入值（必须非负）
 * @return 平方根值
 * @throws PreciseDecimalUnsupported 如果输入为负数
 */
PreciseDecimal sqrt(const PreciseDecimal& val);

/**
 * @brief 幂运算（整数指数）
 * @param base 底数
 * @param exp 指数
 * @return base^exp
 */
PreciseDecimal pow(const PreciseDecimal& base, long long exp);

/**
 * @brief 幂运算（任意指数）
 * @param base 底数
 * @param exp 指数
 * @return base^exp
 */
PreciseDecimal pow(const PreciseDecimal& base, const PreciseDecimal& exp);

/**
 * @brief 向下取整
 * @param val 输入值
 * @return 不大于输入值的最大整数
 */
PreciseDecimal floor(const PreciseDecimal& val);

/**
 * @brief 向上取整
 * @param val 输入值
 * @return 不小于输入值的最小整数
 */
PreciseDecimal ceil(const PreciseDecimal& val);

/**
 * @brief 四舍五入
 * @param val 输入值
 * @return 四舍五入后的整数
 */
PreciseDecimal round(const PreciseDecimal& val);

// ==================== 数学常量 ====================

/**
 * @brief 圆周率 π
 * @return 高精度的 π 值
 */
PreciseDecimal pi();

/**
 * @brief 2π 常量
 * @return 高精度的 2π 值
 */
PreciseDecimal two_pi();

/**
 * @brief π/2 常量
 * @return 高精度的 π/2 值
 */
PreciseDecimal half_pi();

/**
 * @brief 自然对数的底 e
 * @return 高精度的 e 值
 */
PreciseDecimal e();

/**
 * @brief ln(2) 常量
 * @return 高精度的 ln(2) 值
 */
PreciseDecimal ln2();

/**
 * @brief log2(e) 常量
 * @return 高精度的 log2(e) 值
 */
PreciseDecimal log2e();

// ==================== 指数与对数函数 ====================

/**
 * @brief 自然指数函数
 * @param x 指数值
 * @return e^x
 */
PreciseDecimal exp(const PreciseDecimal& x);

/**
 * @brief 2 的指数函数（优化版）
 * @param x 指数值
 * @return 2^x
 */
PreciseDecimal exp2(const PreciseDecimal& x);

/**
 * @brief 自然对数函数
 * @param x 输入值（必须为正数）
 * @return ln(x)
 * @throws std::domain_error 如果 x 非正
 */
PreciseDecimal ln(const PreciseDecimal& x);

/**
 * @brief 以 2 为底的对数函数
 * @param x 输入值（必须为正数）
 * @return log2(x)
 * @throws std::domain_error 如果 x 非正
 */
PreciseDecimal log2(const PreciseDecimal& x);

/**
 * @brief 常用对数函数（以 10 为底）
 * @param x 输入值（必须为正数）
 * @return log10(x)
 */
PreciseDecimal log10(const PreciseDecimal& x);

/**
 * @brief 快速幂运算（优化版）
 * @param base 底数
 * @param exp 指数
 * @return base^exp
 */
PreciseDecimal pow_fast(const PreciseDecimal& base, const PreciseDecimal& exp);

// ==================== 三角函数 ====================

/**
 * @brief 正弦函数
 * @param x 弧度值
 * @return sin(x)
 */
PreciseDecimal sin(const PreciseDecimal& x);

/**
 * @brief 余弦函数
 * @param x 弧度值
 * @return cos(x)
 */
PreciseDecimal cos(const PreciseDecimal& x);

/**
 * @brief 正切函数
 * @param x 弧度值
 * @return tan(x)
 * @throws std::domain_error 如果 x 接近 π/2 + kπ
 */
PreciseDecimal tan(const PreciseDecimal& x);

/**
 * @brief 反正弦函数
 * @param x 输入值（范围 [-1, 1]）
 * @return asin(x)（弧度）
 * @throws std::domain_error 如果 x 超出范围
 */
PreciseDecimal asin(const PreciseDecimal& x);

/**
 * @brief 反余弦函数
 * @param x 输入值（范围 [-1, 1]）
 * @return acos(x)（弧度）
 */
PreciseDecimal acos(const PreciseDecimal& x);

/**
 * @brief 反正切函数
 * @param x 输入值
 * @return atan(x)（弧度）
 */
PreciseDecimal atan(const PreciseDecimal& x);

// ==================== 双曲函数 ====================

/**
 * @brief 双曲正弦函数
 * @param x 输入值
 * @return sinh(x)
 */
PreciseDecimal sinh(const PreciseDecimal& x);

/**
 * @brief 双曲余弦函数
 * @param x 输入值
 * @return cosh(x)
 */
PreciseDecimal cosh(const PreciseDecimal& x);

/**
 * @brief 双曲正切函数
 * @param x 输入值
 * @return tanh(x)
 */
PreciseDecimal tanh(const PreciseDecimal& x);

/**
 * @brief 反双曲正弦函数
 * @param x 输入值
 * @return asinh(x)
 */
PreciseDecimal asinh(const PreciseDecimal& x);

/**
 * @brief 反双曲余弦函数
 * @param x 输入值（必须 >= 1）
 * @return acosh(x)
 * @throws std::domain_error 如果 x < 1
 */
PreciseDecimal acosh(const PreciseDecimal& x);

/**
 * @brief 反双曲正切函数
 * @param x 输入值（范围 (-1, 1)）
 * @return atanh(x)
 * @throws std::domain_error 如果 |x| >= 1
 */
PreciseDecimal atanh(const PreciseDecimal& x);
} // namespace precise

// ============================================================================
// 兼容性接口（保持对旧接口的兼容）
// ============================================================================

/**
 * @brief 高精度加法
 * @param lhs 左操作数
 * @param rhs 右操作数
 * @return 两数之和
 */
PreciseDecimal add_precise_decimal(const PreciseDecimal& lhs, const PreciseDecimal& rhs);

/**
 * @brief 高精度减法
 * @param lhs 左操作数
 * @param rhs 右操作数
 * @return 两数之差
 */
PreciseDecimal subtract_precise_decimal(const PreciseDecimal& lhs, const PreciseDecimal& rhs);

/**
 * @brief 高精度乘法
 * @param lhs 左操作数
 * @param rhs 右操作数
 * @return 两数之积
 */
PreciseDecimal multiply_precise_decimal(const PreciseDecimal& lhs, const PreciseDecimal& rhs);

/**
 * @brief 高精度除法
 * @param lhs 被除数
 * @param rhs 除数
 * @return 商
 */
PreciseDecimal divide_precise_decimal(const PreciseDecimal& lhs, const PreciseDecimal& rhs);

/**
 * @brief 高精度幂运算
 * @param base 底数
 * @param exponent 指数
 * @return base^exponent
 */
PreciseDecimal pow_precise_decimal(const PreciseDecimal& base, long long exponent);

/**
 * @brief 高精度比较
 * @param lhs 左操作数
 * @param rhs 右操作数
 * @return 负数表示 lhs < rhs，零表示相等，正数表示 lhs > rhs
 */
int compare_precise_decimal(const PreciseDecimal& lhs, const PreciseDecimal& rhs);

#endif // TYPES_PRECISE_DECIMAL_H
