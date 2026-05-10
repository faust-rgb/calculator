// ============================================================================
// 数据布局优化 + 10^9 基数快速操作
// ============================================================================
//
// 本文件实现：
// - 内存对齐优化（32/64 字节对齐）
// - 分支预测优化（使用 __builtin_expect）
// - 10^9 基数的逻辑捷径（快速除以 10^k）
// - 符号与绝对值分离
// ============================================================================

#include "precise_decimal.h"

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <cstring>

// 检测编译器特性
#ifdef __GNUC__
#define LIKELY(x) __builtin_expect(!!(x), 1)
#define UNLIKELY(x) __builtin_expect(!!(x), 0)
#else
#define LIKELY(x) (x)
#define UNLIKELY(x) (x)
#endif

namespace precise {

// ============================================================================
// 常量定义
// ============================================================================

constexpr uint32_t kBase = 1000000000;
constexpr uint64_t kBase64 = 1000000000ULL;
constexpr int kBaseDigits = 9;

// 内存对齐常量
constexpr size_t kCacheLineSize = 64;    // CPU 缓存行大小
constexpr size_t kSIMDAlignment = 32;    // AVX2 对齐要求
constexpr size_t kPageSize = 4096;       // 页大小

// ============================================================================
// 对齐内存分配器
// ============================================================================

class AlignedAllocator {
public:
    // 分配对齐的内存
    static uint32_t* allocate(size_t count, size_t alignment = kSIMDAlignment) {
        size_t size = count * sizeof(uint32_t);

        // 确保大小是对齐的倍数
        size = (size + alignment - 1) & ~(alignment - 1);

        void* ptr = nullptr;

#ifdef _MSC_VER
        ptr = _aligned_malloc(size, alignment);
#else
        // C++17 aligned_alloc
        ptr = std::aligned_alloc(alignment, size);
#endif

        if (!ptr) {
            throw std::bad_alloc();
        }

        return static_cast<uint32_t*>(ptr);
    }

    // 释放对齐的内存
    static void deallocate(uint32_t* ptr) {
        if (ptr) {
#ifdef _MSC_VER
            _aligned_free(ptr);
#else
            std::free(ptr);
#endif
        }
    }

    // 检查指针是否对齐
    static bool is_aligned(const void* ptr, size_t alignment) {
        return (reinterpret_cast<uintptr_t>(ptr) & (alignment - 1)) == 0;
    }
};

// ============================================================================
// 对齐的 BigIntData 包装器
// ============================================================================

class AlignedBigInt {
public:
    explicit AlignedBigInt(size_t capacity = 16, size_t alignment = kSIMDAlignment)
        : alignment_(alignment), size_(0), capacity_(capacity) {
        data_ = AlignedAllocator::allocate(capacity_, alignment_);
    }

    ~AlignedBigInt() {
        AlignedAllocator::deallocate(data_);
    }

    // 禁止拷贝
    AlignedBigInt(const AlignedBigInt&) = delete;
    AlignedBigInt& operator=(const AlignedBigInt&) = delete;

    // 允许移动
    AlignedBigInt(AlignedBigInt&& other) noexcept
        : data_(other.data_), alignment_(other.alignment_),
          size_(other.size_), capacity_(other.capacity_) {
        other.data_ = nullptr;
        other.size_ = 0;
        other.capacity_ = 0;
    }

    AlignedBigInt& operator=(AlignedBigInt&& other) noexcept {
        if (this != &other) {
            AlignedAllocator::deallocate(data_);
            data_ = other.data_;
            alignment_ = other.alignment_;
            size_ = other.size_;
            capacity_ = other.capacity_;
            other.data_ = nullptr;
            other.size_ = 0;
            other.capacity_ = 0;
        }
        return *this;
    }

    // 访问元素
    uint32_t& operator[](size_t i) { return data_[i]; }
    const uint32_t& operator[](size_t i) const { return data_[i]; }

    uint32_t* data() { return data_; }
    const uint32_t* data() const { return data_; }

    size_t size() const { return size_; }
    size_t capacity() const { return capacity_; }

    void resize(size_t new_size) {
        if (new_size > capacity_) {
            reserve(new_size);
        }
        size_ = new_size;
    }

    void reserve(size_t new_capacity) {
        if (new_capacity > capacity_) {
            uint32_t* new_data = AlignedAllocator::allocate(new_capacity, alignment_);
            std::copy(data_, data_ + size_, new_data);
            AlignedAllocator::deallocate(data_);
            data_ = new_data;
            capacity_ = new_capacity;
        }
    }

    void push_back(uint32_t val) {
        if (UNLIKELY(size_ == capacity_)) {
            reserve(capacity_ * 2);
        }
        data_[size_++] = val;
    }

    void clear() { size_ = 0; }

    bool is_aligned() const {
        return AlignedAllocator::is_aligned(data_, alignment_);
    }

private:
    uint32_t* data_;
    size_t alignment_;
    size_t size_;
    size_t capacity_;
};

// ============================================================================
// 快速除以 10^k（10^9 基数优化）
// ============================================================================
//
// 对于 10^9 基数，除以 10^k 有特殊优化：
// - 如果 k 是 9 的倍数，直接执行指针偏移
// - 如果 k < 9，只需对每个 limb 进行简单除法
// ============================================================================

// 快速除以 10^k（原地操作）
size_t divide_by_pow10_fast(uint32_t* data, size_t size, int k) {
    if (k <= 0 || size == 0 || (size == 1 && data[0] == 0)) {
        return size;
    }

    // k = chunk_shift * 9 + digit_shift
    int chunk_shift = k / kBaseDigits;
    int digit_shift = k % kBaseDigits;

    // 如果 k >= 9 * size，结果为 0
    if (chunk_shift >= static_cast<int>(size)) {
        data[0] = 0;
        return 1;
    }

    // 第一步：指针偏移（如果 k >= 9）
    size_t new_size = size - chunk_shift;
    if (chunk_shift > 0) {
        std::memmove(data, data + chunk_shift, new_size * sizeof(uint32_t));
    }

    // 第二步：处理剩余的 digit_shift 位
    if (digit_shift > 0) {
        // 计算 10^digit_shift
        uint32_t divisor = 1;
        for (int i = 0; i < digit_shift; ++i) {
            divisor *= 10;
        }

        // 从高位到低位进行除法
        uint64_t remainder = 0;
        for (int i = static_cast<int>(new_size) - 1; i >= 0; --i) {
            uint64_t current = remainder * kBase64 + data[i];
            data[i] = static_cast<uint32_t>(current / divisor);
            remainder = current % divisor;
        }
    }

    // 去除前导零
    while (new_size > 1 && data[new_size - 1] == 0) {
        --new_size;
    }

    return new_size;
}

// 快速乘以 10^k（原地操作）
size_t multiply_by_pow10_fast(uint32_t* data, size_t size, size_t capacity, int k) {
    if (k <= 0 || size == 0 || (size == 1 && data[0] == 0)) {
        return size;
    }

    int chunk_shift = k / kBaseDigits;
    int digit_shift = k % kBaseDigits;

    // 检查容量
    size_t new_size = size + chunk_shift + (digit_shift > 0 ? 1 : 0);
    if (new_size > capacity) {
        return 0;  // 失败，需要扩展
    }

    // 第一步：处理 digit_shift
    if (digit_shift > 0) {
        uint32_t multiplier = 1;
        for (int i = 0; i < digit_shift; ++i) {
            multiplier *= 10;
        }

        // 乘以 multiplier
        uint64_t carry = 0;
        for (size_t i = 0; i < size; ++i) {
            uint64_t cur = static_cast<uint64_t>(data[i]) * multiplier + carry;
            data[i] = static_cast<uint32_t>(cur % kBase);
            carry = cur / kBase;
        }

        if (carry > 0) {
            data[size++] = static_cast<uint32_t>(carry);
        }
    }

    // 第二步：指针偏移（插入 chunk_shift 个零）
    if (chunk_shift > 0) {
        std::memmove(data + chunk_shift, data, size * sizeof(uint32_t));
        std::fill(data, data + chunk_shift, 0);
        size += chunk_shift;
    }

    return size;
}

// 快速检查是否是 10^k 的倍数
bool is_multiple_of_pow10_fast(const uint32_t* data, size_t size, int k) {
    if (k <= 0) {
        return true;
    }
    if (size == 0 || (size == 1 && data[0] == 0)) {
        return true;  // 0 是任何数的倍数
    }

    int chunk_shift = k / kBaseDigits;
    int digit_shift = k % kBaseDigits;

    // 检查低 chunk_shift 个 limb 是否全为零
    for (int i = 0; i < chunk_shift && i < static_cast<int>(size); ++i) {
        if (data[i] != 0) {
            return false;
        }
    }

    // 检查 digit_shift 位
    if (digit_shift > 0 && chunk_shift < static_cast<int>(size)) {
        uint32_t divisor = 1;
        for (int i = 0; i < digit_shift; ++i) {
            divisor *= 10;
        }
        return (data[chunk_shift] % divisor) == 0;
    }

    return true;
}

// ============================================================================
// 快速计算末尾零的数量
// ============================================================================

int count_trailing_zeros_fast(const uint32_t* data, size_t size) {
    if (size == 0 || (size == 1 && data[0] == 0)) {
        return 0;
    }

    int count = 0;

    // 计算完整的零 limb
    size_t i = 0;
    while (i < size && data[i] == 0) {
        count += kBaseDigits;
        ++i;
    }

    // 计算部分零
    if (i < size) {
        uint32_t val = data[i];
        while (val > 0 && val % 10 == 0) {
            ++count;
            val /= 10;
        }
    }

    return count;
}

// ============================================================================
// 符号与绝对值分离的运算
// ============================================================================
//
// 核心运算使用绝对值数组，仅在顶层处理符号
// 避免在每个加法操作内部进行复杂的 if (negative) 分支判断
// ============================================================================

// 绝对值加法（假设 a >= b）
size_t add_abs_inplace(uint32_t* a, size_t an, const uint32_t* b, size_t bn) {
    // 进位概率约为 1/10^9，非常罕见
    uint32_t carry = 0;
    size_t i = 0;

    for (; i < bn; ++i) {
        uint64_t sum = static_cast<uint64_t>(a[i]) + b[i] + carry;
        a[i] = static_cast<uint32_t>(sum % kBase);

        // 分支预测优化：进位非常罕见
        if (UNLIKELY(sum >= kBase)) {
            carry = 1;
        } else {
            carry = 0;
        }
    }

    for (; i < an && LIKELY(carry > 0); ++i) {
        uint64_t sum = static_cast<uint64_t>(a[i]) + carry;
        a[i] = static_cast<uint32_t>(sum % kBase);
        carry = (sum >= kBase) ? 1 : 0;
    }

    if (UNLIKELY(carry > 0)) {
        a[i++] = carry;
    }

    return i;
}

// 绝对值减法（假设 a >= b）
size_t sub_abs_inplace(uint32_t* a, size_t an, const uint32_t* b, size_t bn) {
    int32_t borrow = 0;
    size_t i = 0;

    for (; i < bn; ++i) {
        int64_t diff = static_cast<int64_t>(a[i]) - borrow - b[i];

        // 分支预测优化：借位也相对罕见（约 1/10^9）
        if (UNLIKELY(diff < 0)) {
            diff += kBase;
            borrow = 1;
        } else {
            borrow = 0;
        }
        a[i] = static_cast<uint32_t>(diff);
    }

    for (; i < an && UNLIKELY(borrow > 0); ++i) {
        int64_t diff = static_cast<int64_t>(a[i]) - borrow;
        if (diff < 0) {
            diff += kBase;
            borrow = 1;
        } else {
            borrow = 0;
        }
        a[i] = static_cast<uint32_t>(diff);
    }

    // 去除前导零
    while (an > 1 && a[an - 1] == 0) {
        --an;
    }

    return an;
}

// 带符号加法（返回符号）
// result = |a| * sign_a + |b| * sign_b
int add_signed(const uint32_t* a, size_t an, bool a_neg,
               const uint32_t* b, size_t bn, bool b_neg,
               uint32_t* result, size_t& result_size) {
    if (a_neg == b_neg) {
        // 同号：直接相加
        result_size = add_abs_inplace(result, an, b, bn);
        return a_neg ? -1 : 1;
    }

    // 异号：绝对值相减
    int cmp = compare_raw(a, an, b, bn);

    if (cmp == 0) {
        result[0] = 0;
        result_size = 1;
        return 1;  // 零为正
    }

    if (cmp > 0) {
        // |a| > |b|
        std::copy(a, a + an, result);
        result_size = sub_abs_inplace(result, an, b, bn);
        return a_neg ? -1 : 1;
    } else {
        // |b| > |a|
        std::copy(b, b + bn, result);
        result_size = sub_abs_inplace(result, bn, a, an);
        return b_neg ? -1 : 1;
    }
}

// ============================================================================
// 快速规范化（去除前导零）
// ============================================================================

size_t normalize_fast(uint32_t* data, size_t size) {
    // 从后向前查找第一个非零元素
    while (size > 1 && LIKELY(data[size - 1] == 0)) {
        --size;
    }
    return size;
}

// ============================================================================
// 快速比较
// ============================================================================

int compare_fast(const uint32_t* a, size_t an, const uint32_t* b, size_t bn) {
    // 长度不同时，长度大的值大
    if (an != bn) {
        return an < bn ? -1 : 1;
    }

    // 从高位到低位比较
    for (size_t i = an; i > 0; --i) {
        if (a[i - 1] != b[i - 1]) {
            return a[i - 1] < b[i - 1] ? -1 : 1;
        }
    }

    return 0;
}

// ============================================================================
// 快速拷贝（利用对齐优化）
// ============================================================================

void copy_aligned(uint32_t* dst, const uint32_t* src, size_t count) {
    if (AlignedAllocator::is_aligned(dst, kSIMDAlignment) &&
        AlignedAllocator::is_aligned(src, kSIMDAlignment)) {
        // 使用 SIMD 加速的拷贝
        std::memcpy(dst, src, count * sizeof(uint32_t));
    } else {
        // 普通拷贝
        std::copy(src, src + count, dst);
    }
}

// ============================================================================
// 快速填充零
// ============================================================================

void fill_zero_aligned(uint32_t* dst, size_t count) {
    if (AlignedAllocator::is_aligned(dst, kSIMDAlignment)) {
        std::memset(dst, 0, count * sizeof(uint32_t));
    } else {
        std::fill(dst, dst + count, 0);
    }
}

// ============================================================================
// 导出函数
// ============================================================================

// 快速除以 10^k（BigIntData 版本）
BigIntData divide_bigint_by_pow10_fast(BigIntData v, int n) {
    size_t new_size = divide_by_pow10_fast(v.data(), v.size(), n);
    v.resize(new_size);
    return v;
}

// 快速乘以 10^k（BigIntData 版本）
BigIntData multiply_bigint_by_pow10_fast(BigIntData v, int n) {
    if (n <= 0 || v.empty() || (v.size() == 1 && v[0] == 0)) {
        return v;
    }

    int chunk_shift = n / kBaseDigits;
    int digit_shift = n % kBaseDigits;

    // 先处理 digit_shift
    if (digit_shift > 0) {
        uint32_t multiplier = 1;
        for (int i = 0; i < digit_shift; ++i) {
            multiplier *= 10;
        }
        v = multiply_bigint_by_uint32(v, multiplier);
    }

    // 再处理 chunk_shift
    if (chunk_shift > 0) {
        v.insert(v.begin(), static_cast<size_t>(chunk_shift), 0);
    }

    return v;
}

// 快速检查是否是 10^k 的倍数（BigIntData 版本）
bool is_bigint_multiple_of_pow10_fast(const BigIntData& v, int n) {
    return is_multiple_of_pow10_fast(v.data(), v.size(), n);
}

// 快速计算末尾零的数量（BigIntData 版本）
int count_trailing_zeros_bigint_fast(const BigIntData& v) {
    return count_trailing_zeros_fast(v.data(), v.size());
}

} // namespace precise