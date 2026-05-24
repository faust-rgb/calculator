// ============================================================================
// SIMD 优化的大整数算术操作
// ============================================================================
//
// 本文件实现使用 AVX2/SSE 指令集优化的高精度算术操作：
// - 向量化加法/减法（延迟进位）
// - 内存对齐优化
// - 分支预测优化
//
// 基数为 10^9，每个 limb 占用约 30 位
// ============================================================================

#include "precise_decimal.h"

#include <algorithm>
#include <cstdint>
#include <cstring>

// 检测 CPU 支持的指令集
#ifdef __AVX2__
#include <immintrin.h>
#define HAS_AVX2 1
#else
#define HAS_AVX2 0
#endif

#ifdef __SSE2__
#include <emmintrin.h>
#define HAS_SSE2 1
#else
#define HAS_SSE2 0
#endif

namespace precise {

// ============================================================================
// 常量定义
// ============================================================================

constexpr uint32_t kBase = 1000000000;
constexpr uint64_t kBase64 = 1000000000ULL;
constexpr int kBaseDigits = 9;

// 内存对齐常量
constexpr size_t kAlignment = 32;  // 32 字节对齐（AVX2）
constexpr size_t kAlignmentMask = kAlignment - 1;

// ============================================================================
// 分支预测优化宏
// ============================================================================

// 对于 10^9 基数，进位概率约为 1/10^9，非常罕见
#define LIKELY_NO_CARRY(x) __builtin_expect((x), 0)
#define LIKELY_CARRY(x) __builtin_expect((x), 1)

// ============================================================================
// 对齐内存辅助函数
// ============================================================================

inline bool is_aligned(const void* ptr) {
    return (reinterpret_cast<uintptr_t>(ptr) & kAlignmentMask) == 0;
}

inline uint32_t* aligned_alloc(size_t n) {
    void* ptr = nullptr;
#ifdef _MSC_VER
    ptr = _aligned_malloc(n * sizeof(uint32_t), kAlignment);
#else
    ptr = std::aligned_alloc(kAlignment, ((n * sizeof(uint32_t) + kAlignmentMask) & ~kAlignmentMask));
#endif
    return static_cast<uint32_t*>(ptr);
}

inline void aligned_free(uint32_t* ptr) {
    if (ptr) {
#ifdef _MSC_VER
        _aligned_free(ptr);
#else
        std::free(ptr);
#endif
    }
}

// ============================================================================
// 延迟进位加法 (Delayed Carry Addition)
// ============================================================================
//
// 策略：
// 1. 使用 uint64_t 存储中间结果，允许暂时超过 10^9
// 2. 执行多次加法后再统一进行进位清理
// 3. SIMD 并行处理不带进位的加法
// ============================================================================

// 单次加法，结果存储在 uint64_t 数组中（允许超过 10^9）
void add_raw_delayed_carry(const uint32_t* a, size_t an,
                           const uint32_t* b, size_t bn,
                           uint64_t* res, size_t& res_size) {
    res_size = std::max(an, bn);

#if HAS_AVX2
    // AVX2 向量化：每次处理 4 个 uint64_t（256 位）
    size_t i = 0;
    size_t vec_end = res_size / 4 * 4;

    // 预处理：将 uint32_t 扩展为 uint64_t
    // 使用 SIMD 进行并行加法
    if (is_aligned(a) && is_aligned(b) && is_aligned(res)) {
        for (; i < vec_end; i += 4) {
            __m128i a_lo = _mm_loadu_si128(reinterpret_cast<const __m128i*>(a + i));
            __m128i b_lo = _mm_loadu_si128(reinterpret_cast<const __m128i*>(b + i));

            // 将 32 位扩展为 64 位并相加
            __m256i a_256 = _mm256_cvtepu32_epi64(a_lo);
            __m256i b_256 = _mm256_cvtepu32_epi64(b_lo);
            __m256i sum = _mm256_add_epi64(a_256, b_256);

            _mm256_store_si256(reinterpret_cast<__m256i*>(res + i), sum);
        }
    }
#endif

    // 处理剩余元素（标量）
    for (size_t i = 0; i < res_size; ++i) {
        uint32_t a_val = (i < an) ? a[i] : 0;
        uint32_t b_val = (i < bn) ? b[i] : 0;
        res[i] = static_cast<uint64_t>(a_val) + b_val;
    }
}

// 进位清理：将 uint64_t 数组转换为规范的 uint32_t 数组（基数为 10^9）
size_t propagate_carry(uint64_t* data, size_t size, uint32_t* out) {
    uint64_t carry = 0;
    size_t out_size = 0;

    for (size_t i = 0; i < size; ++i) {
        uint64_t val = data[i] + carry;
        out[i] = static_cast<uint32_t>(val % kBase64);
        carry = val / kBase64;
        out_size = i + 1;
    }

    // 处理剩余进位
    while (carry > 0) {
        if (out_size >= size) {
            // 需要扩展输出数组
            out[out_size] = static_cast<uint32_t>(carry % kBase64);
            carry /= kBase64;
            ++out_size;
        } else {
            out[out_size] = static_cast<uint32_t>(carry % kBase64);
            carry /= kBase64;
            ++out_size;
        }
    }

    // 去除前导零
    while (out_size > 1 && out[out_size - 1] == 0) {
        --out_size;
    }

    return out_size;
}

// ============================================================================
// SIMD 优化的加法（完整实现）
// ============================================================================

size_t add_raw_simd(const uint32_t* a, size_t an,
                    const uint32_t* b, size_t bn,
                    uint32_t* res) {
    // 使用延迟进位策略
    size_t max_size = std::max(an, bn) + 1;

    // 临时 uint64_t 数组存储中间结果
    uint64_t* temp = new uint64_t[max_size];
    size_t temp_size;

    add_raw_delayed_carry(a, an, b, bn, temp, temp_size);

    // 进位清理
    size_t result_size = propagate_carry(temp, temp_size, res);

    delete[] temp;
    return result_size;
}

// ============================================================================
// SIMD 优化的减法
// ============================================================================

size_t sub_raw_simd(const uint32_t* a, size_t an,
                    const uint32_t* b, size_t bn,
                    uint32_t* res) {
    // 减法不适合延迟进位，直接使用标量实现
    // 但可以利用 SIMD 进行并行比较和掩码操作

    int64_t borrow = 0;
    size_t i = 0;

#if HAS_AVX2
    // AVX2 优化：批量处理借位检查
    // 由于借位是串行依赖的，SIMD 优化有限
    // 主要用于批量加载和存储
#endif

    for (; i < bn; ++i) {
        int64_t sub = static_cast<int64_t>(a[i]) - borrow - b[i];
        if (sub < 0) {
            sub += kBase;
            borrow = 1;
        } else {
            borrow = 0;
        }
        res[i] = static_cast<uint32_t>(sub);
    }

    for (; i < an; ++i) {
        int64_t sub = static_cast<int64_t>(a[i]) - borrow;
        if (sub < 0) {
            sub += kBase;
            borrow = 1;
        } else {
            borrow = 0;
        }
        res[i] = static_cast<uint32_t>(sub);
    }

    // 去除前导零
    while (i > 1 && res[i - 1] == 0) {
        --i;
    }

    return i;
}

// ============================================================================
// 批量加法优化（用于矩阵乘法、点积等）
// ============================================================================
//
// 执行多次加法，只在最后进行一次进位清理
// ============================================================================

class BatchAdder {
public:
    BatchAdder(size_t capacity) : capacity_(capacity) {
        // 分配对齐的临时缓冲区
        temp_ = new uint64_t[capacity + 16];
        std::fill(temp_, temp_ + capacity + 16, 0);
        size_ = 0;
    }

    ~BatchAdder() {
        delete[] temp_;
    }

    // 添加一个数（不立即处理进位）
    void add(const uint32_t* data, size_t n) {
        if (n > size_) {
            // 扩展缓冲区
            for (size_t i = size_; i < n; ++i) {
                temp_[i] = 0;
            }
            size_ = n;
        }

        // 累加（允许超过 10^9）
        for (size_t i = 0; i < n; ++i) {
            temp_[i] += data[i];
        }
    }

    // 乘加操作：temp += a * b
    void multiply_add(const uint32_t* a, size_t an, uint32_t b) {
        if (an > size_) {
            for (size_t i = size_; i < an; ++i) {
                temp_[i] = 0;
            }
            size_ = an;
        }

        for (size_t i = 0; i < an; ++i) {
            temp_[i] += static_cast<uint64_t>(a[i]) * b;
        }
    }

    // 最终进位清理，输出结果
    size_t finalize(uint32_t* out) {
        uint64_t carry = 0;
        size_t out_size = 0;

        for (size_t i = 0; i < size_; ++i) {
            uint64_t val = temp_[i] + carry;
            out[i] = static_cast<uint32_t>(val % kBase64);
            carry = val / kBase64;
            out_size = i + 1;
        }

        while (carry > 0) {
            out[out_size++] = static_cast<uint32_t>(carry % kBase64);
            carry /= kBase64;
        }

        // 去除前导零
        while (out_size > 1 && out[out_size - 1] == 0) {
            --out_size;
        }

        // 重置缓冲区
        std::fill(temp_, temp_ + size_, 0);
        size_ = 0;

        return out_size;
    }

private:
    uint64_t* temp_;
    size_t size_;
    size_t capacity_;
};

// ============================================================================
// 快速比较（SIMD 优化）
// ============================================================================

int compare_raw_simd(const uint32_t* a, size_t an, const uint32_t* b, size_t bn) {
    if (an != bn) return an < bn ? -1 : 1;

#if HAS_AVX2
    // 从高位开始比较，使用 SIMD 批量比较
    size_t remaining = an;
    size_t offset = 0;

    while (remaining >= 8) {
        size_t idx = an - offset - 8;
        __m256i a_vec = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(a + idx));
        __m256i b_vec = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(b + idx));

        // 比较是否相等
        __m256i cmp = _mm256_cmpeq_epi32(a_vec, b_vec);
        int mask = _mm256_movemask_epi8(cmp);

        if (mask != -1) {  // 不全相等
            // 找到第一个不相等的位置
            for (int j = 7; j >= 0; --j) {
                if (a[idx + j] != b[idx + j]) {
                    return a[idx + j] < b[idx + j] ? -1 : 1;
                }
            }
        }

        offset += 8;
        remaining -= 8;
    }
#endif

    // 标量处理剩余部分
    for (size_t i = an; i > 0; --i) {
        if (a[i - 1] != b[i - 1]) {
            return a[i - 1] < b[i - 1] ? -1 : 1;
        }
    }

    return 0;
}

// ============================================================================
// 快速规范化（去除前导零）
// ============================================================================

size_t normalize_raw(uint32_t* data, size_t size) {
    // 使用 SIMD 找到最后一个非零元素
#if HAS_AVX2
    if (size >= 8) {
        size_t idx = size;
        while (idx >= 8) {
            idx -= 8;
            __m256i vec = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(data + idx));
            __m256i zero = _mm256_setzero_si256();
            __m256i cmp = _mm256_cmpeq_epi32(vec, zero);
            int mask = _mm256_movemask_epi8(cmp);

            if (mask != -1) {  // 存在非零元素
                // 找到最后一个非零位置
                for (int j = 7; j >= 0; --j) {
                    if (data[idx + j] != 0) {
                        return idx + j + 1;
                    }
                }
            }
        }
        // 检查剩余部分
        while (idx > 0 && data[idx - 1] == 0) {
            --idx;
        }
        return idx > 0 ? idx : 1;
    }
#endif

    // 标量实现
    while (size > 1 && data[size - 1] == 0) {
        --size;
    }
    return size;
}

// ============================================================================
// 导出函数接口
// ============================================================================

size_t add_raw_optimized(const uint32_t* a, size_t an,
                         const uint32_t* b, size_t bn,
                         uint32_t* res) {
#if HAS_AVX2
    return add_raw_simd(a, an, b, bn, res);
#else
    // 回退到标量实现
    constexpr uint32_t kBase = 1000000000;
    uint32_t carry = 0;
    size_t i = 0;
    size_t min_n = std::min(an, bn);
    for (; i < min_n; ++i) {
        uint64_t sum = static_cast<uint64_t>(a[i]) + b[i] + carry;
        res[i] = static_cast<uint32_t>(sum % kBase);
        carry = static_cast<uint32_t>(sum / kBase);
    }
    const uint32_t* longer = (an > bn) ? a : b;
    size_t max_n = std::max(an, bn);
    for (; i < max_n; ++i) {
        uint64_t sum = static_cast<uint64_t>(longer[i]) + carry;
        res[i] = static_cast<uint32_t>(sum % kBase);
        carry = static_cast<uint32_t>(sum / kBase);
    }
    if (carry) {
        res[i++] = carry;
    }
    return i;
#endif
}

size_t sub_raw_optimized(const uint32_t* a, size_t an,
                         const uint32_t* b, size_t bn,
                         uint32_t* res) {
#if HAS_AVX2
    return sub_raw_simd(a, an, b, bn, res);
#else
    // 回退到标量实现
    constexpr uint32_t kBase = 1000000000;
    int32_t borrow = 0;
    size_t i = 0;
    for (; i < bn; ++i) {
        int64_t sub = static_cast<int64_t>(a[i]) - borrow - b[i];
        if (sub < 0) {
            sub += kBase;
            borrow = 1;
        } else {
            borrow = 0;
        }
        res[i] = static_cast<uint32_t>(sub);
    }
    for (; i < an; ++i) {
        int64_t sub = static_cast<int64_t>(a[i]) - borrow;
        if (sub < 0) {
            sub += kBase;
            borrow = 1;
        } else {
            borrow = 0;
        }
        res[i] = static_cast<uint32_t>(sub);
    }
    while (i > 1 && res[i - 1] == 0) i--;
    return i;
#endif
}

int compare_raw_optimized(const uint32_t* a, size_t an,
                          const uint32_t* b, size_t bn) {
#if HAS_AVX2
    return compare_raw_simd(a, an, b, bn);
#else
    if (an != bn) return an < bn ? -1 : 1;
    for (size_t i = an; i > 0; --i) {
        if (a[i - 1] != b[i - 1]) {
            return a[i - 1] < b[i - 1] ? -1 : 1;
        }
    }
    return 0;
#endif
}

} // namespace precise