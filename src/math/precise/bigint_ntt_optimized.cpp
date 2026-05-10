// ============================================================================
// NTT (数论变换) 优化实现
// ============================================================================
//
// 本文件实现高效的 NTT 算法用于大整数乘法：
// - 预计算并缓存 twiddle factors
// - 多模数 NTT + CRT 合并（适配 10^9 基数）
// - 并行执行三次 NTT 变换
// - 内存池优化
//
// 对于 10^9 基数，NTT 后的系数可能达到 N * (10^9)^2
// 使用三个 32 位素数进行 NTT，然后通过 CRT 合并
// ============================================================================

#include "precise_decimal.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <mutex>
#include <thread>
#include <vector>

// 检测 CPU 支持的指令集
#ifdef __AVX2__
#include <immintrin.h>
#define HAS_AVX2 1
#else
#define HAS_AVX2 0
#endif

namespace precise {

// ============================================================================
// 常量定义
// ============================================================================

constexpr uint32_t kBase = 1000000000;
constexpr int kBaseDigits = 9;

// NTT 模数选择（支持 NTT 的特殊素数）
// 这些素数满足：p = k * 2^n + 1，且原根为 3
namespace ntt_mods {

// 模数 1: 998244353 = 119 * 2^23 + 1，支持最大长度 2^23
constexpr uint32_t P1 = 998244353;
constexpr uint32_t G1 = 3;  // 原根

// 模数 2: 1004535809 = 479 * 2^21 + 1，支持最大长度 2^21
constexpr uint32_t P2 = 1004535809;
constexpr uint32_t G2 = 3;

// 模数 3: 469762049 = 7 * 2^26 + 1，支持最大长度 2^26
constexpr uint32_t P3 = 469762049;
constexpr uint32_t G3 = 3;

// 模数乘积的位数上限
// P1 * P2 * P3 ≈ 4.7 * 10^26，足够容纳 N * (10^9)^2
constexpr uint64_t M12 = static_cast<uint64_t>(P1) * P2;
// M12 * P3 ≈ 4.7 * 10^26

} // namespace ntt_mods

// ============================================================================
// Twiddle Factors 缓存
// ============================================================================
//
// 预计算并缓存不同长度的 twiddle factors，避免重复计算
// 使用 thread_local 确保线程安全
// ============================================================================

template<uint32_t P, uint32_t G>
class TwiddleCache {
public:
    struct CacheEntry {
        std::vector<uint32_t> forward;   // 正向变换的 twiddle factors
        std::vector<uint32_t> inverse;   // 逆向变换的 twiddle factors
        size_t length;                   // 支持的长度
    };

    static TwiddleCache& instance() {
        static thread_local TwiddleCache cache;
        return cache;
    }

    // 获取指定长度的 twiddle factors
    const CacheEntry& get(size_t n) {
        // 查找缓存
        for (const auto& entry : cache_) {
            if (entry.length == n) {
                return entry;
            }
        }

        // 未找到，创建新的缓存条目
        cache_.push_back(create_entry(n));
        return cache_.back();
    }

    // 清除缓存
    void clear() {
        cache_.clear();
    }

private:
    TwiddleCache() = default;

    std::vector<CacheEntry> cache_;

    static uint32_t power(uint32_t a, uint32_t b, uint32_t mod) {
        uint32_t res = 1;
        a %= mod;
        while (b > 0) {
            if (b & 1) {
                res = static_cast<uint32_t>((static_cast<uint64_t>(res) * a) % mod);
            }
            a = static_cast<uint32_t>((static_cast<uint64_t>(a) * a) % mod);
            b >>= 1;
        }
        return res;
    }

    static uint32_t inv(uint32_t n, uint32_t mod) {
        return power(n, mod - 2, mod);
    }

    CacheEntry create_entry(size_t n) {
        CacheEntry entry;
        entry.length = n;
        entry.forward.reserve(n / 2);
        entry.inverse.reserve(n / 2);

        // 预计算所有层级的 twiddle factors
        for (size_t len = 2; len <= n; len <<= 1) {
            uint32_t wlen = power(G, (P - 1) / static_cast<uint32_t>(len), P);
            uint32_t wlen_inv = inv(wlen, P);

            // 存储该层级的 twiddle factors
            uint32_t w = 1;
            uint32_t w_inv = 1;
            for (size_t j = 0; j < len / 2; ++j) {
                entry.forward.push_back(w);
                entry.inverse.push_back(w_inv);
                w = static_cast<uint32_t>((static_cast<uint64_t>(w) * wlen) % P);
                w_inv = static_cast<uint32_t>((static_cast<uint64_t>(w_inv) * wlen_inv) % P);
            }
        }

        return entry;
    }
};

// ============================================================================
// NTT 变换核心实现
// ============================================================================

template<uint32_t P, uint32_t G>
class NTTCore {
public:
    static uint32_t power(uint32_t a, uint32_t b) {
        uint32_t res = 1;
        a %= P;
        while (b > 0) {
            if (b & 1) {
                res = static_cast<uint32_t>((static_cast<uint64_t>(res) * a) % P);
            }
            a = static_cast<uint32_t>((static_cast<uint64_t>(a) * a) % P);
            b >>= 1;
        }
        return res;
    }

    static uint32_t inv(uint32_t n) {
        return power(n, P - 2);
    }

    // 位反转排列（Bit-reversal permutation）
    static void bit_reverse(uint32_t* a, size_t n) {
        for (size_t i = 1, j = 0; i < n; ++i) {
            size_t bit = n >> 1;
            for (; j & bit; bit >>= 1) {
                j ^= bit;
            }
            j ^= bit;
            if (i < j) {
                std::swap(a[i], a[j]);
            }
        }
    }

    // NTT 正向变换（使用缓存的 twiddle factors）
    static void transform(uint32_t* a, size_t n, bool use_cache = true) {
        bit_reverse(a, n);

        if (use_cache) {
            // 使用缓存的 twiddle factors
            auto& cache = TwiddleCache<P, G>::instance().get(n);
            size_t idx = 0;

            for (size_t len = 2; len <= n; len <<= 1) {
                for (size_t i = 0; i < n; i += len) {
                    for (size_t j = 0; j < len / 2; ++j) {
                        uint32_t u = a[i + j];
                        uint32_t v = static_cast<uint32_t>(
                            (static_cast<uint64_t>(a[i + j + len / 2]) * cache.forward[idx + j]) % P);
                        a[i + j] = (u + v < P) ? (u + v) : (u + v - P);
                        a[i + j + len / 2] = (u >= v) ? (u - v) : (u + P - v);
                    }
                }
                idx += len / 2;
            }
        } else {
            // 不使用缓存（动态计算）
            for (size_t len = 2; len <= n; len <<= 1) {
                uint32_t wlen = power(G, (P - 1) / static_cast<uint32_t>(len));
                for (size_t i = 0; i < n; i += len) {
                    uint32_t w = 1;
                    for (size_t j = 0; j < len / 2; ++j) {
                        uint32_t u = a[i + j];
                        uint32_t v = static_cast<uint32_t>(
                            (static_cast<uint64_t>(a[i + j + len / 2]) * w) % P);
                        a[i + j] = (u + v < P) ? (u + v) : (u + v - P);
                        a[i + j + len / 2] = (u >= v) ? (u - v) : (u + P - v);
                        w = static_cast<uint32_t>((static_cast<uint64_t>(w) * wlen) % P);
                    }
                }
            }
        }
    }

    // NTT 逆向变换
    static void inverse_transform(uint32_t* a, size_t n, bool use_cache = true) {
        bit_reverse(a, n);

        if (use_cache) {
            auto& cache = TwiddleCache<P, G>::instance().get(n);
            size_t idx = 0;

            for (size_t len = 2; len <= n; len <<= 1) {
                for (size_t i = 0; i < n; i += len) {
                    for (size_t j = 0; j < len / 2; ++j) {
                        uint32_t u = a[i + j];
                        uint32_t v = static_cast<uint32_t>(
                            (static_cast<uint64_t>(a[i + j + len / 2]) * cache.inverse[idx + j]) % P);
                        a[i + j] = (u + v < P) ? (u + v) : (u + v - P);
                        a[i + j + len / 2] = (u >= v) ? (u - v) : (u + P - v);
                    }
                }
                idx += len / 2;
            }
        } else {
            for (size_t len = 2; len <= n; len <<= 1) {
                uint32_t wlen = inv(power(G, (P - 1) / static_cast<uint32_t>(len)));
                for (size_t i = 0; i < n; i += len) {
                    uint32_t w = 1;
                    for (size_t j = 0; j < len / 2; ++j) {
                        uint32_t u = a[i + j];
                        uint32_t v = static_cast<uint32_t>(
                            (static_cast<uint64_t>(a[i + j + len / 2]) * w) % P);
                        a[i + j] = (u + v < P) ? (u + v) : (u + v - P);
                        a[i + j + len / 2] = (u >= v) ? (u - v) : (u + P - v);
                        w = static_cast<uint32_t>((static_cast<uint64_t>(w) * wlen) % P);
                    }
                }
            }
        }

        // 乘以 n 的逆元
        uint32_t n_inv = inv(static_cast<uint32_t>(n));
        for (size_t i = 0; i < n; ++i) {
            a[i] = static_cast<uint32_t>((static_cast<uint64_t>(a[i]) * n_inv) % P);
        }
    }
};

// ============================================================================
// AVX2 优化的 NTT Butterfly 操作
// ============================================================================

#if HAS_AVX2

template<uint32_t P>
inline void butterfly_avx2(uint32_t* a, size_t len, uint32_t wlen) {
    // 使用 AVX2 进行 8 个 butterfly 操作
    // 注意：由于模运算的特殊性，需要小心处理溢出

    __m256i p_vec = _mm256_set1_epi32(P);

    for (size_t i = 0; i < len / 2; i += 8) {
        // 加载 u 和 v
        __m256i u = _mm256_loadu_si256(reinterpret_cast<__m256i*>(a + i));
        __m256i v = _mm256_loadu_si256(reinterpret_cast<__m256i*>(a + i + len / 2));

        // 计算 v * w（需要逐元素处理，因为 w 不同）
        // 这里简化处理，假设 w 相同（对于某些情况）
        // 实际实现需要更复杂的处理

        // 加法：u + v
        __m256i sum = _mm256_add_epi32(u, v);
        // 减法：u - v
        __m256i diff = _mm256_sub_epi32(u, v);

        // 模运算（简化版，需要更精确的实现）
        // 使用条件选择处理溢出

        _mm256_storeu_si256(reinterpret_cast<__m256i*>(a + i), sum);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(a + i + len / 2), diff);
    }
}

#endif

// ============================================================================
// CRT (中国剩余定理) 合并
// ============================================================================
//
// 将三个模数下的结果合并为单一结果
// 针对 10^9 基数进行优化
// ============================================================================

BigIntData crt_merge_optimized(const uint32_t* r1, const uint32_t* r2,
                               const uint32_t* r3, size_t n) {
    using namespace ntt_mods;

    BigIntData res;
    res.reserve(n + 4);

    // 预计算 CRT 常数
    const uint32_t invP1_P2 = NTTCore<P2, G2>::inv(P1);
    const uint32_t invM12_P3 = NTTCore<P3, G3>::inv(static_cast<uint32_t>(M12 % P3));

    // 使用 unsigned __int128 处理大数运算
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
    unsigned __int128 carry = 0;
#pragma GCC diagnostic pop

    for (size_t i = 0; i < n; ++i) {
        // 第一步：合并 P1 和 P2
        // k1 = (r2 - r1) * inv(P1) mod P2
        uint32_t k1 = static_cast<uint32_t>(
            (static_cast<uint64_t>(r2[i] + P2 - r1[i]) * invP1_P2) % P2);

        // r12 = r1 + k1 * P1
        uint64_t r12 = r1[i] + static_cast<uint64_t>(k1) * P1;

        // 第二步：合并 P3
        // k2 = (r3 - r12 mod P3) * inv(M12) mod P3
        uint32_t r12_mod_p3 = static_cast<uint32_t>(r12 % P3);
        uint32_t k2 = static_cast<uint32_t>(
            (static_cast<uint64_t>(r3[i] + P3 - r12_mod_p3) * invM12_P3) % P3);

        // 最终值：r12 + k2 * M12
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
        unsigned __int128 val = carry + r12 + (unsigned __int128)k2 * M12;
#pragma GCC diagnostic pop

        // 转换为 10^9 基数
        res.push_back(static_cast<uint32_t>(val % kBase));
        carry = val / kBase;
    }

    // 处理剩余进位
    while (carry > 0) {
        res.push_back(static_cast<uint32_t>(carry % kBase));
        carry /= kBase;
    }

    // 去除前导零
    while (res.size() > 1 && res.back() == 0) {
        res.pop_back();
    }

    return res;
}

// ============================================================================
// 并行 NTT 变换
// ============================================================================

void parallel_ntt_transform(uint32_t* fa1, uint32_t* fa2, uint32_t* fa3,
                            uint32_t* fb1, uint32_t* fb2, uint32_t* fb3,
                            size_t n) {
    using namespace ntt_mods;

    // 如果数据量足够大，使用多线程并行执行
    if (n >= 4096) {
        std::thread t1([&]() {
            NTTCore<P1, G1>::transform(fa1, n);
            NTTCore<P1, G1>::transform(fb1, n);
        });

        std::thread t2([&]() {
            NTTCore<P2, G2>::transform(fa2, n);
            NTTCore<P2, G2>::transform(fb2, n);
        });

        // 当前线程处理第三个模数
        NTTCore<P3, G3>::transform(fa3, n);
        NTTCore<P3, G3>::transform(fb3, n);

        t1.join();
        t2.join();
    } else {
        // 小数据量，单线程顺序执行
        NTTCore<P1, G1>::transform(fa1, n);
        NTTCore<P1, G1>::transform(fb1, n);
        NTTCore<P2, G2>::transform(fa2, n);
        NTTCore<P2, G2>::transform(fb2, n);
        NTTCore<P3, G3>::transform(fa3, n);
        NTTCore<P3, G3>::transform(fb3, n);
    }
}

void parallel_ntt_inverse(uint32_t* fa1, uint32_t* fa2, uint32_t* fa3, size_t n) {
    using namespace ntt_mods;

    if (n >= 4096) {
        std::thread t1([&]() { NTTCore<P1, G1>::inverse_transform(fa1, n); });
        std::thread t2([&]() { NTTCore<P2, G2>::inverse_transform(fa2, n); });
        NTTCore<P3, G3>::inverse_transform(fa3, n);
        t1.join();
        t2.join();
    } else {
        NTTCore<P1, G1>::inverse_transform(fa1, n);
        NTTCore<P2, G2>::inverse_transform(fa2, n);
        NTTCore<P3, G3>::inverse_transform(fa3, n);
    }
}

// ============================================================================
// 优化的 NTT 乘法入口函数
// ============================================================================

BigIntData multiply_bigint_ntt_optimized(const BigIntData& lhs, const BigIntData& rhs) {
    using namespace ntt_mods;

    if ((lhs.size() == 1 && lhs[0] == 0) || (rhs.size() == 1 && rhs[0] == 0)) {
        return {0};
    }

    // 计算所需的 NTT 长度（必须是 2 的幂）
    size_t n = 1;
    while (n < lhs.size() + rhs.size()) {
        n <<= 1;
    }

    // 使用内存池
    static thread_local BigIntArena arena;
    const size_t required_words = n * 6 + lhs.size() + rhs.size() + 64;
    if (arena.capacity() < required_words) {
        arena = BigIntArena(required_words);
    }
    arena.reset();

    // 分配内存
    uint32_t* fa1 = arena.allocate(n, 0);
    uint32_t* fa2 = arena.allocate(n, 0);
    uint32_t* fa3 = arena.allocate(n, 0);
    uint32_t* fb1 = arena.allocate(n, 0);
    uint32_t* fb2 = arena.allocate(n, 0);
    uint32_t* fb3 = arena.allocate(n, 0);

    // 复制数据
    std::copy(lhs.begin(), lhs.end(), fa1);
    std::copy(lhs.begin(), lhs.end(), fa2);
    std::copy(lhs.begin(), lhs.end(), fa3);
    std::copy(rhs.begin(), rhs.end(), fb1);
    std::copy(rhs.begin(), rhs.end(), fb2);
    std::copy(rhs.begin(), rhs.end(), fb3);

    // 并行 NTT 变换
    parallel_ntt_transform(fa1, fa2, fa3, fb1, fb2, fb3, n);

    // 点乘（逐元素乘法）
    for (size_t i = 0; i < n; ++i) {
        fa1[i] = static_cast<uint32_t>((static_cast<uint64_t>(fa1[i]) * fb1[i]) % P1);
        fa2[i] = static_cast<uint32_t>((static_cast<uint64_t>(fa2[i]) * fb2[i]) % P2);
        fa3[i] = static_cast<uint32_t>((static_cast<uint64_t>(fa3[i]) * fb3[i]) % P3);
    }

    // 并行逆变换
    parallel_ntt_inverse(fa1, fa2, fa3, n);

    // CRT 合并
    BigIntData result = crt_merge_optimized(fa1, fa2, fa3, n);

    return result;
}

// ============================================================================
// 自适应 NTT 长度选择
// ============================================================================
//
// 根据输入大小选择最优的 NTT 长度
// 考虑缓存命中率和计算效率
// ============================================================================

size_t select_optimal_ntt_length(size_t required_size) {
    // 预定义的优化长度（与缓存对齐）
    static const size_t optimal_lengths[] = {
        256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072
    };

    for (size_t len : optimal_lengths) {
        if (len >= required_size) {
            return len;
        }
    }

    // 如果超出预定义范围，计算下一个 2 的幂
    size_t n = 1;
    while (n < required_size) {
        n <<= 1;
    }
    return n;
}

// ============================================================================
// 缓存预热
// ============================================================================
//
// 在程序启动时预计算常用的 twiddle factors
// ============================================================================

void warmup_ntt_cache() {
    using namespace ntt_mods;

    // 预热常用长度的缓存
    static const size_t warmup_lengths[] = {256, 512, 1024, 2048, 4096};

    for (size_t n : warmup_lengths) {
        TwiddleCache<P1, G1>::instance().get(n);
        TwiddleCache<P2, G2>::instance().get(n);
        TwiddleCache<P3, G3>::instance().get(n);
    }
}

} // namespace precise