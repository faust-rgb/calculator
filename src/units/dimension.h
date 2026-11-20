// ============================================================================
// dimension.h - 物理量纲定义与量纲代数
// ============================================================================

#ifndef UNITS_DIMENSION_H
#define UNITS_DIMENSION_H

#include <string>
#include <array>
#include <cstdint>

namespace units {

/**
 * @struct Dimension
 * @brief 7 个 SI 国际单位制基础量纲的指数表示
 *
 * 基础量纲：
 * 0: L (长度 Length, meter)
 * 1: M (质量 Mass, kilogram)
 * 2: T (时间 Time, second)
 * 3: I (电流 Electric Current, ampere)
 * 4: Theta (热力学温度 Temperature, kelvin)
 * 5: N (物质的量 Amount of Substance, mole)
 * 6: J (发光强度 Luminous Intensity, candela)
 * 7: Angle (平面角 Angle, radian)
 */
class Dimension {
public:
    static constexpr std::size_t kNumDimensions = 8;

    enum BaseDim : std::size_t {
        L = 0,      // Length
        M = 1,      // Mass
        T = 2,      // Time
        I = 3,      // Electric Current
        Theta = 4,  // Temperature
        N = 5,      // Amount of Substance
        J = 6,      // Luminous Intensity
        Angle = 7   // Angle (rad)
    };

    std::array<int8_t, kNumDimensions> exponents{};

    // 默认构造：无量纲
    constexpr Dimension() = default;

    // 单一基础量纲构造
    constexpr explicit Dimension(BaseDim dim, int8_t power = 1) {
        exponents[static_cast<std::size_t>(dim)] = power;
    }

    // 指定全部指数
    constexpr explicit Dimension(std::array<int8_t, kNumDimensions> exps)
        : exponents(exps) {}

    // 常用预设量纲
    static constexpr Dimension Dimensionless() { return Dimension(); }
    static constexpr Dimension Length() { return Dimension(L, 1); }
    static constexpr Dimension Mass() { return Dimension(M, 1); }
    static constexpr Dimension Time() { return Dimension(T, 1); }
    static constexpr Dimension Current() { return Dimension(I, 1); }
    static constexpr Dimension Temperature() { return Dimension(Theta, 1); }
    static constexpr Dimension Amount() { return Dimension(N, 1); }
    static constexpr Dimension Luminosity() { return Dimension(J, 1); }

    // 常用复合量纲
    static constexpr Dimension Area() { return Dimension(L, 2); }
    static constexpr Dimension Volume() { return Dimension(L, 3); }
    static constexpr Dimension Velocity() {
        Dimension d; d.exponents[L] = 1; d.exponents[T] = -1; return d;
    }
    static constexpr Dimension Acceleration() {
        Dimension d; d.exponents[L] = 1; d.exponents[T] = -2; return d;
    }
    static constexpr Dimension Force() {
        Dimension d; d.exponents[M] = 1; d.exponents[L] = 1; d.exponents[T] = -2; return d;
    }
    static constexpr Dimension Energy() {
        Dimension d; d.exponents[M] = 1; d.exponents[L] = 2; d.exponents[T] = -2; return d;
    }
    static constexpr Dimension Power() {
        Dimension d; d.exponents[M] = 1; d.exponents[L] = 2; d.exponents[T] = -3; return d;
    }
    static constexpr Dimension Pressure() {
        Dimension d; d.exponents[M] = 1; d.exponents[L] = -1; d.exponents[T] = -2; return d;
    }
    static constexpr Dimension Charge() {
        Dimension d; d.exponents[I] = 1; d.exponents[T] = 1; return d;
    }
    static constexpr Dimension Voltage() {
        Dimension d; d.exponents[M] = 1; d.exponents[L] = 2; d.exponents[T] = -3; d.exponents[I] = -1; return d;
    }
    static constexpr Dimension Resistance() {
        Dimension d; d.exponents[M] = 1; d.exponents[L] = 2; d.exponents[T] = -3; d.exponents[I] = -2; return d;
    }
    static constexpr Dimension Capacitance() {
        Dimension d; d.exponents[M] = -1; d.exponents[L] = -2; d.exponents[T] = 4; d.exponents[I] = 2; return d;
    }
    static constexpr Dimension Frequency() {
        Dimension d; d.exponents[T] = -1; return d;
    }

    // 运算
    bool is_dimensionless() const;
    Dimension operator*(const Dimension& other) const;
    Dimension operator/(const Dimension& other) const;
    Dimension pow(int power) const;
    bool operator==(const Dimension& other) const;
    bool operator!=(const Dimension& other) const;

    // 格式化输出，如 "[L M T^-2]" 或 "[1]"
    std::string to_string() const;
};

} // namespace units

#endif // UNITS_DIMENSION_H
