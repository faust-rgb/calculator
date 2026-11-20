// ============================================================================
// test_units.cpp - 物理单位与量纲分析测试套件
// ============================================================================

#include "calculator.h"
#include "units/dimension.h"
#include "units/unit_registry.h"
#include "units/quantity.h"
#include "units/physical_constants.h"
#include <iostream>
#include <cassert>
#include <cmath>

namespace test_suites {

int run_units_tests(int& passed, int& failed) {
    std::cout << "Running Units & Dimensional Analysis tests..." << std::endl;
    Calculator calc;

    // 1. 量纲代数测试
    try {
        units::Dimension length = units::Dimension::Length();
        units::Dimension time = units::Dimension::Time();
        units::Dimension velocity = length / time;
        units::Dimension acceleration = velocity / time;
        units::Dimension mass = units::Dimension::Mass();
        units::Dimension force = mass * acceleration;

        if (force == units::Dimension::Force() &&
            velocity == units::Dimension::Velocity() &&
            acceleration == units::Dimension::Acceleration()) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: Dimension algebra mismatch\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: Dimension test threw: " << ex.what() << "\n";
    }

    // 2. 单位注册与复合单位解析
    try {
        units::Dimension dim;
        mymath::Scalar scale;
        bool ok = units::UnitRegistry::instance().parse_unit_expression("kg*m/s^2", &dim, &scale);
        if (ok && dim == units::Dimension::Force() && scale == mymath::Scalar(1.0L)) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: parse_unit_expression kg*m/s^2\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: Unit registry test threw: " << ex.what() << "\n";
    }

    // 3. 物理量转换 (km/h -> m/s)
    try {
        std::string res = calc.evaluate_for_display("convert_unit(100, \"km/h\", \"m/s\")", false);
        // 100 km/h = 27.7777777778 m/s
        if (res.find("27.7777") != std::string::npos && res.find("m/s") != std::string::npos) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: convert_unit km/h -> m/s: " << res << "\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: convert_unit test threw: " << ex.what() << "\n";
    }

    // 4. 温度转换 (degC -> degF, degC -> K)
    try {
        std::string res1 = calc.evaluate_for_display("convert_unit(100, \"degC\", \"degF\")", false);
        std::string res2 = calc.evaluate_for_display("convert_unit(0, \"degC\", \"K\")", false);
        if (res1.find("212") != std::string::npos && res2.find("273.15") != std::string::npos) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: temperature conversion: " << res1 << ", " << res2 << "\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: temperature conversion test threw: " << ex.what() << "\n";
    }

    // 5. 物理量运算 (乘法与除法生成新量纲)
    try {
        std::string res_mul = calc.evaluate_for_display("unit_mul(\"50 m/s\", \"10 s\")", false);
        std::string res_div = calc.evaluate_for_display("unit_div(\"100 J\", \"5 s\")", false);
        if (res_mul.find("500") != std::string::npos && res_div.find("20") != std::string::npos &&
            res_div.find("W") != std::string::npos) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: unit_mul/div: " << res_mul << ", " << res_div << "\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: unit_mul/div test threw: " << ex.what() << "\n";
    }

    // 6. 量纲不匹配报错保护 (10 m + 5 s)
    try {
        calc.evaluate_for_display("unit_add(\"10 m\", \"5 s\")", false);
        ++failed;
        std::cout << "FAIL: unit_add with incompatible dimensions should throw\n";
    } catch (const std::exception&) {
        ++passed;
    }

    // 7. 物理常数查询
    try {
        std::string c_val = calc.evaluate_for_display("c_light()", false);
        std::string g_val = calc.evaluate_for_display("g_earth()", false);
        if (c_val.find("299792458") != std::string::npos && g_val.find("9.80665") != std::string::npos) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: physical constants: " << c_val << ", " << g_val << "\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: physical constants test threw: " << ex.what() << "\n";
    }

    // 8. 命令测试 (:convert, :constants, :units)
    try {
        std::string out_conv = calc.process_line(":convert 100 km/h to m/s", false);
        std::string out_const = calc.process_line(":constants", false);
        std::string out_units = calc.process_line(":units", false);
        if (out_conv.find("27.7777") != std::string::npos &&
            out_const.find("Speed of light") != std::string::npos &&
            out_units.find("meter") != std::string::npos) {
            ++passed;
        } else {
            ++failed;
            std::cout << "FAIL: units commands output unexpected\n"
                      << "  out_conv: " << out_conv << "\n"
                      << "  out_const: " << out_const << "\n"
                      << "  out_units: " << out_units << "\n";
        }
    } catch (const std::exception& ex) {
        ++failed;
        std::cout << "FAIL: units commands test threw: " << ex.what() << "\n";
    }

    return 0;
}

} // namespace test_suites
