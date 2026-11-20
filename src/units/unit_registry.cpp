// ============================================================================
// unit_registry.cpp - 单位注册表实现与单位表达式解析
// ============================================================================

#include "units/unit_registry.h"
#include "core/services/string_utils.h"
#include <algorithm>
#include <cctype>
#include <sstream>
#include <stdexcept>

namespace units {

UnitRegistry::UnitRegistry() {
    initialize_standard_units();
}

UnitRegistry& UnitRegistry::instance() {
    static UnitRegistry reg;
    return reg;
}

void UnitRegistry::register_unit(UnitInfo unit) {
    const std::size_t idx = units_list_.size();
    units_list_.push_back(unit);
    symbol_map_[unit.symbol] = idx;
    if (!unit.name.empty()) {
        symbol_map_[unit.name] = idx;
    }
}

void UnitRegistry::initialize_standard_units() {
    // ------------------------------------------------------------------------
    // 1. SI 基础单位
    // ------------------------------------------------------------------------
    register_unit({"meter", "m", "Length", Dimension::Length(), 1.0L, 0.0L, "SI base unit of length"});
    register_unit({"kilogram", "kg", "Mass", Dimension::Mass(), 1.0L, 0.0L, "SI base unit of mass"});
    register_unit({"second", "s", "Time", Dimension::Time(), 1.0L, 0.0L, "SI base unit of time"});
    register_unit({"ampere", "A", "Current", Dimension::Current(), 1.0L, 0.0L, "SI base unit of electric current"});
    register_unit({"kelvin", "K", "Temperature", Dimension::Temperature(), 1.0L, 0.0L, "SI base unit of thermodynamic temperature"});
    register_unit({"mole", "mol", "Amount", Dimension::Amount(), 1.0L, 0.0L, "SI base unit of amount of substance"});
    register_unit({"candela", "cd", "Luminosity", Dimension::Luminosity(), 1.0L, 0.0L, "SI base unit of luminous intensity"});
    register_unit({"radian", "rad", "Angle", Dimension(Dimension::Angle, 1), 1.0L, 0.0L, "SI unit of plane angle"});
    register_unit({"degree", "deg", "Angle", Dimension(Dimension::Angle, 1), Scalar(3.141592653589793238462643383279502884L) / Scalar(180.0L), 0.0L, "Angle degree"});

    // ------------------------------------------------------------------------
    // 2. 长度单位 (Length)
    // ------------------------------------------------------------------------
    register_unit({"kilometer", "km", "Length", Dimension::Length(), 1000.0L, 0.0L, "1000 meters"});
    register_unit({"centimeter", "cm", "Length", Dimension::Length(), 0.01L, 0.0L, "1/100 meter"});
    register_unit({"millimeter", "mm", "Length", Dimension::Length(), 0.001L, 0.0L, "1/1000 meter"});
    register_unit({"micrometer", "um", "Length", Dimension::Length(), 1e-6L, 0.0L, "Micrometer / micron"});
    register_unit({"nanometer", "nm", "Length", Dimension::Length(), 1e-9L, 0.0L, "Nanometer"});
    register_unit({"picometer", "pm", "Length", Dimension::Length(), 1e-12L, 0.0L, "Picometer"});
    register_unit({"angstrom", "angstrom", "Length", Dimension::Length(), 1e-10L, 0.0L, "Angstrom (10^-10 m)"});
    register_unit({"inch", "in", "Length", Dimension::Length(), 0.0254L, 0.0L, "Inch (2.54 cm)"});
    register_unit({"foot", "ft", "Length", Dimension::Length(), 0.3048L, 0.0L, "Foot (12 inches)"});
    register_unit({"yard", "yd", "Length", Dimension::Length(), 0.9144L, 0.0L, "Yard (3 feet)"});
    register_unit({"mile", "mi", "Length", Dimension::Length(), 1609.344L, 0.0L, "Statute mile"});
    register_unit({"nautical_mile", "nmi", "Length", Dimension::Length(), 1852.0L, 0.0L, "Nautical mile"});
    register_unit({"astronomical_unit", "au", "Length", Dimension::Length(), 149597870700.0L, 0.0L, "Astronomical Unit"});
    register_unit({"light_year", "ly", "Length", Dimension::Length(), 9460730472580800.0L, 0.0L, "Light year"});
    register_unit({"parsec", "pc", "Length", Dimension::Length(), 3.085677581491367e16L, 0.0L, "Parsec"});

    // ------------------------------------------------------------------------
    // 3. 质量单位 (Mass)
    // ------------------------------------------------------------------------
    register_unit({"gram", "g", "Mass", Dimension::Mass(), 0.001L, 0.0L, "Gram (1/1000 kg)"});
    register_unit({"milligram", "mg", "Mass", Dimension::Mass(), 1e-6L, 0.0L, "Milligram"});
    register_unit({"microgram", "ug", "Mass", Dimension::Mass(), 1e-9L, 0.0L, "Microgram"});
    register_unit({"tonne", "t", "Mass", Dimension::Mass(), 1000.0L, 0.0L, "Metric ton (1000 kg)"});
    register_unit({"pound", "lb", "Mass", Dimension::Mass(), 0.45359237L, 0.0L, "Avoirdupois pound"});
    register_unit({"ounce", "oz", "Mass", Dimension::Mass(), 0.028349523125L, 0.0L, "Ounce (1/16 lb)"});
    register_unit({"ton_us", "ton", "Mass", Dimension::Mass(), 907.18474L, 0.0L, "US short ton (2000 lb)"});
    register_unit({"slug", "slug", "Mass", Dimension::Mass(), 14.5939029L, 0.0L, "Slug mass unit"});

    // ------------------------------------------------------------------------
    // 4. 时间单位 (Time)
    // ------------------------------------------------------------------------
    register_unit({"millisecond", "ms", "Time", Dimension::Time(), 1e-3L, 0.0L, "Millisecond"});
    register_unit({"microsecond", "us", "Time", Dimension::Time(), 1e-6L, 0.0L, "Microsecond"});
    register_unit({"nanosecond", "ns", "Time", Dimension::Time(), 1e-9L, 0.0L, "Nanosecond"});
    register_unit({"minute", "min", "Time", Dimension::Time(), 60.0L, 0.0L, "Minute (60 s)"});
    register_unit({"hour", "h", "Time", Dimension::Time(), 3600.0L, 0.0L, "Hour (3600 s)"});
    register_unit({"hour_alias", "hr", "Time", Dimension::Time(), 3600.0L, 0.0L, "Hour alias"});
    register_unit({"day", "d", "Time", Dimension::Time(), 86400.0L, 0.0L, "Day (86400 s)"});
    register_unit({"day_alias", "day", "Time", Dimension::Time(), 86400.0L, 0.0L, "Day alias"});
    register_unit({"week", "week", "Time", Dimension::Time(), 604800.0L, 0.0L, "Week (7 days)"});
    register_unit({"year", "yr", "Time", Dimension::Time(), 31557600.0L, 0.0L, "Julian year (365.25 d)"});
    register_unit({"year_alias", "year", "Time", Dimension::Time(), 31557600.0L, 0.0L, "Julian year alias"});

    // ------------------------------------------------------------------------
    // 5. 频率与力学衍生单位
    // ------------------------------------------------------------------------
    register_unit({"hertz", "Hz", "Frequency", Dimension::Frequency(), 1.0L, 0.0L, "Hertz (1/s)"});
    register_unit({"kilohertz", "kHz", "Frequency", Dimension::Frequency(), 1e3L, 0.0L, "Kilohertz"});
    register_unit({"megahertz", "MHz", "Frequency", Dimension::Frequency(), 1e6L, 0.0L, "Megahertz"});
    register_unit({"gigahertz", "GHz", "Frequency", Dimension::Frequency(), 1e9L, 0.0L, "Gigahertz"});
    register_unit({"terahertz", "THz", "Frequency", Dimension::Frequency(), 1e12L, 0.0L, "Terahertz"});
    register_unit({"rpm", "rpm", "Frequency", Dimension::Frequency(), Scalar(1.0L) / Scalar(60.0L), 0.0L, "Revolutions per minute"});

    register_unit({"newton", "N", "Force", Dimension::Force(), 1.0L, 0.0L, "Newton (kg*m/s^2)"});
    register_unit({"kilonewton", "kN", "Force", Dimension::Force(), 1e3L, 0.0L, "Kilonewton"});
    register_unit({"meganewton", "MN", "Force", Dimension::Force(), 1e6L, 0.0L, "Meganewton"});
    register_unit({"dyne", "dyn", "Force", Dimension::Force(), 1e-5L, 0.0L, "Dyne (10^-5 N)"});
    register_unit({"pound_force", "lbf", "Force", Dimension::Force(), 4.4482216152605L, 0.0L, "Pound force"});
    register_unit({"kilogram_force", "kgf", "Force", Dimension::Force(), 9.80665L, 0.0L, "Kilogram force"});

    // ------------------------------------------------------------------------
    // 6. 压强与应力 (Pressure)
    // ------------------------------------------------------------------------
    register_unit({"pascal", "Pa", "Pressure", Dimension::Pressure(), 1.0L, 0.0L, "Pascal (N/m^2)"});
    register_unit({"kilopascal", "kPa", "Pressure", Dimension::Pressure(), 1e3L, 0.0L, "Kilopascal"});
    register_unit({"megapascal", "MPa", "Pressure", Dimension::Pressure(), 1e6L, 0.0L, "Megapascal"});
    register_unit({"gigapascal", "GPa", "Pressure", Dimension::Pressure(), 1e9L, 0.0L, "Gigapascal"});
    register_unit({"bar", "bar", "Pressure", Dimension::Pressure(), 1e5L, 0.0L, "Bar (100 kPa)"});
    register_unit({"millibar", "mbar", "Pressure", Dimension::Pressure(), 100.0L, 0.0L, "Millibar (1 hPa)"});
    register_unit({"atmosphere", "atm", "Pressure", Dimension::Pressure(), 101325.0L, 0.0L, "Standard atmosphere"});
    register_unit({"psi", "psi", "Pressure", Dimension::Pressure(), 6894.757293168L, 0.0L, "Pounds per square inch"});
    register_unit({"torr", "torr", "Pressure", Dimension::Pressure(), 101325.0L / 760.0L, 0.0L, "Torr (1/760 atm)"});
    register_unit({"mmHg", "mmHg", "Pressure", Dimension::Pressure(), 133.322387415L, 0.0L, "Millimeter of mercury"});

    // ------------------------------------------------------------------------
    // 7. 能量与功率 (Energy & Power)
    // ------------------------------------------------------------------------
    register_unit({"joule", "J", "Energy", Dimension::Energy(), 1.0L, 0.0L, "Joule (N*m)"});
    register_unit({"kilojoule", "kJ", "Energy", Dimension::Energy(), 1e3L, 0.0L, "Kilojoule"});
    register_unit({"megajoule", "MJ", "Energy", Dimension::Energy(), 1e6L, 0.0L, "Megajoule"});
    register_unit({"gigajoule", "GJ", "Energy", Dimension::Energy(), 1e9L, 0.0L, "Gigajoule"});
    register_unit({"calorie", "cal", "Energy", Dimension::Energy(), 4.184L, 0.0L, "Thermochemical calorie"});
    register_unit({"kilocalorie", "kcal", "Energy", Dimension::Energy(), 4184.0L, 0.0L, "Kilocalorie (food calorie)"});
    register_unit({"watt_hour", "Wh", "Energy", Dimension::Energy(), 3600.0L, 0.0L, "Watt hour (3600 J)"});
    register_unit({"kilowatt_hour", "kWh", "Energy", Dimension::Energy(), 3.6e6L, 0.0L, "Kilowatt hour"});
    register_unit({"megawatt_hour", "MWh", "Energy", Dimension::Energy(), 3.6e9L, 0.0L, "Megawatt hour"});
    register_unit({"electronvolt", "eV", "Energy", Dimension::Energy(), 1.602176634e-19L, 0.0L, "Electronvolt"});
    register_unit({"kiloelectronvolt", "keV", "Energy", Dimension::Energy(), 1.602176634e-16L, 0.0L, "keV"});
    register_unit({"megaelectronvolt", "MeV", "Energy", Dimension::Energy(), 1.602176634e-13L, 0.0L, "MeV"});
    register_unit({"gigaelectronvolt", "GeV", "Energy", Dimension::Energy(), 1.602176634e-10L, 0.0L, "GeV"});
    register_unit({"teraelectronvolt", "TeV", "Energy", Dimension::Energy(), 1.602176634e-7L, 0.0L, "TeV"});
    register_unit({"btu", "BTU", "Energy", Dimension::Energy(), 1055.05585262L, 0.0L, "British Thermal Unit"});
    register_unit({"erg", "erg", "Energy", Dimension::Energy(), 1e-7L, 0.0L, "Erg (10^-7 J)"});

    register_unit({"watt", "W", "Power", Dimension::Power(), 1.0L, 0.0L, "Watt (J/s)"});
    register_unit({"milliwatt", "mW", "Power", Dimension::Power(), 1e-3L, 0.0L, "Milliwatt"});
    register_unit({"kilowatt", "kW", "Power", Dimension::Power(), 1e3L, 0.0L, "Kilowatt"});
    register_unit({"megawatt", "MW", "Power", Dimension::Power(), 1e6L, 0.0L, "Megawatt"});
    register_unit({"gigawatt", "GW", "Power", Dimension::Power(), 1e9L, 0.0L, "Gigawatt"});
    register_unit({"horsepower", "hp", "Power", Dimension::Power(), 745.69987158227L, 0.0L, "Mechanical horsepower (550 ft*lbf/s)"});

    // ------------------------------------------------------------------------
    // 8. 电磁学衍生单位 (Electromagnetism)
    // ------------------------------------------------------------------------
    register_unit({"coulomb", "C", "Charge", Dimension::Charge(), 1.0L, 0.0L, "Coulomb (A*s)"});
    register_unit({"milliampere_hour", "mAh", "Charge", Dimension::Charge(), 3.6L, 0.0L, "Milliampere hour"});
    register_unit({"ampere_hour", "Ah", "Charge", Dimension::Charge(), 3600.0L, 0.0L, "Ampere hour"});

    register_unit({"volt", "V", "Voltage", Dimension::Voltage(), 1.0L, 0.0L, "Volt (W/A)"});
    register_unit({"millivolt", "mV", "Voltage", Dimension::Voltage(), 1e-3L, 0.0L, "Millivolt"});
    register_unit({"kilovolt", "kV", "Voltage", Dimension::Voltage(), 1e3L, 0.0L, "Kilovolt"});

    register_unit({"ohm", "ohm", "Resistance", Dimension::Resistance(), 1.0L, 0.0L, "Ohm (V/A)"});
    register_unit({"milliohm", "mohm", "Resistance", Dimension::Resistance(), 1e-3L, 0.0L, "Milliohm"});
    register_unit({"kiloohm", "kohm", "Resistance", Dimension::Resistance(), 1e3L, 0.0L, "Kiloohm"});
    register_unit({"megaohm", "Mohm", "Resistance", Dimension::Resistance(), 1e6L, 0.0L, "Megaohm"});

    register_unit({"farad", "F", "Capacitance", Dimension::Capacitance(), 1.0L, 0.0L, "Farad (C/V)"});
    register_unit({"microfarad", "uF", "Capacitance", Dimension::Capacitance(), 1e-6L, 0.0L, "Microfarad"});
    register_unit({"nanofarad", "nF", "Capacitance", Dimension::Capacitance(), 1e-9L, 0.0L, "Nanofarad"});
    register_unit({"picofarad", "pF", "Capacitance", Dimension::Capacitance(), 1e-12L, 0.0L, "Picofarad"});

    register_unit({"henry", "H", "Inductance", Dimension(std::array<int8_t, 8>{2, 1, -2, -2, 0, 0, 0, 0}), 1.0L, 0.0L, "Henry (V*s/A)"});
    register_unit({"millihenry", "mH", "Inductance", Dimension(std::array<int8_t, 8>{2, 1, -2, -2, 0, 0, 0, 0}), 1e-3L, 0.0L, "Millihenry"});
    register_unit({"microhenry", "uH", "Inductance", Dimension(std::array<int8_t, 8>{2, 1, -2, -2, 0, 0, 0, 0}), 1e-6L, 0.0L, "Microhenry"});

    register_unit({"weber", "Wb", "MagneticFlux", Dimension(std::array<int8_t, 8>{2, 1, -2, -1, 0, 0, 0, 0}), 1.0L, 0.0L, "Weber (V*s)"});
    register_unit({"tesla", "T", "MagneticField", Dimension(std::array<int8_t, 8>{0, 1, -2, -1, 0, 0, 0, 0}), 1.0L, 0.0L, "Tesla (Wb/m^2)"});
    register_unit({"gauss", "G", "MagneticField", Dimension(std::array<int8_t, 8>{0, 1, -2, -1, 0, 0, 0, 0}), 1e-4L, 0.0L, "Gauss (10^-4 T)"});

    // ------------------------------------------------------------------------
    // 9. 面积与体积 (Area & Volume)
    // ------------------------------------------------------------------------
    register_unit({"square_meter", "m2", "Area", Dimension::Area(), 1.0L, 0.0L, "Square meter"});
    register_unit({"square_centimeter", "cm2", "Area", Dimension::Area(), 1e-4L, 0.0L, "Square centimeter"});
    register_unit({"square_millimeter", "mm2", "Area", Dimension::Area(), 1e-6L, 0.0L, "Square millimeter"});
    register_unit({"square_kilometer", "km2", "Area", Dimension::Area(), 1e6L, 0.0L, "Square kilometer"});
    register_unit({"hectare", "ha", "Area", Dimension::Area(), 10000.0L, 0.0L, "Hectare (10000 m^2)"});
    register_unit({"acre", "acre", "Area", Dimension::Area(), 4046.8564224L, 0.0L, "Acre (43560 sq ft)"});
    register_unit({"square_foot", "sqft", "Area", Dimension::Area(), 0.09290304L, 0.0L, "Square foot"});

    register_unit({"cubic_meter", "m3", "Volume", Dimension::Volume(), 1.0L, 0.0L, "Cubic meter"});
    register_unit({"cubic_centimeter", "cm3", "Volume", Dimension::Volume(), 1e-6L, 0.0L, "Cubic centimeter (cc / mL)"});
    register_unit({"liter", "L", "Volume", Dimension::Volume(), 0.001L, 0.0L, "Liter (1 dm^3)"});
    register_unit({"milliliter", "mL", "Volume", Dimension::Volume(), 1e-6L, 0.0L, "Milliliter"});
    register_unit({"gallon_us", "gal", "Volume", Dimension::Volume(), 0.003785411784L, 0.0L, "US liquid gallon"});
    register_unit({"fluid_ounce", "fl_oz", "Volume", Dimension::Volume(), 2.95735295625e-5L, 0.0L, "US fluid ounce"});

    // ------------------------------------------------------------------------
    // 10. 速度与加速度 (Speed & Acceleration)
    // ------------------------------------------------------------------------
    register_unit({"meter_per_second", "m/s", "Speed", Dimension::Velocity(), 1.0L, 0.0L, "Meter per second"});
    register_unit({"kilometer_per_hour", "km/h", "Speed", Dimension::Velocity(), Scalar(1000.0L) / Scalar(3600.0L), 0.0L, "Kilometer per hour"});
    register_unit({"mile_per_hour", "mph", "Speed", Dimension::Velocity(), Scalar(1609.344L) / Scalar(3600.0L), 0.0L, "Miles per hour"});
    register_unit({"knot", "knot", "Speed", Dimension::Velocity(), Scalar(1852.0L) / Scalar(3600.0L), 0.0L, "Nautical miles per hour"});
    register_unit({"feet_per_second", "ft/s", "Speed", Dimension::Velocity(), 0.3048L, 0.0L, "Feet per second"});

    register_unit({"standard_gravity", "g_force", "Acceleration", Dimension::Acceleration(), 9.80665L, 0.0L, "Standard earth acceleration"});

    // ------------------------------------------------------------------------
    // 11. 温度单位 (带仿射偏移) (Temperature)
    // ------------------------------------------------------------------------
    register_unit({"celsius", "degC", "Temperature", Dimension::Temperature(), 1.0L, 273.15L, "Degree Celsius"});
    register_unit({"fahrenheit", "degF", "Temperature", Dimension::Temperature(), Scalar(5.0L) / Scalar(9.0L), Scalar(273.15L) - Scalar(32.0L) * Scalar(5.0L) / Scalar(9.0L), "Degree Fahrenheit"});
    register_unit({"rankine", "degR", "Temperature", Dimension::Temperature(), Scalar(5.0L) / Scalar(9.0L), 0.0L, "Rankine"});

    // 映射代表性量纲到经典具名单位
    dimension_canonical_unit_[Dimension::Force().to_string()] = "N";
    dimension_canonical_unit_[Dimension::Energy().to_string()] = "J";
    dimension_canonical_unit_[Dimension::Power().to_string()] = "W";
    dimension_canonical_unit_[Dimension::Pressure().to_string()] = "Pa";
    dimension_canonical_unit_[Dimension::Voltage().to_string()] = "V";
    dimension_canonical_unit_[Dimension::Charge().to_string()] = "C";
    dimension_canonical_unit_[Dimension::Resistance().to_string()] = "ohm";
    dimension_canonical_unit_[Dimension::Capacitance().to_string()] = "F";
    dimension_canonical_unit_[Dimension::Frequency().to_string()] = "Hz";
    dimension_canonical_unit_[Dimension::Velocity().to_string()] = "m/s";
    dimension_canonical_unit_[Dimension::Acceleration().to_string()] = "m/s^2";
    dimension_canonical_unit_[Dimension::Area().to_string()] = "m^2";
    dimension_canonical_unit_[Dimension::Volume().to_string()] = "m^3";
}

const UnitInfo* UnitRegistry::find_unit(const std::string& symbol) const {
    auto it = symbol_map_.find(symbol);
    if (it != symbol_map_.end()) {
        return &units_list_[it->second];
    }
    return nullptr;
}

std::string UnitRegistry::simplify_dimension(const Dimension& dim) const {
    auto it = dimension_canonical_unit_.find(dim.to_string());
    if (it != dimension_canonical_unit_.end()) {
        return it->second;
    }
    return dim.to_string();
}

namespace {

struct Token {
    enum Type { Symbol, Mul, Div, Pow, OpenParen, CloseParen, Number, End } type;
    std::string text;
    int int_val = 0;
};

std::vector<Token> tokenize_unit(const std::string& str) {
    std::vector<Token> tokens;
    std::size_t i = 0;
    while (i < str.size()) {
        char c = str[i];
        if (std::isspace(static_cast<unsigned char>(c))) {
            ++i;
            continue;
        }
        if (c == '*' || c == '.') {
            tokens.push_back({Token::Mul, "*", 0});
            ++i;
        } else if (c == '/') {
            tokens.push_back({Token::Div, "/", 0});
            ++i;
        } else if (c == '^') {
            tokens.push_back({Token::Pow, "^", 0});
            ++i;
        } else if (c == '(') {
            tokens.push_back({Token::OpenParen, "(", 0});
            ++i;
        } else if (c == ')') {
            tokens.push_back({Token::CloseParen, ")", 0});
            ++i;
        } else if (std::isdigit(static_cast<unsigned char>(c)) || (c == '-' && i + 1 < str.size() && std::isdigit(static_cast<unsigned char>(str[i+1])))) {
            std::size_t start = i;
            if (c == '-') ++i;
            while (i < str.size() && std::isdigit(static_cast<unsigned char>(str[i]))) {
                ++i;
            }
            std::string num_str = str.substr(start, i - start);
            tokens.push_back({Token::Number, num_str, std::stoi(num_str)});
        } else if (std::isalpha(static_cast<unsigned char>(c)) || c == '_' || c == '%') {
            std::size_t start = i;
            while (i < str.size() && (std::isalnum(static_cast<unsigned char>(str[i])) || str[i] == '_' || str[i] == '%')) {
                ++i;
            }
            tokens.push_back({Token::Symbol, str.substr(start, i - start), 0});
        } else {
            ++i;
        }
    }
    tokens.push_back({Token::End, "", 0});
    return tokens;
}

} // namespace

bool UnitRegistry::parse_unit_expression(const std::string& expr,
                                        Dimension* out_dimension,
                                        Scalar* out_scale_to_si,
                                        Scalar* out_offset_to_si) const {
    std::string trimmed = trim_copy(expr);
    if (trimmed.empty()) {
        if (out_dimension) *out_dimension = Dimension::Dimensionless();
        if (out_scale_to_si) *out_scale_to_si = 1.0L;
        if (out_offset_to_si) *out_offset_to_si = 0.0L;
        return true;
    }

    // 单一已知单位快速查找
    if (const auto* direct = find_unit(trimmed)) {
        if (out_dimension) *out_dimension = direct->dimension;
        if (out_scale_to_si) *out_scale_to_si = direct->scale_to_si;
        if (out_offset_to_si) *out_offset_to_si = direct->offset_to_si;
        return true;
    }

    if (out_offset_to_si) *out_offset_to_si = 0.0L;

    // 解析复合表达式
    auto tokens = tokenize_unit(trimmed);
    std::size_t pos = 0;

    Dimension total_dim;
    Scalar total_scale = 1.0L;

    // 辅助解析单一项：symbol 或 (expr) 附带可选 ^n
    auto parse_factor = [&](auto& self, Dimension& f_dim, Scalar& f_scale) -> bool {
        if (pos >= tokens.size()) return false;
        if (tokens[pos].type == Token::OpenParen) {
            ++pos;
            Dimension p_dim;
            Scalar p_scale = 1.0L;
            // 简单括号内乘除
            bool in_div = false;
            while (pos < tokens.size() && tokens[pos].type != Token::CloseParen && tokens[pos].type != Token::End) {
                if (tokens[pos].type == Token::Mul) {
                    in_div = false; ++pos; continue;
                }
                if (tokens[pos].type == Token::Div) {
                    in_div = true; ++pos; continue;
                }
                Dimension item_dim;
                Scalar item_scale = 1.0L;
                if (!self(self, item_dim, item_scale)) return false;
                if (!in_div) {
                    p_dim = p_dim * item_dim;
                    p_scale = p_scale * item_scale;
                } else {
                    p_dim = p_dim / item_dim;
                    p_scale = p_scale / item_scale;
                }
            }
            if (pos < tokens.size() && tokens[pos].type == Token::CloseParen) {
                ++pos;
            }
            f_dim = p_dim;
            f_scale = p_scale;
        } else if (tokens[pos].type == Token::Symbol) {
            std::string sym = tokens[pos].text;
            ++pos;
            const auto* u = find_unit(sym);
            if (!u) {
                return false;
            }
            f_dim = u->dimension;
            f_scale = u->scale_to_si;
        } else if (tokens[pos].type == Token::Number && tokens[pos].int_val == 1) {
            f_dim = Dimension::Dimensionless();
            f_scale = 1.0L;
            ++pos;
        } else {
            return false;
        }

        // 检查是否有 ^exponent
        if (pos < tokens.size() && tokens[pos].type == Token::Pow) {
            ++pos;
            if (pos >= tokens.size() || tokens[pos].type != Token::Number) return false;
            int exp = tokens[pos].int_val;
            ++pos;
            f_dim = f_dim.pow(exp);
            Scalar s = 1.0L;
            if (exp > 0) {
                for (int k = 0; k < exp; ++k) s = s * f_scale;
            } else if (exp < 0) {
                for (int k = 0; k < -exp; ++k) s = s / f_scale;
            }
            f_scale = s;
        }
        return true;
    };

    bool is_div = false;
    while (pos < tokens.size() && tokens[pos].type != Token::End) {
        if (tokens[pos].type == Token::Mul) {
            is_div = false;
            ++pos;
            continue;
        }
        if (tokens[pos].type == Token::Div) {
            is_div = true;
            ++pos;
            continue;
        }

        Dimension cur_dim;
        Scalar cur_scale = 1.0L;
        if (!parse_factor(parse_factor, cur_dim, cur_scale)) {
            return false;
        }

        if (!is_div) {
            total_dim = total_dim * cur_dim;
            total_scale = total_scale * cur_scale;
        } else {
            total_dim = total_dim / cur_dim;
            total_scale = total_scale / cur_scale;
        }
    }

    if (out_dimension) *out_dimension = total_dim;
    if (out_scale_to_si) *out_scale_to_si = total_scale;
    return true;
}

} // namespace units
