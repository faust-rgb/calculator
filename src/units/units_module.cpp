// ============================================================================
// units_module.cpp - 物理单位与量纲分析模块实现
// ============================================================================

#include "units/units_module.h"
#include "units/quantity.h"
#include "units/physical_constants.h"
#include "core/services/format_utils.h"
#include "core/services/string_utils.h"
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace units {

namespace {

Quantity parse_quantity_arg(const StoredValue& val) {
    if (val.is_string_type()) {
        Quantity q;
        if (Quantity::try_parse(val.get_string_value(), &q)) {
            return q;
        }
        throw std::runtime_error("invalid quantity string: " + val.get_string_value());
    }
    if (val.is_scalar_type()) {
        return Quantity(val.as_scalar(), Dimension::Dimensionless(), 1.0L, "");
    }
    throw std::runtime_error("expected quantity string or scalar number");
}

} // namespace

std::map<std::string, std::function<StoredValue(const std::vector<StoredValue>&)>> UnitsModule::get_functions_map() const {
    std::map<std::string, std::function<StoredValue(const std::vector<StoredValue>&)>> map;

    // qty(value, "unit") 或 qty("50 km/h")
    auto qty_fn = [](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.empty() || args.size() > 2) {
            throw std::runtime_error("qty requires 1 or 2 arguments: qty(val, \"unit\") or qty(\"50 km/h\")");
        }
        if (args.size() == 1) {
            Quantity q = parse_quantity_arg(args[0]);
            return StoredValue(q.to_string());
        }
        Scalar val = args[0].as_scalar();
        std::string u = args[1].as_string();
        Quantity q(val, u);
        return StoredValue(q.to_string());
    };
    map["qty"] = qty_fn;
    map["unit"] = qty_fn;

    // convert_unit(qty, "to_unit") 或 convert_unit(val, "from_unit", "to_unit")
    auto convert_fn = [](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.size() == 2) {
            Quantity q = parse_quantity_arg(args[0]);
            std::string to_u = args[1].as_string();
            Quantity conv = q.convert_to(to_u);
            return StoredValue(conv.to_string());
        }
        if (args.size() == 3) {
            Scalar val = args[0].as_scalar();
            std::string from_u = args[1].as_string();
            std::string to_u = args[2].as_string();
            Quantity q(val, from_u);
            Quantity conv = q.convert_to(to_u);
            return StoredValue(conv.to_string());
        }
        throw std::runtime_error("convert_unit requires 2 or 3 arguments: convert_unit(qty, \"to_u\") or convert_unit(val, \"from_u\", \"to_u\")");
    };
    map["convert_unit"] = convert_fn;
    map["to_unit"] = convert_fn;
    map["in_unit"] = convert_fn;

    // unit_val("50 km/h") -> 50
    map["unit_val"] = [](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.size() != 1) throw std::runtime_error("unit_val requires 1 argument");
        Quantity q = parse_quantity_arg(args[0]);
        return StoredValue(q.value);
    };

    // unit_dim("50 km/h") -> "[L T^-1]"
    map["unit_dim"] = [](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.size() != 1) throw std::runtime_error("unit_dim requires 1 argument");
        if (args[0].is_string_type()) {
            std::string str = args[0].get_string_value();
            Quantity q;
            if (Quantity::try_parse(str, &q)) {
                return StoredValue(q.dimension.to_string());
            }
            Dimension dim;
            Scalar scale;
            if (UnitRegistry::instance().parse_unit_expression(str, &dim, &scale)) {
                return StoredValue(dim.to_string());
            }
        }
        Quantity q = parse_quantity_arg(args[0]);
        return StoredValue(q.dimension.to_string());
    };

    // unit_simplify("1000 g*m/s^2") -> "1 N"
    map["unit_simplify"] = [](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.size() != 1) throw std::runtime_error("unit_simplify requires 1 argument");
        Quantity q = parse_quantity_arg(args[0]);
        std::string sym = UnitRegistry::instance().simplify_dimension(q.dimension);
        Quantity simplified = q.convert_to(sym);
        return StoredValue(simplified.to_string());
    };

    // 算术运算: unit_add, unit_sub, unit_mul, unit_div, unit_pow
    map["unit_add"] = [](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.size() != 2) throw std::runtime_error("unit_add requires 2 arguments");
        Quantity q1 = parse_quantity_arg(args[0]);
        Quantity q2 = parse_quantity_arg(args[1]);
        return StoredValue((q1 + q2).to_string());
    };

    map["unit_sub"] = [](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.size() != 2) throw std::runtime_error("unit_sub requires 2 arguments");
        Quantity q1 = parse_quantity_arg(args[0]);
        Quantity q2 = parse_quantity_arg(args[1]);
        return StoredValue((q1 - q2).to_string());
    };

    map["unit_mul"] = [](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.size() != 2) throw std::runtime_error("unit_mul requires 2 arguments");
        if (args[0].is_scalar_type() && args[1].is_string_type()) {
            Quantity q2 = parse_quantity_arg(args[1]);
            return StoredValue((q2 * args[0].as_scalar()).to_string());
        }
        if (args[0].is_string_type() && args[1].is_scalar_type()) {
            Quantity q1 = parse_quantity_arg(args[0]);
            return StoredValue((q1 * args[1].as_scalar()).to_string());
        }
        Quantity q1 = parse_quantity_arg(args[0]);
        Quantity q2 = parse_quantity_arg(args[1]);
        return StoredValue((q1 * q2).to_string());
    };

    map["unit_div"] = [](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.size() != 2) throw std::runtime_error("unit_div requires 2 arguments");
        if (args[0].is_string_type() && args[1].is_scalar_type()) {
            Quantity q1 = parse_quantity_arg(args[0]);
            return StoredValue((q1 / args[1].as_scalar()).to_string());
        }
        Quantity q1 = parse_quantity_arg(args[0]);
        Quantity q2 = parse_quantity_arg(args[1]);
        return StoredValue((q1 / q2).to_string());
    };

    map["unit_pow"] = [](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.size() != 2) throw std::runtime_error("unit_pow requires 2 arguments: unit_pow(qty, exp)");
        Quantity q = parse_quantity_arg(args[0]);
        int exp = static_cast<int>(args[1].as_scalar());
        return StoredValue(q.pow(exp).to_string());
    };

    // 物理常数
    for (const auto& c : PhysicalConstants::get_all()) {
        std::string q_str = c.quantity.to_string();
        map[c.name] = [q_str](const std::vector<StoredValue>&) -> StoredValue {
            return StoredValue(q_str);
        };
    }

    return map;
}

std::vector<std::string> UnitsModule::get_function_names() const {
    std::vector<std::string> names = {
        "qty", "unit", "convert_unit", "to_unit", "in_unit", "unit_val", "unit_dim", "unit_simplify",
        "unit_add", "unit_sub", "unit_mul", "unit_div", "unit_pow"
    };
    for (const auto& c : PhysicalConstants::get_all()) {
        names.push_back(c.name);
    }
    return names;
}

std::vector<std::string> UnitsModule::get_commands() const {
    return {":units", ":constants", ":convert"};
}

std::string UnitsModule::get_help_snippet(const std::string& topic) const {
    if (topic == "commands" || topic == "units") {
        return "Physical Units & Constants:\n"
               "  :units              List all supported physical units\n"
               "  :constants          List fundamental physical constants\n"
               "  :convert <val> <from_u> [to] <to_u>  Convert between units\n"
               "Functions:\n"
               "  qty(val, \"u\"), convert_unit(qty, \"to_u\"), unit_dim(qty), unit_simplify(qty)\n"
               "  unit_add(q1, q2), unit_sub(q1, q2), unit_mul(q1, q2), unit_div(q1, q2)\n"
               "Constants: c_light, h_planck, hbar, G_gravity, g_earth, k_boltzmann, e_charge, R_gas";
    }
    return "";
}

std::string UnitsModule::execute_args_view(std::string_view command,
                                           const std::vector<std::string_view>& args,
                                           ModuleServices&) {
    if (command == ":units") {
        const auto& all = UnitRegistry::instance().get_all_units();
        std::ostringstream oss;
        oss << "Supported Physical Units (" << all.size() << " units):\n";
        oss << "------------------------------------------------------------\n";
        oss << std::left << std::setw(12) << "Symbol"
            << std::setw(20) << "Name"
            << std::setw(16) << "Category"
            << "Description\n";
        oss << "------------------------------------------------------------\n";
        for (const auto& u : all) {
            oss << std::left << std::setw(12) << u.symbol
                << std::setw(20) << u.name
                << std::setw(16) << u.category
                << u.description << "\n";
        }
        return oss.str();
    }

    if (command == ":constants") {
        const auto& all = PhysicalConstants::get_all();
        std::ostringstream oss;
        oss << "Fundamental Physical Constants:\n";
        oss << "--------------------------------------------------------------------------------\n";
        oss << std::left << std::setw(16) << "Name"
            << std::setw(10) << "Symbol"
            << std::setw(30) << "Value & Unit"
            << "Description\n";
        oss << "--------------------------------------------------------------------------------\n";
        for (const auto& c : all) {
            oss << std::left << std::setw(16) << c.name
                << std::setw(10) << c.symbol
                << std::setw(30) << c.quantity.to_string()
                << c.description << "\n";
        }
        return oss.str();
    }

    if (command == ":convert") {
        std::string raw_line;
        for (auto a : args) {
            if (!raw_line.empty() && a != "/" && a != "^" && raw_line.back() != '/' && raw_line.back() != '^') {
                raw_line += " ";
            }
            raw_line += std::string(a);
        }

        std::string from_part, to_part;
        std::size_t to_pos = raw_line.find(" to ");
        if (to_pos != std::string::npos) {
            from_part = trim_copy(raw_line.substr(0, to_pos));
            to_part = trim_copy(raw_line.substr(to_pos + 4));
        } else if ((to_pos = raw_line.find(" in ")) != std::string::npos) {
            from_part = trim_copy(raw_line.substr(0, to_pos));
            to_part = trim_copy(raw_line.substr(to_pos + 4));
        } else if ((to_pos = raw_line.find(" -> ")) != std::string::npos) {
            from_part = trim_copy(raw_line.substr(0, to_pos));
            to_part = trim_copy(raw_line.substr(to_pos + 4));
        } else {
            std::size_t last_space = raw_line.rfind(' ');
            if (last_space != std::string::npos) {
                from_part = trim_copy(raw_line.substr(0, last_space));
                to_part = trim_copy(raw_line.substr(last_space + 1));
            }
        }

        if (from_part.empty() || to_part.empty()) {
            return "Usage: :convert <value> <from_unit> to <to_unit>  (e.g., :convert 100 km/h to m/s)";
        }
        try {
            Quantity q;
            if (!Quantity::try_parse(from_part, &q)) {
                return "Failed to parse quantity: " + from_part;
            }
            Quantity res = q.convert_to(to_part);
            return res.to_string();
        } catch (const std::exception& ex) {
            return std::string("Conversion error: ") + ex.what();
        }
    }

    return "";
}

} // namespace units
