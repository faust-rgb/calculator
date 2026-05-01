import re

cpp_file = "src/analysis/calculator_integration.cpp"
with open(cpp_file, "r") as f:
    content = f.read()

# Replace the coord system parsing part
old_parse = """    std::string coord_system = "rect";
    if (arguments.size() >= 2) {
        std::string last = utils::trim_copy(arguments.back());
        if (last == "\\"polar\\"" || last == "polar" || last == "\\"cyl\\"" || last == "cyl" || last == "\\"sph\\"" || last == "sph" || last == "\\"rect\\"" || last == "rect") {
            coord_system = last; if (coord_system.front() == '"') coord_system = coord_system.substr(1, coord_system.size() - 2);
            arguments.pop_back();
        }
    }"""

new_parse = """    std::string coord_system = "rect";
    std::string method = "simpson";
    double tol = 1e-6;
    bool has_options = true;
    while (arguments.size() >= 2 && has_options) {
        has_options = false;
        std::string last = utils::trim_copy(arguments.back());
        if (last.front() == '"' && last.back() == '"') last = last.substr(1, last.size() - 2);
        
        if (last == "polar" || last == "cyl" || last == "sph" || last == "rect") {
            coord_system = last;
            arguments.pop_back();
            has_options = true;
        } else if (last.find("method=") == 0) {
            method = last.substr(7);
            if (method.front() == '"' && method.back() == '"') method = method.substr(1, method.size() - 2);
            arguments.pop_back();
            has_options = true;
        } else if (last.find("tol=") == 0) {
            try { tol = ctx.parse_decimal(last.substr(4)); } catch (...) {}
            arguments.pop_back();
            has_options = true;
        } else if (last == "adaptive" || last == "monte_carlo" || last == "sparse_grid" || last == "simpson") {
            method = last;
            arguments.pop_back();
            has_options = true;
        }
    }"""

content = content.replace(old_parse, new_parse)

# We need to change double_integral and triple_integral to use method and tol.
# But it's easier to implement the logic inside double_integral itself if we pass method and tol down.
# Let's change the signatures of double_integral and triple_integral.

