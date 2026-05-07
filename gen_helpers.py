import re

with open("src/symbolic/integral/symbolic_expression_integral_helpers.cpp", "r") as f:
    lines = f.readlines()

header_lines = [
    "#ifndef SYMBOLIC_EXPRESSION_INTEGRAL_HELPERS_H",
    "#define SYMBOLIC_EXPRESSION_INTEGRAL_HELPERS_H",
    "#include \"symbolic/symbolic_expression_internal.h\"",
    "namespace symbolic_expression_internal {",
    ""
]

inside_anon_ns = False
for line in lines:
    if line.startswith("namespace {"):
        inside_anon_ns = True
        continue
        
    if inside_anon_ns:
        # Match function signatures like: bool try_...(...) {
        m = re.match(r'^([a-zA-Z0-9_:\<\>]+(?:[\s\*&]+[a-zA-Z0-9_:\<\>]+)*)\s+([a-zA-Z0-9_]+)\s*\((.*?)\)\s*\{', line.strip())
        if m:
            ret_type = m.group(1)
            name = m.group(2)
            args = m.group(3)
            # Filter out things that look like if/for/while
            if name not in ["if", "for", "while", "switch", "catch"]:
                header_lines.append(f"{ret_type} {name}({args});")

header_lines.append("")
header_lines.append("} // namespace symbolic_expression_internal")
header_lines.append("#endif")

with open("src/symbolic/integral/symbolic_expression_integral_helpers.h", "w") as f:
    f.write("\n".join(header_lines))
