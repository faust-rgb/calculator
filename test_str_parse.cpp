#include <iostream>
#include <string_view>
#include "src/core/command_parser.h"

int main() {
    CommandParser p("\"abc\" + \"def\"");
    CommandASTNode n = p.parse();
    if (n.kind == CommandKind::kStringLiteral) {
        std::cout << "StringLiteral: " << n.get_string_literal() << std::endl;
    } else if (n.kind == CommandKind::kExpression) {
        std::cout << "Expression" << std::endl;
    }
    return 0;
}
