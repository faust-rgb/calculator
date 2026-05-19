// ============================================================================
// main.cpp - 计算器应用程序入口
// ============================================================================

#include "core/api/calculator.h"
#include "app/repl_application.h"
#include <iostream>
#include <string>
#include <vector>
#include <unistd.h>

int main(int argc, char* argv[]) {
    Calculator calculator;
    bool plain_mode = !isatty(STDIN_FILENO);
    std::string script_path;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--plain") {
            plain_mode = true;
        } else if (arg == "--version") {
            std::cout << "Calculator Version 3.0.0 (Clean Architecture)\n";
            return 0;
        } else if (arg.size() > 5 && arg.substr(arg.size() - 5) == ".calc") {
            script_path = arg;
        }
    }

    if (!script_path.empty()) {
        try {
            std::cout << calculator.execute_script_file(script_path, false) << "\n";
            return 0;
        } catch (const std::exception& ex) {
            std::cerr << "Error: " << ex.what() << "\n";
            return 1;
        }
    }

    ReplApplication app(calculator);
    return app.run(plain_mode);
}
