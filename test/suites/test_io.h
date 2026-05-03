// ============================================================================
// IO Module Tests
// ============================================================================

#ifndef TEST_IO_H
#define TEST_IO_H

#include "core/calculator.h"
#include <filesystem>
#include <fstream>
#include <iostream>

inline void test_io_file_operations(int& passed, int& failed) {
    Calculator calc;
    std::string output;

    // Test exists function
    try {
        calc.try_process_function_command("exists(\"/tmp/test_io_calc.txt\")", &output);
        if (output.find("0") != std::string::npos) {
            passed++;
        } else {
            std::cout << "FAIL: exists should return 0 for non-existent file\n";
            failed++;
        }
    } catch (const std::exception& e) {
        std::cout << "FAIL: exists threw: " << e.what() << "\n";
        failed++;
    }

    // Test open/write/close
    try {
        calc.try_process_function_command("fd = open(\"/tmp/test_io_calc.txt\", \"w\")", &output);
        calc.try_process_function_command("write(fd, \"Hello, World!\")", &output);
        calc.try_process_function_command("close(fd)", &output);

        if (std::filesystem::exists("/tmp/test_io_calc.txt")) {
            passed++;
        } else {
            std::cout << "FAIL: file was not created\n";
            failed++;
        }
    } catch (const std::exception& e) {
        std::cout << "FAIL: open/write/close threw: " << e.what() << "\n";
        failed++;
    }

    // Test read
    try {
        calc.try_process_function_command("fd = open(\"/tmp/test_io_calc.txt\", \"r\")", &output);
        calc.try_process_function_command("content = read(fd)", &output);
        calc.try_process_function_command("close(fd)", &output);

        // Check the content variable value
        calc.try_process_function_command("content", &output);
        if (output.find("Hello") != std::string::npos || output.find("World") != std::string::npos) {
            passed++;
        } else {
            std::cout << "FAIL: read content mismatch: " << output << "\n";
            failed++;
        }
    } catch (const std::exception& e) {
        std::cout << "FAIL: read threw: " << e.what() << "\n";
        failed++;
    }

    // Test delete
    try {
        calc.try_process_function_command("delete(\"/tmp/test_io_calc.txt\")", &output);
        if (!std::filesystem::exists("/tmp/test_io_calc.txt")) {
            passed++;
        } else {
            std::cout << "FAIL: file was not deleted\n";
            failed++;
        }
    } catch (const std::exception& e) {
        std::cout << "FAIL: delete threw: " << e.what() << "\n";
        failed++;
    }
}

inline void test_io_csv_operations(int& passed, int& failed) {
    Calculator calc;
    std::string output;

    // Test write_csv and read_csv
    try {
        calc.try_process_function_command("m = [1, 2, 3; 4, 5, 6]", &output);
        calc.try_process_function_command("write_csv(\"/tmp/test_matrix.csv\", m)", &output);

        if (std::filesystem::exists("/tmp/test_matrix.csv")) {
            passed++;
        } else {
            std::cout << "FAIL: CSV file was not created\n";
            failed++;
        }

        calc.try_process_function_command("loaded = read_csv(\"/tmp/test_matrix.csv\")", &output);

        // Verify matrix content
        calc.try_process_function_command("loaded", &output);
        if (output.find("1") != std::string::npos && output.find("6") != std::string::npos) {
            passed++;
        } else {
            std::cout << "FAIL: CSV read content mismatch: " << output << "\n";
            failed++;
        }

        std::filesystem::remove("/tmp/test_matrix.csv");
    } catch (const std::exception& e) {
        std::cout << "FAIL: CSV operations threw: " << e.what() << "\n";
        failed++;
    }
}

inline void test_io_json_operations(int& passed, int& failed) {
    Calculator calc;
    std::string output;

    // Test JSON with matrix
    try {
        calc.try_process_function_command("m = [1, 2; 3, 4]", &output);
        calc.try_process_function_command("write_json(\"/tmp/test.json\", m)", &output);

        if (std::filesystem::exists("/tmp/test.json")) {
            passed++;
        } else {
            std::cout << "FAIL: JSON file was not created\n";
            failed++;
        }

        calc.try_process_function_command("loaded = read_json(\"/tmp/test.json\")", &output);
        calc.try_process_function_command("loaded", &output);

        if (output.find("1") != std::string::npos && output.find("4") != std::string::npos) {
            passed++;
        } else {
            std::cout << "FAIL: JSON matrix read mismatch: " << output << "\n";
            failed++;
        }

        std::filesystem::remove("/tmp/test.json");
    } catch (const std::exception& e) {
        std::cout << "FAIL: JSON matrix operations threw: " << e.what() << "\n";
        failed++;
    }

    // Test JSON with dict
    try {
        calc.try_process_function_command("data = {\"name\": \"test\", \"value\": 42}", &output);
        calc.try_process_function_command("write_json(\"/tmp/test_dict.json\", data)", &output);

        calc.try_process_function_command("loaded = read_json(\"/tmp/test_dict.json\")", &output);
        calc.try_process_function_command("loaded", &output);

        if (output.find("name") != std::string::npos && output.find("test") != std::string::npos) {
            passed++;
        } else {
            std::cout << "FAIL: JSON dict read mismatch: " << output << "\n";
            failed++;
        }

        std::filesystem::remove("/tmp/test_dict.json");
    } catch (const std::exception& e) {
        std::cout << "FAIL: JSON dict operations threw: " << e.what() << "\n";
        failed++;
    }
}

inline void test_io_seek_tell(int& passed, int& failed) {
    Calculator calc;
    std::string output;

    try {
        // Create test file
        std::ofstream out("/tmp/test_seek.txt");
        out << "Hello, World!";
        out.close();

        calc.try_process_function_command("fd = open(\"/tmp/test_seek.txt\", \"r\")", &output);

        // Test tell
        calc.try_process_function_command("pos = tell(fd)", &output);
        if (output.find("0") != std::string::npos) {
            passed++;
        } else {
            std::cout << "FAIL: tell should return 0 at start\n";
            failed++;
        }

        // Test seek
        calc.try_process_function_command("seek(fd, 7)", &output);
        calc.try_process_function_command("pos = tell(fd)", &output);
        if (output.find("7") != std::string::npos) {
            passed++;
        } else {
            std::cout << "FAIL: tell should return 7 after seek\n";
            failed++;
        }

        // Test readline after seek
        calc.try_process_function_command("line = readline(fd)", &output);
        if (output.find("World") != std::string::npos) {
            passed++;
        } else {
            std::cout << "FAIL: readline after seek mismatch: " << output << "\n";
            failed++;
        }

        calc.try_process_function_command("close(fd)", &output);
        std::filesystem::remove("/tmp/test_seek.txt");
    } catch (const std::exception& e) {
        std::cout << "FAIL: seek/tell operations threw: " << e.what() << "\n";
        failed++;
    }
}

inline void run_io_tests(int& passed, int& failed) {
    test_io_file_operations(passed, failed);
    test_io_csv_operations(passed, failed);
    test_io_json_operations(passed, failed);
    test_io_seek_tell(passed, failed);
}

#endif // TEST_IO_H
