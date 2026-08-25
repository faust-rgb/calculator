// ============================================================================
// IO Module Tests Implementation
// ============================================================================

#include "test_io.h"
#include "test_helpers.h"
#include "core/api/calculator.h"
#include <filesystem>
#include <fstream>
#include <iostream>

static void test_io_file_operations(int& passed, int& failed) {
    Calculator calc;
    std::string output;
    test_helpers::TempFileGuard guard(test_helpers::make_test_path("test_io_calc.txt"));
    const std::string path_str = guard.path.string();

    // Test exists function on non-existent file
    try {
        calc.try_process_function_command("exists(\"" + path_str + "\")", &output);
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
        calc.try_process_function_command("fd = open(\"" + path_str + "\", \"w\")", &output);
        calc.try_process_function_command("write(fd, \"Hello, World!\")", &output);
        calc.try_process_function_command("close(fd)", &output);

        if (std::filesystem::exists(guard.path)) {
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
        calc.try_process_function_command("fd = open(\"" + path_str + "\", \"r\")", &output);
        calc.try_process_function_command("content = read(fd)", &output);
        calc.try_process_function_command("close(fd)", &output);

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
        calc.try_process_function_command("delete(\"" + path_str + "\")", &output);
        if (!std::filesystem::exists(guard.path)) {
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

static void test_io_csv_operations(int& passed, int& failed) {
    Calculator calc;
    std::string output;
    test_helpers::TempFileGuard guard(test_helpers::make_test_path("test_matrix.csv"));
    const std::string path_str = guard.path.string();

    try {
        calc.try_process_function_command("m = [1, 2, 3; 4, 5, 6]", &output);
        calc.try_process_function_command("write_csv(\"" + path_str + "\", m)", &output);

        if (std::filesystem::exists(guard.path)) {
            passed++;
        } else {
            std::cout << "FAIL: CSV file was not created\n";
            failed++;
        }

        calc.try_process_function_command("loaded = read_csv(\"" + path_str + "\")", &output);

        calc.try_process_function_command("loaded", &output);
        if (output.find("1") != std::string::npos && output.find("6") != std::string::npos) {
            passed++;
        } else {
            std::cout << "FAIL: CSV read content mismatch: " << output << "\n";
            failed++;
        }
    } catch (const std::exception& e) {
        std::cout << "FAIL: CSV operations threw: " << e.what() << "\n";
        failed++;
    }
}

static void test_io_json_operations(int& passed, int& failed) {
    Calculator calc;
    std::string output;
    test_helpers::TempFileGuard guard_mat(test_helpers::make_test_path("test_matrix.json"));
    test_helpers::TempFileGuard guard_dict(test_helpers::make_test_path("test_dict.json"));

    // Test JSON with matrix
    try {
        calc.try_process_function_command("m = [1, 2; 3, 4]", &output);
        calc.try_process_function_command("write_json(\"" + guard_mat.path.string() + "\", m)", &output);

        if (std::filesystem::exists(guard_mat.path)) {
            passed++;
        } else {
            std::cout << "FAIL: JSON file was not created\n";
            failed++;
        }

        calc.try_process_function_command("loaded = read_json(\"" + guard_mat.path.string() + "\")", &output);
        calc.try_process_function_command("loaded", &output);

        if (output.find("1") != std::string::npos && output.find("4") != std::string::npos) {
            passed++;
        } else {
            std::cout << "FAIL: JSON matrix read mismatch: " << output << "\n";
            failed++;
        }
    } catch (const std::exception& e) {
        std::cout << "FAIL: JSON matrix operations threw: " << e.what() << "\n";
        failed++;
    }

    // Test JSON with dict
    try {
        calc.try_process_function_command("data = {\"name\": \"test\", \"value\": 42}", &output);
        calc.try_process_function_command("write_json(\"" + guard_dict.path.string() + "\", data)", &output);

        calc.try_process_function_command("loaded = read_json(\"" + guard_dict.path.string() + "\")", &output);
        calc.try_process_function_command("loaded", &output);

        if (output.find("name") != std::string::npos && output.find("test") != std::string::npos) {
            passed++;
        } else {
            std::cout << "FAIL: JSON dict read mismatch: " << output << "\n";
            failed++;
        }
    } catch (const std::exception& e) {
        std::cout << "FAIL: JSON dict operations threw: " << e.what() << "\n";
        failed++;
    }
}

static void test_io_seek_tell(int& passed, int& failed) {
    Calculator calc;
    std::string output;
    test_helpers::TempFileGuard guard(test_helpers::make_test_path("test_seek.txt"));
    const std::string path_str = guard.path.string();

    try {
        std::ofstream out(guard.path);
        out << "Hello, World!";
        out.close();

        calc.try_process_function_command("fd = open(\"" + path_str + "\", \"r\")", &output);

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
    } catch (const std::exception& e) {
        std::cout << "FAIL: seek/tell operations threw: " << e.what() << "\n";
        failed++;
    }
}

static void test_io_advanced_features(int& passed, int& failed) {
    Calculator calc;
    std::string output;
    test_helpers::TempFileGuard guard_txt(test_helpers::make_test_path("test_adv_io.txt"));
    test_helpers::TempFileGuard guard_csv(test_helpers::make_test_path("test_adv_csv.csv"));
    test_helpers::TempFileGuard guard_json(test_helpers::make_test_path("test_adv_json.json"));

    // 1. Test read with count & seek with whence
    try {
        calc.try_process_function_command("fd = open(\"" + guard_txt.path.string() + "\", \"w+\")", &output);
        calc.try_process_function_command("write(fd, \"0123456789ABCDEF\")", &output);
        calc.try_process_function_command("seek(fd, 4, 0)", &output); // whence 0 = beg
        calc.try_process_function_command("chunk = read(fd, 4)", &output);
        calc.try_process_function_command("chunk", &output);
        if (output.find("4567") != std::string::npos) {
            passed++;
        } else {
            std::cout << "FAIL: read with count mismatch: " << output << "\n";
            failed++;
        }

        // seek from end
        calc.try_process_function_command("seek(fd, -4, 2)", &output); // whence 2 = end
        calc.try_process_function_command("tail = read(fd, 4)", &output);
        calc.try_process_function_command("tail", &output);
        if (output.find("CDEF") != std::string::npos) {
            passed++;
        } else {
            std::cout << "FAIL: seek end mismatch: " << output << "\n";
            failed++;
        }

        calc.try_process_function_command("close(fd)", &output);
    } catch (const std::exception& e) {
        std::cout << "FAIL: advanced read/seek threw: " << e.what() << "\n";
        failed++;
    }

    // 2. Test quoted CSV parsing
    try {
        std::ofstream csv_out(guard_csv.path);
        csv_out << "\"10.5\", \"20.25\", 30\n\"40\", 50, \"60.75\"\n";
        csv_out.close();

        calc.try_process_function_command("mat = read_csv(\"" + guard_csv.path.string() + "\")", &output);
        calc.try_process_function_command("mat", &output);
        if (output.find("10.5") != std::string::npos && output.find("60.75") != std::string::npos) {
            passed++;
        } else {
            std::cout << "FAIL: quoted CSV parse mismatch: " << output << "\n";
            failed++;
        }
    } catch (const std::exception& e) {
        std::cout << "FAIL: quoted CSV threw: " << e.what() << "\n";
        failed++;
    }

    // 3. Test JSON Unicode & Escape handling
    try {
        std::ofstream json_out(guard_json.path);
        json_out << "{\"greeting\": \"Hello\\nWorld\", \"symbol\": \"\\u0041\\u0042\"}";
        json_out.close();

        calc.try_process_function_command("jdata = read_json(\"" + guard_json.path.string() + "\")", &output);
        calc.try_process_function_command("jdata", &output);
        if (output.find("Hello") != std::string::npos && output.find("AB") != std::string::npos) {
            passed++;
        } else {
            std::cout << "FAIL: JSON escape read mismatch: " << output << "\n";
            failed++;
        }
    } catch (const std::exception& e) {
        std::cout << "FAIL: JSON escape threw: " << e.what() << "\n";
        failed++;
    }
}

void run_io_tests(int& passed, int& failed) {
    test_io_file_operations(passed, failed);
    test_io_csv_operations(passed, failed);
    test_io_json_operations(passed, failed);
    test_io_seek_tell(passed, failed);
    test_io_advanced_features(passed, failed);
}
