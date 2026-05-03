// ============================================================================
// Script Feature Tests (match-case, for-in, list/dict, logical operators)
// ============================================================================

#ifndef TEST_SCRIPT_FEATURES_H
#define TEST_SCRIPT_FEATURES_H

#include "core/calculator.h"
#include <iostream>

inline void test_script_match_case(int& passed, int& failed) {
    Calculator calc;
    std::string output;

    // Test basic match-case
    try {
        output = calc.execute_script(R"(
value = 2
result = 0
match value:
    case 1:
        result = 10
    case 2:
        result = 20
    case _:
        result = 30
print("match result: ", result)
)", false);

        if (output.find("20") != std::string::npos) {
            passed++;
        } else {
            std::cout << "FAIL: match-case basic: " << output << "\n";
            failed++;
        }
    } catch (const std::exception& e) {
        std::cout << "FAIL: match-case basic threw: " << e.what() << "\n";
        failed++;
    }

    // Test match with default case
    try {
        output = calc.execute_script(R"(
value = 100
result = 0
match value:
    case 1:
        result = 10
    case 2:
        result = 20
    case _:
        result = 99
)", false);

        if (output.find("99") != std::string::npos) {
            passed++;
        } else {
            std::cout << "FAIL: match-case default: " << output << "\n";
            failed++;
        }
    } catch (const std::exception& e) {
        std::cout << "FAIL: match-case default threw: " << e.what() << "\n";
        failed++;
    }
}

inline void test_script_for_in_list(int& passed, int& failed) {
    Calculator calc;
    std::string output;

    // Test for-in with list
    try {
        output = calc.execute_script(R"(
total = 0
for item in [1, 2, 3, 4, 5]:
    total = total + item
print("list sum: ", total)
)", false);

        if (output.find("15") != std::string::npos) {
            passed++;
        } else {
            std::cout << "FAIL: for-in list: " << output << "\n";
            failed++;
        }
    } catch (const std::exception& e) {
        std::cout << "FAIL: for-in list threw: " << e.what() << "\n";
        failed++;
    }
}

inline void test_script_for_in_matrix(int& passed, int& failed) {
    Calculator calc;
    std::string output;

    // Test for-in with matrix (row iteration)
    try {
        output = calc.execute_script(R"(
m = [1, 2; 3, 4; 5, 6]
row_count = 0
for row in m:
    row_count = row_count + 1
print("row count: ", row_count)
)", false);

        if (output.find("3") != std::string::npos) {
            passed++;
        } else {
            std::cout << "FAIL: for-in matrix: " << output << "\n";
            failed++;
        }
    } catch (const std::exception& e) {
        std::cout << "FAIL: for-in matrix threw: " << e.what() << "\n";
        failed++;
    }
}

inline void test_script_for_in_string(int& passed, int& failed) {
    Calculator calc;
    std::string output;

    // Test for-in with string
    try {
        output = calc.execute_script(R"(
count = 0
for ch in "hello":
    count = count + 1
print("char count: ", count)
)", false);

        if (output.find("5") != std::string::npos) {
            passed++;
        } else {
            std::cout << "FAIL: for-in string: " << output << "\n";
            failed++;
        }
    } catch (const std::exception& e) {
        std::cout << "FAIL: for-in string threw: " << e.what() << "\n";
        failed++;
    }
}

inline void test_script_list_operations(int& passed, int& failed) {
    Calculator calc;
    std::string output;

    // Test list indexing
    try {
        output = calc.execute_script(R"(
lst = [10, 20, 30, 40, 50]
first = lst[0]
last = lst[-1]
print("first: ", first, " last: ", last)
)", false);

        if (output.find("10") != std::string::npos && output.find("50") != std::string::npos) {
            passed++;
        } else {
            std::cout << "FAIL: list indexing: " << output << "\n";
            failed++;
        }
    } catch (const std::exception& e) {
        std::cout << "FAIL: list indexing threw: " << e.what() << "\n";
        failed++;
    }

    // Test list slice
    try {
        output = calc.execute_script(R"(
lst = [1, 2, 3, 4, 5]
sub = lst[1:4]
print("slice: ", sub)
)", false);

        if (output.find("2") != std::string::npos && output.find("3") != std::string::npos) {
            passed++;
        } else {
            std::cout << "FAIL: list slice: " << output << "\n";
            failed++;
        }
    } catch (const std::exception& e) {
        std::cout << "FAIL: list slice threw: " << e.what() << "\n";
        failed++;
    }

    // Test list assignment
    try {
        output = calc.execute_script(R"(
lst = [1, 2, 3]
lst[1] = 20
print("modified: ", lst)
)", false);

        if (output.find("20") != std::string::npos) {
            passed++;
        } else {
            std::cout << "FAIL: list assignment: " << output << "\n";
            failed++;
        }
    } catch (const std::exception& e) {
        std::cout << "FAIL: list assignment threw: " << e.what() << "\n";
        failed++;
    }
}

inline void test_script_dict_operations(int& passed, int& failed) {
    Calculator calc;
    std::string output;

    // Test dict creation and access
    try {
        output = calc.execute_script(R"(
d = {"a": 1, "b": 2, "c": 3}
val_a = d["a"]
val_b = d["b"]
print("a: ", val_a, " b: ", val_b)
)", false);

        if (output.find("1") != std::string::npos && output.find("2") != std::string::npos) {
            passed++;
        } else {
            std::cout << "FAIL: dict access: " << output << "\n";
            failed++;
        }
    } catch (const std::exception& e) {
        std::cout << "FAIL: dict access threw: " << e.what() << "\n";
        failed++;
    }

    // Test dict assignment
    try {
        output = calc.execute_script(R"(
d = {"x": 10}
d["y"] = 20
print("y: ", d["y"])
)", false);

        if (output.find("20") != std::string::npos) {
            passed++;
        } else {
            std::cout << "FAIL: dict assignment: " << output << "\n";
            failed++;
        }
    } catch (const std::exception& e) {
        std::cout << "FAIL: dict assignment threw: " << e.what() << "\n";
        failed++;
    }
}

inline void test_script_logical_operators(int& passed, int& failed) {
    Calculator calc;
    std::string output;

    // Test && operator
    try {
        output = calc.execute_script(R"(
x = 5
result = 0
if x > 0 && x < 10:
    result = 1
else:
    result = 0
print("and result: ", result)
)", false);

        if (output.find("1") != std::string::npos) {
            passed++;
        } else {
            std::cout << "FAIL: && operator: " << output << "\n";
            failed++;
        }
    } catch (const std::exception& e) {
        std::cout << "FAIL: && operator threw: " << e.what() << "\n";
        failed++;
    }

    // Test || operator
    try {
        output = calc.execute_script(R"(
x = 15
result = 0
if x < 0 || x > 10:
    result = 1
else:
    result = 0
print("or result: ", result)
)", false);

        if (output.find("1") != std::string::npos) {
            passed++;
        } else {
            std::cout << "FAIL: || operator: " << output << "\n";
            failed++;
        }
    } catch (const std::exception& e) {
        std::cout << "FAIL: || operator threw: " << e.what() << "\n";
        failed++;
    }

    // Test short-circuit evaluation
    try {
        output = calc.execute_script(R"(
x = 0
result = 1
if x != 0 && 10 / x > 1:
    result = 0
else:
    result = 1
print("short-circuit: ", result)
)", false);

        if (output.find("1") != std::string::npos) {
            passed++;
        } else {
            std::cout << "FAIL: short-circuit: " << output << "\n";
            failed++;
        }
    } catch (const std::exception& e) {
        std::cout << "FAIL: short-circuit threw: " << e.what() << "\n";
        failed++;
    }
}

inline void test_script_matrix_indexing(int& passed, int& failed) {
    Calculator calc;
    std::string output;

    // Test matrix indexing
    try {
        output = calc.execute_script(R"(
m = [1, 2, 3; 4, 5, 6]
val = m[0, 1]
print("m[0,1]: ", val)
)", false);

        if (output.find("2") != std::string::npos) {
            passed++;
        } else {
            std::cout << "FAIL: matrix indexing: " << output << "\n";
            failed++;
        }
    } catch (const std::exception& e) {
        std::cout << "FAIL: matrix indexing threw: " << e.what() << "\n";
        failed++;
    }

    // Test matrix assignment
    try {
        output = calc.execute_script(R"(
m = [1, 2; 3, 4]
m[0, 0] = 10
print("m[0,0]: ", m[0, 0])
)", false);

        if (output.find("10") != std::string::npos) {
            passed++;
        } else {
            std::cout << "FAIL: matrix assignment: " << output << "\n";
            failed++;
        }
    } catch (const std::exception& e) {
        std::cout << "FAIL: matrix assignment threw: " << e.what() << "\n";
        failed++;
    }

    // Test negative indexing
    try {
        output = calc.execute_script(R"(
m = [1, 2, 3; 4, 5, 6]
val = m[-1, -1]
print("m[-1,-1]: ", val)
)", false);

        if (output.find("6") != std::string::npos) {
            passed++;
        } else {
            std::cout << "FAIL: negative indexing: " << output << "\n";
            failed++;
        }
    } catch (const std::exception& e) {
        std::cout << "FAIL: negative indexing threw: " << e.what() << "\n";
        failed++;
    }
}

inline void run_script_feature_tests(int& passed, int& failed) {
    test_script_match_case(passed, failed);
    test_script_for_in_list(passed, failed);
    test_script_for_in_matrix(passed, failed);
    test_script_for_in_string(passed, failed);
    test_script_list_operations(passed, failed);
    test_script_dict_operations(passed, failed);
    test_script_logical_operators(passed, failed);
    test_script_matrix_indexing(passed, failed);
}

#endif // TEST_SCRIPT_FEATURES_H
