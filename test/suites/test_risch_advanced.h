#ifndef TEST_RISCH_ADVANCED_H
#define TEST_RISCH_ADVANCED_H

namespace test_suites {

void test_risch_advanced_independence();
void test_risch_nested();
void test_risch_special_part();
void test_risch_decision_procedure();
void test_parametric_rde();

inline void run_risch_advanced_tests() {
    test_risch_advanced_independence();
    test_risch_nested();
    test_risch_special_part();
    test_risch_decision_procedure();
    test_parametric_rde();
}

} // namespace test_suites

#endif // TEST_RISCH_ADVANCED_H
