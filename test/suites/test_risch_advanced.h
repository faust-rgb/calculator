#ifndef TEST_RISCH_ADVANCED_H
#define TEST_RISCH_ADVANCED_H

namespace test_suites {

void test_risch_advanced_independence();
void test_risch_nested();
void test_risch_special_part();
void test_risch_decision_procedure();
void test_parametric_rde();
void test_risch_high_degree_rational();
void test_risch_algebraic_curve_divisor();

inline void run_risch_advanced_tests() {
    test_risch_advanced_independence();
    test_risch_nested();
    test_risch_special_part();
    test_risch_decision_procedure();
    test_parametric_rde();
    test_risch_high_degree_rational();
    test_risch_algebraic_curve_divisor();
}

} // namespace test_suites

#endif // TEST_RISCH_ADVANCED_H
