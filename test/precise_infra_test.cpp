#include "precise/precise_decimal.h"
#include <iostream>
#include <cassert>

void test_precise_operators() {
    PreciseDecimal a(10LL);
    PreciseDecimal b(20LL);
    
    // 加法
    assert((a + b).to_string() == "30");
    
    // 减法
    assert((a - b).to_string() == "-10");
    
    // 乘法
    assert((a * b).to_string() == "200");
    
    // 除法
    PreciseDecimal c(1LL);
    PreciseDecimal d(3LL);
    std::string div_res = (c / d).to_string();
    assert(div_res.substr(0, 10) == "0.33333333");
    
    // 混合运算
    PreciseDecimal e(0.5);
    assert((a * e).to_string() == "5");
    
    // 比较
    assert(a < b);
    assert(b > a);
    assert(a != b);
    assert(a == PreciseDecimal(10LL));
    
    // Sqrt
    PreciseDecimal two(2LL);
    PreciseDecimal s2 = precise::sqrt(two);
    std::string s2_str = s2.to_string();
    assert(s2_str.substr(0, 10) == "1.41421356");
    
    std::cout << "PreciseDecimal Infrastructure Tests Passed!" << std::endl;
}

int main() {
    try {
        test_precise_operators();
    } catch (const std::exception& e) {
        std::cerr << "Test failed: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
