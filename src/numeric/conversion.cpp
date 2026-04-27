#include "conversion.h"
#include "decimal.h"

#include <cstdio>
#include <cstdlib>
#include <string>

namespace numeric {

double to_double(const Number& value) {
    BigDecimal dec = value.to_decimal();
    std::string s = dec.to_string();
    return std::strtod(s.c_str(), nullptr);
}

bool can_convert_to_double(const Number& value) {
    BigDecimal dec = value.to_decimal();
    std::string s = dec.to_string();

    if (s == "inf" || s == "-inf" || s == "nan") {
        return false;
    }

    std::size_t e_pos = s.find('e');
    if (e_pos != std::string::npos) {
        int exp = std::stoi(s.substr(e_pos + 1));
        if (exp < -308 || exp > 308) {
            return false;
        }
    }

    return true;
}

Number from_double(double value) {
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%.17g", value);
    return Number(BigDecimal::from_string(buf));
}

long double to_long_double(const Number& value) {
    BigDecimal dec = value.to_decimal();
    std::string s = dec.to_string();
    return std::strtold(s.c_str(), nullptr);
}

}  // namespace numeric
