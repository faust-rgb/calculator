#include "unit_conversions.h"
#include "mymath.h"

long double degrees_to_radians(long double value) {
    return value * mymath::kPi / 180.0L;
}

long double radians_to_degrees(long double value) {
    return value * 180.0L / mymath::kPi;
}

long double celsius_to_fahrenheit(long double value) {
    return value * 9.0 / 5.0 + 32.0;
}

long double fahrenheit_to_celsius(long double value) {
    return (value - 32.0) * 5.0 / 9.0;
}
