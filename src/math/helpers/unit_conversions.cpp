#include "unit_conversions.h"
#include "mymath.h"

double degrees_to_radians(double value) {
    return value * mymath::kPi / 180.0;
}

double radians_to_degrees(double value) {
    return value * 180.0 / mymath::kPi;
}

double celsius_to_fahrenheit(double value) {
    return value * 9.0 / 5.0 + 32.0;
}

double fahrenheit_to_celsius(double value) {
    return (value - 32.0) * 5.0 / 9.0;
}
