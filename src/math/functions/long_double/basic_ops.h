#pragma once

namespace mymath {

long double abs(long double x);
long double abs_long_double(long double x);
long double floor(long double x);
long double ceil(long double x);
long double round(long double x);
long double trunc(long double x);
long double modf(long double x, long double* integer_part);
long double clamp(long double value, long double low, long double high);
long double fmod(long double x, long double y);
long double remainder(long double x, long double y);
long double normalize_angle(long double x);

}
