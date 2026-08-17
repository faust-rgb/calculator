#include "math/functions/integer/absolute_value.h"

namespace mymath {

int abs(int x) { return x < 0 ? -x : x; }
long abs(long x) { return x < 0 ? -x : x; }
long long abs(long long x) { return x < 0 ? -x : x; }

}  // namespace mymath
