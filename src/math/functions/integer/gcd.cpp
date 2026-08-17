#include "math/functions/integer/gcd.h"

namespace mymath {

long long gcd(long long a, long long b) {
    while (b != 0) {
        const long long remainder = a % b;
        a = b;
        b = remainder;
    }
    return a < 0 ? -a : a;
}

}  // namespace mymath
