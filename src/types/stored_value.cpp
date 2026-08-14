// ============================================================================
// StoredValue 代理实现 (src/types/stored_value.cpp)
// ============================================================================

#include "types/stored_value.h"
#include "matrix.h"

StoredValue::StoredValue(matrix::Matrix val)
    : data(std::make_shared<matrix::Matrix>(std::move(val))) {
}

StoredValue::IsMatrixProxy::operator bool() const {
    return owner()->is_matrix_type();
}

StoredValue::IsMatrixProxy& StoredValue::IsMatrixProxy::operator=(bool set_mat) {
    if (set_mat && !owner()->is_matrix_type()) {
        owner()->data = std::make_shared<matrix::Matrix>();
    } else if (!set_mat && owner()->is_matrix_type()) {
        owner()->data = std::monostate{};
    }
    return *this;
}
