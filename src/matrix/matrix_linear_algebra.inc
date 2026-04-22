Matrix inverse(const Matrix& matrix) {
    if (!matrix.is_square()) {
        throw std::runtime_error("inverse requires a square matrix");
    }

    // 通过对 [A | I] 做 Gauss-Jordan 消元，把左半边化成 I，
    // 右半边就会变成 A 的逆矩阵。
    const std::size_t n = matrix.rows;
    Matrix augmented(n, n * 2, 0.0);
    for (std::size_t row = 0; row < n; ++row) {
        for (std::size_t col = 0; col < n; ++col) {
            augmented.at(row, col) = matrix.at(row, col);
            augmented.at(row, n + col) = row == col ? 1.0 : 0.0;
        }
    }

    for (std::size_t col = 0; col < n; ++col) {
        std::size_t best_row = col;
        double best_value = mymath::abs(augmented.at(best_row, col));
        for (std::size_t row = col + 1; row < n; ++row) {
            const double current = mymath::abs(augmented.at(row, col));
            if (current > best_value) {
                best_value = current;
                best_row = row;
            }
        }

        if (mymath::is_near_zero(best_value, kMatrixEps)) {
            throw std::runtime_error("matrix is singular and cannot be inverted");
        }

        swap_rows(&augmented, col, best_row);
        const double pivot = augmented.at(col, col);
        for (std::size_t current_col = 0; current_col < augmented.cols; ++current_col) {
            augmented.at(col, current_col) /= pivot;
        }

        for (std::size_t row = 0; row < n; ++row) {
            if (row == col) {
                continue;
            }
            const double factor = augmented.at(row, col);
            if (mymath::is_near_zero(factor, kMatrixEps)) {
                continue;
            }
            for (std::size_t current_col = 0; current_col < augmented.cols; ++current_col) {
                augmented.at(row, current_col) -= factor * augmented.at(col, current_col);
                if (mymath::is_near_zero(augmented.at(row, current_col), kMatrixEps)) {
                    augmented.at(row, current_col) = 0.0;
                }
            }
        }
    }

    Matrix result(n, n, 0.0);
    for (std::size_t row = 0; row < n; ++row) {
        for (std::size_t col = 0; col < n; ++col) {
            result.at(row, col) = augmented.at(row, n + col);
        }
    }
    return result;
}

Matrix pseudo_inverse(const Matrix& matrix) {
    const ReducedSvd svd = compute_reduced_svd(matrix);
    Matrix sigma_pinv(svd.s.cols, svd.s.rows, 0.0);
    const std::size_t diagonal =
        svd.s.rows < svd.s.cols ? svd.s.rows : svd.s.cols;
    for (std::size_t i = 0; i < diagonal; ++i) {
        const double sigma = svd.s.at(i, i);
        if (!mymath::is_near_zero(sigma, kMatrixEps)) {
            sigma_pinv.at(i, i) = 1.0 / sigma;
        }
    }
    return multiply(multiply(transpose(svd.vt), sigma_pinv), transpose(svd.u));
}

Matrix nullspace(const Matrix& matrix) {
    return nullspace_basis(matrix);
}

Matrix least_squares(const Matrix& coefficients, const Matrix& rhs) {
    if (rhs.cols != 1 && rhs.rows != 1) {
        throw std::runtime_error("least_squares currently requires the right-hand side to be a vector");
    }

    const std::size_t rhs_size = rhs.rows == 1 ? rhs.cols : rhs.rows;
    if (rhs_size != coefficients.rows) {
        throw std::runtime_error("least_squares requires rhs to match the number of rows in A");
    }

    Matrix rhs_column(coefficients.rows, 1, 0.0);
    for (std::size_t row = 0; row < coefficients.rows; ++row) {
        rhs_column.at(row, 0) = rhs.rows == 1 ? rhs.at(0, row) : rhs.at(row, 0);
    }
    return multiply(pseudo_inverse(coefficients), rhs_column);
}

Matrix qr_q(const Matrix& matrix) {
    return qr_decompose(matrix).first;
}

Matrix qr_r(const Matrix& matrix) {
    return qr_decompose(matrix).second;
}

Matrix lu_l(const Matrix& matrix) {
    return lu_decompose(matrix).first;
}

Matrix lu_u(const Matrix& matrix) {
    return lu_decompose(matrix).second;
}

Matrix svd_u(const Matrix& matrix) {
    return compute_reduced_svd(matrix).u;
}

Matrix svd_s(const Matrix& matrix) {
    return compute_reduced_svd(matrix).s;
}

Matrix svd_vt(const Matrix& matrix) {
    return compute_reduced_svd(matrix).vt;
}

Matrix solve(const Matrix& coefficients, const Matrix& rhs) {
    if (!coefficients.is_square()) {
        throw std::runtime_error("solve requires a square coefficient matrix");
    }
    if (rhs.cols != 1 && rhs.rows != 1) {
        throw std::runtime_error("solve currently requires the right-hand side to be a vector");
    }

    const std::size_t n = coefficients.rows;
    const std::size_t rhs_size = rhs.rows == 1 ? rhs.cols : rhs.rows;
    if (rhs_size != n) {
        throw std::runtime_error("solve requires rhs to match the coefficient matrix dimension");
    }

    // 对增广矩阵 [A | b] 做 Gauss-Jordan 消元，右端最终就是解向量。
    Matrix augmented(n, n + 1, 0.0);
    for (std::size_t row = 0; row < n; ++row) {
        for (std::size_t col = 0; col < n; ++col) {
            augmented.at(row, col) = coefficients.at(row, col);
        }
        augmented.at(row, n) = rhs.rows == 1 ? rhs.at(0, row) : rhs.at(row, 0);
    }

    for (std::size_t col = 0; col < n; ++col) {
        std::size_t best_row = col;
        double best_value = mymath::abs(augmented.at(best_row, col));
        for (std::size_t row = col + 1; row < n; ++row) {
            const double current = mymath::abs(augmented.at(row, col));
            if (current > best_value) {
                best_value = current;
                best_row = row;
            }
        }

        if (mymath::is_near_zero(best_value, kMatrixEps)) {
            throw std::runtime_error("linear system has no unique solution");
        }

        swap_rows(&augmented, col, best_row);
        const double pivot = augmented.at(col, col);
        for (std::size_t current_col = 0; current_col < augmented.cols; ++current_col) {
            augmented.at(col, current_col) /= pivot;
        }

        for (std::size_t row = 0; row < n; ++row) {
            if (row == col) {
                continue;
            }
            const double factor = augmented.at(row, col);
            if (mymath::is_near_zero(factor, kMatrixEps)) {
                continue;
            }
            for (std::size_t current_col = 0; current_col < augmented.cols; ++current_col) {
                augmented.at(row, current_col) -= factor * augmented.at(col, current_col);
                if (mymath::is_near_zero(augmented.at(row, current_col), kMatrixEps)) {
                    augmented.at(row, current_col) = 0.0;
                }
            }
        }
    }

    Matrix result(n, 1, 0.0);
    for (std::size_t row = 0; row < n; ++row) {
        result.at(row, 0) = augmented.at(row, n);
    }
    return result;
}

Matrix power(Matrix base, long long exponent) {
    if (!base.is_square()) {
        throw std::runtime_error("matrix powers require a square matrix");
    }

    if (exponent < 0) {
        base = inverse(base);
        exponent = -exponent;
    }

    Matrix result = Matrix::identity(base.rows);
    while (exponent > 0) {
        if ((exponent & 1LL) != 0) {
            result = multiply(result, base);
        }
        base = multiply(base, base);
        exponent >>= 1LL;
    }
    return result;
}

double condition_number(const Matrix& matrix) {
    const std::size_t effective_rank =
        static_cast<std::size_t>(rank(matrix));
    const std::size_t full_rank = matrix.rows < matrix.cols ? matrix.rows : matrix.cols;
    if (effective_rank < full_rank) {
        return std::numeric_limits<double>::infinity();
    }

    const Matrix singular_values = svd_s(matrix);
    double largest = 0.0;
    double smallest = std::numeric_limits<double>::infinity();
    const std::size_t diagonal =
        singular_values.rows < singular_values.cols ? singular_values.rows : singular_values.cols;
    for (std::size_t i = 0; i < diagonal; ++i) {
        const double sigma = mymath::abs(singular_values.at(i, i));
        if (sigma > largest) {
            largest = sigma;
        }
        if (!mymath::is_near_zero(sigma, kMatrixEps) && sigma < smallest) {
            smallest = sigma;
        }
    }
    if (mymath::is_near_zero(largest, kMatrixEps) ||
        smallest == std::numeric_limits<double>::infinity()) {
        return std::numeric_limits<double>::infinity();
    }
    return largest / smallest;
}

Matrix cholesky(const Matrix& matrix) {
    if (!matrix.is_square()) {
        throw std::runtime_error("cholesky requires a square matrix");
    }
    Matrix result(matrix.rows, matrix.cols, 0.0);
    for (std::size_t i = 0; i < matrix.rows; ++i) {
        for (std::size_t j = 0; j <= i; ++j) {
            double sum = matrix.at(i, j);
            for (std::size_t k = 0; k < j; ++k) {
                sum -= result.at(i, k) * result.at(j, k);
            }
            if (i == j) {
                if (sum <= 0.0) {
                    throw std::runtime_error("cholesky requires a positive-definite matrix");
                }
                result.at(i, j) = mymath::sqrt(sum);
            } else {
                result.at(i, j) = sum / result.at(j, j);
            }
        }
    }
    return result;
}

Matrix hessenberg(const Matrix& matrix) {
    if (!matrix.is_square()) {
        throw std::runtime_error("hessenberg requires a square matrix");
    }

    Matrix h = matrix;
    const std::size_t n = matrix.rows;
    if (n <= 2) {
        return h;
    }

    for (std::size_t k = 0; k + 2 < n; ++k) {
        std::vector<double> x(n - k - 1, 0.0);
        for (std::size_t i = k + 1; i < n; ++i) {
            x[i - k - 1] = h.at(i, k);
        }

        double norm_x = 0.0;
        for (double value : x) {
            norm_x += value * value;
        }
        norm_x = mymath::sqrt(norm_x);
        if (mymath::is_near_zero(norm_x, kMatrixEps)) {
            continue;
        }

        x[0] += (x[0] >= 0.0 ? norm_x : -norm_x);
        double norm_v = 0.0;
        for (double value : x) {
            norm_v += value * value;
        }
        norm_v = mymath::sqrt(norm_v);
        if (mymath::is_near_zero(norm_v, kMatrixEps)) {
            continue;
        }
        for (double& value : x) {
            value /= norm_v;
        }

        for (std::size_t j = k; j < n; ++j) {
            double projection = 0.0;
            for (std::size_t i = 0; i < x.size(); ++i) {
                projection += x[i] * h.at(k + 1 + i, j);
            }
            projection *= 2.0;
            for (std::size_t i = 0; i < x.size(); ++i) {
                h.at(k + 1 + i, j) -= projection * x[i];
            }
        }

        for (std::size_t i = 0; i < n; ++i) {
            double projection = 0.0;
            for (std::size_t j = 0; j < x.size(); ++j) {
                projection += h.at(i, k + 1 + j) * x[j];
            }
            projection *= 2.0;
            for (std::size_t j = 0; j < x.size(); ++j) {
                h.at(i, k + 1 + j) -= projection * x[j];
            }
        }

        for (std::size_t i = k + 2; i < n; ++i) {
            if (mymath::is_near_zero(h.at(i, k), kMatrixEps)) {
                h.at(i, k) = 0.0;
            }
        }
    }
    return h;
}

Matrix schur(const Matrix& matrix) {
    if (!matrix.is_square()) {
        throw std::runtime_error("schur requires a square matrix");
    }
    Matrix current = hessenberg(matrix);
    for (int iteration = 0; iteration < 256; ++iteration) {
        const auto qr = qr_decompose(current);
        current = multiply(qr.second, qr.first);
        if (off_diagonal_magnitude(current) < 1e-8) {
            break;
        }
    }
    return current;
}

double get(const Matrix& matrix, std::size_t row, std::size_t col) {
    return matrix.at(row, col);
}

double get(const Matrix& matrix, std::size_t index) {
    if (!matrix.is_vector()) {
        throw std::runtime_error("single-index get only works on vectors");
    }
    if (matrix.rows == 1) {
        return matrix.at(0, index);
    }
    return matrix.at(index, 0);
}

Matrix set(Matrix matrix, std::size_t row, std::size_t col, double value) {
    if (row >= matrix.rows || col >= matrix.cols) {
        const std::size_t new_rows = row < matrix.rows ? matrix.rows : row + 1;
        const std::size_t new_cols = col < matrix.cols ? matrix.cols : col + 1;
        matrix.resize(new_rows, new_cols, 0.0);
    }
    matrix.at(row, col) = value;
    return matrix;
}

Matrix set(Matrix matrix, std::size_t index, double value) {
    if (!matrix.is_vector()) {
        throw std::runtime_error("single-index set only works on vectors");
    }
    if (matrix.rows == 1) {
        if (index >= matrix.cols) {
            matrix.resize(1, index + 1, 0.0);
        }
        matrix.at(0, index) = value;
    } else {
        if (index >= matrix.rows) {
            matrix.resize(index + 1, 1, 0.0);
        }
        matrix.at(index, 0) = value;
    }
    return matrix;
}

double norm(const Matrix& matrix) {
    double sum = 0.0;
    for (double value : matrix.data) {
        sum += value * value;
    }
    return mymath::sqrt(sum);
}

double trace(const Matrix& matrix) {
    if (!matrix.is_square()) {
        throw std::runtime_error("trace requires a square matrix");
    }
    double sum = 0.0;
    for (std::size_t i = 0; i < matrix.rows; ++i) {
        sum += matrix.at(i, i);
    }
    return sum;
}

double determinant(const Matrix& matrix) {
    if (!matrix.is_square()) {
        throw std::runtime_error("determinant requires a square matrix");
    }

    // 通过带部分主元选取的上三角消元求行列式。
    // 行交换会翻转符号，对角线乘积给出最终结果。
    Matrix reduced = matrix;
    double result = 1.0;
    int swap_count = 0;

    for (std::size_t col = 0; col < reduced.cols; ++col) {
        std::size_t best_row = col;
        double best_value = mymath::abs(reduced.at(best_row, col));
        for (std::size_t row = col + 1; row < reduced.rows; ++row) {
            const double current = mymath::abs(reduced.at(row, col));
            if (current > best_value) {
                best_value = current;
                best_row = row;
            }
        }

        if (mymath::is_near_zero(best_value, kMatrixEps)) {
            return 0.0;
        }

        if (best_row != col) {
            swap_rows(&reduced, best_row, col);
            ++swap_count;
        }

        const double pivot = reduced.at(col, col);
        result *= pivot;
        for (std::size_t row = col + 1; row < reduced.rows; ++row) {
            const double factor = reduced.at(row, col) / pivot;
            for (std::size_t current_col = col; current_col < reduced.cols; ++current_col) {
                reduced.at(row, current_col) -= factor * reduced.at(col, current_col);
            }
        }
    }

    return (swap_count % 2 == 0 ? 1.0 : -1.0) * result;
}

double rank(const Matrix& matrix) {
    Matrix reduced = matrix;
    return static_cast<double>(rref_in_place(&reduced).size());
}

Matrix rref(Matrix matrix) {
    rref_in_place(&matrix);
    return matrix;
}

Matrix eigenvalues(const Matrix& matrix) {
    if (!matrix.is_square()) {
        throw std::runtime_error("eigvals requires a square matrix");
    }

    if (matrix.rows == 0) {
        return Matrix::vector(std::vector<double>());
    }

    if (matrix.rows == 1) {
        return Matrix::vector({matrix.at(0, 0)});
    }

    if (matrix.rows == 2) {
        // 2x2 情况直接用特征多项式闭式解，结果更稳定也更快。
        const double a = 1.0;
        const double b = -(matrix.at(0, 0) + matrix.at(1, 1));
        const double c = determinant(matrix);
        const double discriminant = b * b - 4.0 * a * c;
        if (discriminant < -kMatrixEps) {
            throw std::runtime_error("eigvals only supports matrices with real eigenvalues");
        }
        const double root = mymath::sqrt(discriminant < 0.0 ? 0.0 : discriminant);
        return Matrix::vector({(-b + root) / 2.0, (-b - root) / 2.0});
    }

    Matrix current = matrix;
    // 更高阶情况走未移位的 QR 迭代。
    // 这不是最强健的工业实现，但足够覆盖当前项目里的中小型实矩阵场景。
    for (int iteration = 0; iteration < 128; ++iteration) {
        const auto qr = qr_decompose(current);
        current = multiply(qr.second, qr.first);
        if (off_diagonal_magnitude(current) < 1e-8) {
            break;
        }
    }

    std::vector<double> values;
    values.reserve(current.rows);
    for (std::size_t i = 0; i < current.rows; ++i) {
        values.push_back(current.at(i, i));
    }
    return Matrix::vector(values);
}

Matrix eigenvectors(const Matrix& matrix) {
    if (!matrix.is_square()) {
        throw std::runtime_error("eigvecs requires a square matrix");
    }

    // 先求特征值，再逐个解 (A - lambda I)v = 0。
    // 返回结果按“列向量矩阵”组织：每一列是一个特征向量。
    const Matrix values = eigenvalues(matrix);
    Matrix vectors(matrix.rows, matrix.cols, 0.0);
    for (std::size_t col = 0; col < values.cols; ++col) {
        Matrix shifted = matrix;
        for (std::size_t i = 0; i < shifted.rows; ++i) {
            shifted.at(i, i) -= values.at(0, col);
        }

        const std::vector<double> basis = nullspace_vector(shifted);
        for (std::size_t row = 0; row < basis.size(); ++row) {
            vectors.at(row, col) = basis[row];
        }
    }
    return vectors;
}
