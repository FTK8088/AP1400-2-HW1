#include "hw1.h"

namespace algebra {

    Matrix CreatMatrix(size_t n, size_t m, double v) {
        if (n <= 0 || m <= 0) {
            throw std::logic_error("Matrix dimensions must be positive.");
        }
        Matrix x;
        for (int i = 0; i < n; ++i) {
            std::vector<double> vec(m, v);
            x.push_back(vec);
        }
        return x;
    }

    Matrix zeros(size_t n, size_t m) {
        return CreatMatrix(n, m, 0);
    }

    Matrix ones(size_t n, size_t m) {
        return CreatMatrix(n, m, 1);
    }

    Matrix random(size_t n, size_t m, double min, double max) {
        if (min > max) {
            throw std::logic_error("Min can't be greater than max\n");
        }
        std::default_random_engine rng;
        std::uniform_real_distribution<double> dist(min, max);
        Matrix x = zeros(n, m);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                x[i][j] = dist(rng);
            }
        }
        return x;
    }

    void show(const Matrix& matrix) {
        int row = matrix.size();
        int column = matrix[0].size();
        for (int i = 0; i < row; i++) {
            for (int j = 0; j < column; j++) {
                // std::cout << matrix[i][j] << std::endl;
                std::cout << std::setw(7) << std::fixed << std::setprecision(3) << matrix[i][j];
            }
            std::cout << std::endl;
        }
    }

    Matrix multiply(const Matrix& matrix, double c) {
        int row = matrix.size();
        int column = matrix[0].size();
        if (c == 0) {
            return CreatMatrix(row, column, 0);
        }

        Matrix x = CreatMatrix(row, column, 0);
        for (int i = 0; i < row; i++) {
            for (int j = 0; j < column; j++) {
                x[i][j] = matrix[i][j] * c;
            }
        }
        return x;
    }

    Matrix multiply(const Matrix& matrix1, const Matrix& matrix2) {
        if (matrix1.empty() && matrix2.empty()) {
            return Matrix{};
        }

        if (matrix1.empty() || matrix2.empty()) {
            throw logic_error("calculation error\n");
        }

        if (matrix1[0].size() != matrix2.size()) {
            throw logic_error("dimenssion errors\n");
        }

        int row = matrix1.size();
        int column = matrix2[0].size();
        Matrix x = zeros(row, column);

        for (int i = 0; i < row; i++) {
            for (int j = 0; j < column; j++) {
                for (int k = 0; k < matrix2.size(); k++) {
                    x[i][j] += matrix1[i][k] * matrix2[k][j];
                }
            }
        }
        return x;
    }
    
    Matrix sum(const Matrix& matrix, double c) {
        if (matrix.empty()) {
            return Matrix{};
        }

        int row = matrix.size();
        int column = matrix[0].size();
        Matrix x = zeros(row, column);
        
        for (int i = 0; i < row; i++ ){
            for (int j = 0; j < column; j++) {
                x[i][j] = matrix[i][j] + c;
            }
        }
        return x;
    }

    Matrix sum(const Matrix& matrix1, const Matrix& matrix2) {
        if (matrix1.empty() && matrix2.empty()) {
            return Matrix{};
        }

        if (matrix1.empty() || matrix2.empty()) {
            throw logic_error("calculation error\n");
        }

        if (matrix1.size() != matrix2.size() || matrix1[0].size() != matrix2[0].size()) {
            throw logic_error("dimenssion error\n");
        }

        int row = matrix1.size();
        int column = matrix1[0].size();
        Matrix x = zeros(row, column);
        for (int i = 0; i < row; i++) {
            for (int j = 0; j < column; j++) {
                x[i][j] = matrix1[i][j] + matrix2[i][j];
            }
        }
        return x;
    }

    Matrix transpose(const Matrix& matrix) {
        if (matrix.empty()) {
            return Matrix{};
        }
        int row = matrix.size();
        int column = matrix[0].size();
        Matrix x = zeros(column, row);
        for (int i = 0; i < row; i++) {
            for (int j = 0; j < column; j++) {
                x[j][i] = matrix[i][j];
            }
        }
        return x;
    }

    Matrix minor(const Matrix& matrix, size_t n, size_t m) {
        if (n < 0 || n >= matrix.size() || m < 0 || m >= matrix[0].size()) {
            throw logic_error("out of bound\n");
        }

        int row = matrix.size();
        int column = matrix[0].size();
        Matrix x{};
        for (int i = 0; i < row; i++) {
            if (i == n) {
                continue;
            }
            std::vector<double> v;
            for (int j = 0; j < column; j++) {
                if (j == m) {
                    continue;
                }
                v.push_back(matrix[i][j]);
            }
            x.push_back(v);
        }
        return x;
    }

    double determinant(const Matrix& matrix) {
        if (matrix.empty()) {
            return 1;
        }
        if (matrix.size() != matrix[0].size()) {
            throw logic_error("not a square matrix");
        }
        if (matrix.size() == 1) {
            return matrix[0][0];
        } else if (matrix.size() == 2) {
            return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
        } else {
            double res = 0;
            for (int i = 0; i < matrix.size(); i++) {
                res += matrix[i][0] * pow(-1, i) * determinant(minor(matrix, i, 0));
            }
            return res;
        }
        return 0;
    }

    Matrix inverse(const Matrix& matrix) {
        if (matrix.empty()) {
            return Matrix{};
        }
        if (matrix.size() != matrix[0].size() || determinant(matrix) == 0) {
            throw logic_error("The inverse of matrix doesn't exist\n");
        }
        Matrix x = zeros(matrix.size(), matrix.size());
        for (int i = 0; i < matrix.size(); i++) {
            for (int j = 0; j < matrix.size(); j++) {
                x[i][j] = determinant(minor(matrix, i, j)) * pow(-1, i + j);
            }
        }
    
        x = transpose(x);
    
        return multiply(x, 1.0 / determinant(matrix));
    }

    Matrix concatenate(const Matrix& matrix1, const Matrix& matrix2, int axis) {
        if (axis == 0) {
            if (matrix1[0].size() != matrix2[0].size()) {
                throw logic_error("the number of columns has to be consistent\n");
            }
            Matrix x;
            for (int i = 0; i < matrix1.size(); i++) {
                x.push_back(matrix1[i]);
            }
            for (int j = 0; j < matrix2.size(); j++) {
                x.push_back(matrix2[j]);
            }
            return x;
        } else if (axis == 1) {
            if (matrix1.size() != matrix2.size()) {
                throw logic_error("the number of rows has to be consistent\n");
            }
            Matrix x = zeros(matrix1.size(), matrix1[0].size() + matrix2[0].size());
            for (int i = 0; i < matrix1.size(); i++) {
                for (int j = 0; j < x[0].size(); j++) {
                    if (j < matrix1[0].size()) {
                        x[i][j] = matrix1[i][j];
                    } else {
                        x[i][j] = matrix2[i][j - matrix1[0].size()];
                    }
                }
            }
            return x;
        } else {
            throw logic_error("axis must equal to 0 or 1");
        }
    }

    Matrix ero_swap(const Matrix& matrix, size_t r1, size_t r2) {
        
    }

}