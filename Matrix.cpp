#include <algorithm>
#include <iostream>
#include <vector>
#include <stack>

template<typename T>
class MatrixIt;

template<typename T>
class Matrix {
 private:
    std::vector<std::vector<T>> data;

 public:
    Matrix(const std::vector<std::vector<T>> &d) : data(d) {
        SetMatrix(d);
    }

    void SetMatrix(const std::vector<std::vector<T>> &d) {
        data.resize(d.size());
        for (size_t i = 0; i != d.size(); i++) {
            data[i].resize(d[0].size());
            for (size_t j = 0; j != d[i].size(); j++) {
                data[i][j] = d[i][j];
            }
        }
    }

    std::pair<size_t, size_t> size() const {
        std::pair<size_t, size_t> s;
        s.first = data.size();
        s.second = data[0].size();
        return s;
    }

    T &operator()(size_t i, size_t j) {
        return data[i][j];
    }

    const T &operator()(size_t i, size_t j) const {
        return data[i][j];
    }

    Matrix &operator+=(const Matrix &other) {
        for (size_t i = 0; i != other.size().first; i++) {
            for (size_t j = 0; j != other.size().second; j++) {
                data[i][j] += other.data[i][j];
            }
        }
        return *this;
    }

    Matrix operator+(const Matrix &other) const {
        auto C = *this;
        C += other;
        return C;
    }

    Matrix &operator*=(const int &a) {
        for (size_t i = 0; i != data.size(); i++) {
            for (size_t j = 0; j != data[0].size(); j++) {
                data[i][j] *= a;
            }
        }
        return *this;
    }

    Matrix &operator*=(const long int &a) {
        for (size_t i = 0; i != data.size(); i++) {
            for (size_t j = 0; j != data[0].size(); j++) {
                data[i][j] *= a;
            }
        }
        return *this;
    }

    Matrix &operator*=(const double &a) {
        for (size_t i = 0; i != data.size(); i++) {
            for (size_t j = 0; j != data[0].size(); j++) {
                data[i][j] *= a;
            }
        }
        return *this;
    }

    Matrix &operator*=(const float &a) {
        for (size_t i = 0; i != data.size(); i++) {
            for (size_t j = 0; j != data[0].size(); j++) {
                data[i][j] *= a;
            }
        }
        return *this;
    }

    Matrix operator*(const long int &a) const {
        auto C = *this;
        C *= a;
        return C;
    }

    Matrix operator*(const int &a) const {
        auto C = *this;
        C *= a;
        return C;
    }

    Matrix operator*(const double &a) const {
        auto C = *this;
        C *= a;
        return C;
    }

    Matrix operator*(const float &a) const {
        auto C = *this;
        C *= a;
        return C;
    }

    const Matrix transposed() const {
        std::vector<std::vector<T>> vec(data[0].size(), std::vector<T>(data.size(), 0));
        for (size_t i = 0; i != data.size(); i++) {
            for (size_t j = 0; j != data[0].size(); j++) {
                vec[j][i] = data[i][j];
            }
        }
        Matrix<T> C(vec);
        return C;
    }

    Matrix &transpose() {
        auto C = (*this).transposed();
        *this = C;
        return *this;
    }

    const Matrix transpose() const {
        auto C = (*this).transposed();
        return C;
    }

    Matrix &operator*=(const Matrix &other) {
        std::vector<std::vector<T>> vec(data.size(), std::vector<T>(other.size().second, 0));
        Matrix<T> C(vec);
        for (size_t i = 0; i != C.size().first; i++) {
            for (size_t j = 0; j != C.size().second; j++) {
                for (size_t k = 0; k != other.size().first; k++) {
                    C(i, j) += data[i][k] * other(k, j);
                }
            }
        }
        *this = C;
        return *this;
    }

    const Matrix operator*(const Matrix &other) const {
        auto C = *this;
        C *= other;
        return C;
    }

    MatrixIt<T> begin() {
        return MatrixIt<T>(this, 0, 0);
    }

    const MatrixIt<T> begin() const {
        return MatrixIt<T>(this, 0, 0);
    }

    MatrixIt<T> end() {
        return MatrixIt<T>(this);
    }

    const MatrixIt<T> end() const {
        return MatrixIt<T>(this);
    }

    void swap(size_t i, size_t j) {
        T x;
        for (size_t g = 0; g != data[0].size(); g++) {
            x = data[i][g];
            data[i][g] = data[j][g];
            data[j][g] = x;
        }
    }

    template<typename U>
    std::vector<U> solve(const std::vector<U> &b) {
        std::vector<std::vector<U>> vec(data.size(), std::vector<U>(data.size() + 1, 0));
        size_t N = data.size();
        for (size_t i = 0; i != N; i++) {
            for (size_t j = 0; j != N; j++)
                vec[i][j] = static_cast<U> (data[i][j]);
            vec[i][N] = b[i];
        }
        Matrix<U> newMatrix(vec);
        for (size_t column = 0; column != N; column++) {
            bool nonzero_found = true;
            size_t k = column;
            size_t row = k + 1;
            while ((row < N) && (nonzero_found)) {
                if ((newMatrix(row, column) != 0) && (nonzero_found)) {
                    newMatrix.swap(row, k);
                    nonzero_found = false;
                }
                row++;
            }
            for (size_t g = k + 1; g != N; g++) {
                U c = newMatrix(g, column) / newMatrix(column, column);
                for (size_t f = column; f != N + 1; f++) {
                    newMatrix(g, f) -= newMatrix(column, f) * c;
                }
            }
        }
        for (size_t g = 0; g != N; g++) {
            U q = newMatrix(g, g);
            for (size_t f = g; f != N + 1; f++) {
                newMatrix(g, f) /= q;
            }
        }
        for (int i = 0; i != N; i++) {
            for (int j = 0; j != N; j++) {
                if (i == j)
                    continue;
                U q = newMatrix(j, i) / newMatrix(i, i);
                for (int k = 0; k < N + 1; k++) {
                    newMatrix(j, k) -= newMatrix(i, k) * q;
                }
            }
        }
        std::vector<U> result(N, 0);
        for (size_t i = 0; i != N; i++) {
            result[i] = newMatrix(i, N);
        }
        return result;
    }
};

template<typename T>
class MatrixIt {
 private:
    const Matrix<T> *matrix;
    size_t row, column;

 public:
    MatrixIt(Matrix<T> *imatrix) :
            matrix(imatrix), row(imatrix->size().first), column(0) {
    }

    MatrixIt(Matrix<T> *imatrix, size_t irow, size_t icolumn = 0) :
            matrix(imatrix), row(irow), column(icolumn) {
    }

    MatrixIt(const Matrix<T> *imatrix, size_t irow, size_t icolumn = 0) :
            matrix(imatrix), row(irow), column(icolumn) {
    }

    MatrixIt(const Matrix<T> *imatrix) :
            matrix(imatrix), row(imatrix->size().first), column(0) {
    }

    const T &operator*() const {
        return (*matrix)(row, column);
    }

    MatrixIt operator++() {
        ++column;
        if (column == matrix->size().second) {
            column = 0;
            row++;
        }
        return *this;
    }

    MatrixIt operator++(int) {
        auto temp = *this;
        ++*this;
        return temp;
    }

    bool operator==(MatrixIt other) const {
        return ((row == other.row) && (column == other.column));
    }

    bool operator!=(MatrixIt other) const {
        return !(*this == other);
    }
};


template<typename T>
std::ostream &operator<<(std::ostream &out, const Matrix<T> &m) {
    for (size_t i = 0; i != m.size().first; i++) {
        for (size_t j = 0; j != m.size().second; j++) {
            out << m(i, j);
            if ((j != m.size().second - 1) || (m.size().second == 1)) {
                out << '\t';
            }
        }
        if (i != m.size().first - 1) {
            out << '\n';
        }
    }
    return out;
}
