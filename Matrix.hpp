//
// Created by Stepan Didurenko on 20.04.2024.
//

#ifndef INC_23_HOME_2_MATRIX_HPP
#define INC_23_HOME_2_MATRIX_HPP


#include <cstdlib> // for size_t
#include <iostream> // for std::ostream, std::istream
#include <array> // for shape
#include <stdexcept> // for exceptions
#include <functional> // for std::function
#include <vector> // for std::vector constructor

template<typename T>
class Matrix {
private:
    size_t rows = 0;
    size_t cols = 0;
    T **data = nullptr;
public:
    // Constructors
    Matrix() noexcept {
        rows = 0;
        cols = 0;
        data = nullptr;
    }

    Matrix(size_t rows, size_t cols) {
        if (rows == 0 || cols == 0)
            throw std::invalid_argument("Matrix dimensions must be positive");
        this->rows = rows;
        this->cols = cols;
        data = new T *[rows];
        for (size_t i = 0; i < rows; i++)
            data[i] = new T[cols];
    }

    explicit Matrix(size_t size) : Matrix(size, size) {};

    Matrix(const Matrix<T> &other) {
        rows = other.rows;
        cols = other.cols;
        data = new T *[rows];
        for (size_t i = 0; i < rows; i++) {
            data[i] = new T[cols];
            for (size_t j = 0; j < cols; j++)
                data[i][j] = other.data[i][j];
        }
    }

    explicit Matrix(const std::vector<std::vector<T>> &v) {
        if (rows == 0 || cols == 0)
            throw std::invalid_argument("Matrix dimensions must be positive");
        rows = v.size();
        cols = v[0].size();
        data = new T *[rows];
        for (size_t i = 0; i < rows; i++) {
            data[i] = new T[cols];
            for (size_t j = 0; j < cols; j++)
                data[i][j] = v[i][j];
        }
    }

    static Matrix<T> identity(size_t size, T value = 1) {
        Matrix<T> result(size);
        for (size_t i = 0; i < size; i++)
            result.set(i, i, value);
        return result;
    }

    static Matrix<T> zero(size_t rows, size_t cols, T value = 0) {
        Matrix<T> result(rows, cols);
        for (size_t i = 0; i < rows; i++)
            for (size_t j = 0; j < cols; j++)
                result.set(i, j, value);

        return result;
    }

    // Destructor
    ~Matrix() {
        if (data != nullptr) {
            for (size_t i = 0; i < rows; i++)
                delete[] data[i];
            delete[] data;
        }
        rows = 0;
        cols = 0;
        data = nullptr;
    }

    // Getters
    size_t getRows() const noexcept { return rows; }

    size_t getCols() const noexcept { return cols; }

    const T &get(size_t row, size_t col) const {
        if (!isValidIndex(row, col))
            throw std::out_of_range("Index out of range");
        return data[row][col];
    }

    // Setters
    void set(size_t row, size_t col, const T &value) {
        if (!isValidIndex(row, col))
            throw std::out_of_range("Index out of range");
        data[row][col] = value;
    }

    // Utils
    bool isValidIndex(size_t row, size_t col) const noexcept {
        return row < rows && col < cols;
    }

    Matrix<T> elementwise(const Matrix<T> &other, std::function<T(const T &, const T &)> operation) const {
        if (!sameShape(other))
            throw std::invalid_argument("Matrices must have the same shape");
        Matrix<T> result(rows, cols);
        for (size_t i = 0; i < rows; i++)
            for (size_t j = 0; j < cols; j++)
                result.set(i, j, operation(data[i][j], other.data[i][j]));
        return result;
    }

    bool sameShape(const Matrix<T> &other) const {
        return rows == other.rows && cols == other.cols;
    }

    std::array<size_t, 2> shape() const {
        return {rows, cols};
    }

    // Equality
    bool operator==(const Matrix<T> &other) const {
        if (!sameShape(other))
            return false;
        for (size_t i = 0; i < rows; i++)
            for (size_t j = 0; j < cols; j++)
                if (data[i][j] != other.data[i][j])
                    return false;
        return true;
    }

    bool operator!=(const Matrix<T> &other) const {
        return !(*this == other);
    }

    bool operator==(const T &value) const {
        return *this == Matrix<T>::identity(rows, value);
    }

    bool operator!=(const T &value) const {
        return !(*this == value);
    }

    // Assignment
    Matrix<T> &operator=(const Matrix<T> &other) {
        if (this == &other)
            return *this;
        else if (rows != other.rows || cols != other.cols) { // shape differs
            // delete old data
            for (size_t i = 0; i < rows; i++)
                delete[] data[i];
            delete[] data;

            // copy new data
            rows = other.rows;
            cols = other.cols;
            data = new T *[rows];
            for (size_t i = 0; i < rows; i++) {
                data[i] = new T[cols];
                for (size_t j = 0; j < cols; j++)
                    data[i][j] = other.data[i][j];
            }
        } else {
            // copy data
            for (size_t i = 0; i < rows; i++)
                for (size_t j = 0; j < cols; j++)
                    data[i][j] = other.data[i][j];
        }
        return *this;
    }

    // Arithmetic
    Matrix<T> operator+(const Matrix<T> &other) const {
        return elementwise(other, [](T a, T b) { return a + b; });
    }

    Matrix<T> operator-(const Matrix<T> &other) const {
        return elementwise(other, [](T a, T b) { return a - b; });
    }

    Matrix<T> operator*(const Matrix<T> &other) const {
        return elementwise(other, [](T a, T b) { return a * b; });
    }

    Matrix<T> operator*(const T &value) const {
        return elementwise(*this, [value](T a, T b) { return a * value; });
    }

    Matrix<T> operator/(const T &value) const {
        return elementwise(*this, [value](T a, T b) { return a / value; });
    }

    Matrix<T> operator-() const {
        return elementwise(*this, [](T a, T b) { return -a; });
    }

    // Transpose
    Matrix<T> transpose() const {
        Matrix<T> result(cols, rows);
        for (size_t i = 0; i < rows; i++)
            for (size_t j = 0; j < cols; j++)
                result.set(j, i, data[i][j]);
        return result;
    }

    // Elementary row operations
    void swapRows(size_t row1, size_t row2) {
        if (row1 == row2)
            return;
        if (row1 >= rows || row2 >= rows)
            throw std::out_of_range("Index out of range");
        T *temp = data[row1];
        data[row1] = data[row2];
        data[row2] = temp;
    }

    // TODO: check if value is zero for arithmetic types
    void multiplyRow(size_t row, const T &value) {
        if (row >= rows)
            throw std::out_of_range("Index out of range");
        for (size_t j = 0; j < cols; j++)
            data[row][j] *= value;
    }

    void addRow(size_t row1, size_t row2, const T &value) {
        if (row1 >= rows || row2 >= rows)
            throw std::out_of_range("Index out of range");
        if (row1 == row2)
            throw std::invalid_argument("Row indices must be different");

        for (size_t j = 0; j < cols; j++)
            data[row1][j] += data[row2][j] * value;
    }

    // Determinant
    // TODO: implement
    T determinant() const {

    }

    // Inverse
    // TODO: implement
    Matrix<T> operator!() const;
};


// Input and Output
// TODO: implement
template<typename T>
std::istream &operator>>(std::istream &in, Matrix<T> &mtx);

template<typename T>
std::ostream &operator<<(std::ostream &out, const Matrix<T> &mtx) {
    std::array<size_t, 2> shape = mtx.shape();
    out << shape[0] << " " << shape[1] << std::endl;
    for (size_t i = 0; i < shape[0]; i++) {
        for (size_t j = 0; j < shape[1]; j++)
            out << mtx.get(i, j) << " ";
        out << std::endl;
    }
    return out;
}

// LHS multiplication
template<typename T>
Matrix<T> operator*(const T &value, const Matrix<T> &mtx) {
    return mtx.elementwise(mtx, [value](T a, T b) { return value * a; });
}


#endif //INC_23_HOME_2_MATRIX_HPP
