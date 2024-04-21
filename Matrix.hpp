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

    Matrix(size_t newRows, size_t newCols, T **newElements) {
        if (newRows <= 0 || newCols <= 0)
            throw std::invalid_argument("Rows and newCols must be positive numbers");
        this->rows = newRows;
        this->cols = newCols;
        this->data = new T *[newRows];
        for (int rowIdx = 0; rowIdx < newRows; rowIdx++) {
            this->data[rowIdx] = new T[newCols];
            for (int colIdx = 0; colIdx < newCols; colIdx++)
                this->data[rowIdx][colIdx] = newElements[rowIdx][colIdx];
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

    void divideRow(size_t row, const T &value) {
        if (row >= rows)
            throw std::out_of_range("Index out of range");
        for (size_t j = 0; j < cols; j++)
            data[row][j] /= value;
    }

    void addRow(size_t row1, size_t row2, const T &value) {
        if (row1 >= rows || row2 >= rows)
            throw std::out_of_range("Index out of range");
        if (row1 == row2)
            throw std::invalid_argument("Row indices must be different");

        for (size_t j = 0; j < cols; j++)
            data[row1][j] += data[row2][j] * value;
    }

    // Find matrix determinant
    // Convert into diagonal matrix (like in Gauss method) and return product of all diagonal elements
    // https://cp-algorithms.com/linear_algebra/determinant-gauss.html
    T det() const {
        if (this->getRows() != this->getCols())
            throw std::invalid_argument("Matrix must be square");

        const size_t size = this->getCols();
        if (size == 2) {
            T a = this->get(0, 0),
                    b = this->get(0, 1),
                    c = this->get(1, 0),
                    d = this->get(1, 1);
            return a * d - b * c;
        }

        Matrix mtx(*this);
        T result = 1;
        for (size_t stepIdx = 0; stepIdx < size; stepIdx++) {
            // Choose row with maximal absolute value of stepIdx_th element
            size_t maxByAbsRowIdx = stepIdx;
            T maxByAbsValue = mtx.get(maxByAbsRowIdx, stepIdx);
            for (size_t maxByAbsRowCandidate = stepIdx + 1; maxByAbsRowCandidate < size; maxByAbsRowCandidate++) {
                T candidateValue = std::abs(mtx.get(maxByAbsRowCandidate, stepIdx));
                if (candidateValue > std::abs(maxByAbsValue)) {
                    maxByAbsRowIdx = maxByAbsRowCandidate;
                    maxByAbsValue = mtx.get(maxByAbsRowIdx, stepIdx);
                }
            }
            if (std::abs(maxByAbsValue) == 0)
                return 0;

            // Swap current step row and max abs row
            if (stepIdx != maxByAbsRowIdx)
                result = -result;
            mtx.swapRows(stepIdx, maxByAbsRowIdx);

            // Divide max abs row by first stepIdx element
            result *= maxByAbsValue;
            mtx.divideRow(stepIdx, maxByAbsValue);

            // Nullify all elements under [stepIdx][stepIdx]
            for (size_t lowerRowIdx = stepIdx + 1; lowerRowIdx < size; lowerRowIdx++) {
                if (std::abs(mtx.get(lowerRowIdx, stepIdx)) == 0)
                    continue;
                T coefficient = -mtx.get(lowerRowIdx, stepIdx);
                mtx.addRow(lowerRowIdx, stepIdx, coefficient);
            }
        }

        return result;
    }

    // Inverse matrix
    // Convert into identity matrix (like in Gauss method) and do same operation on identity matrix
    // https://cp-algorithms.com/linear_algebra/linear-system-gauss.html
    // Returns 0 matrix if matrix is not invertible
    Matrix operator!() const {
        if (this->getRows() != this->getCols())
            throw std::invalid_argument("Matrix must be square");

        const size_t size = this->getCols();
        const T det = this->det();
        if (det == 0)
            return {size, size};

        if (size == 2) {
            T a = this->get(0, 0),
                    b = this->get(0, 1),
                    c = this->get(1, 0),
                    d = this->get(1, 1);

            T **resultElements;
            resultElements = new T *[2];
            resultElements[0] = new T[2]{d, -b};
            resultElements[1] = new T[2]{-c, a};

            return Matrix(2, 2, data) / det;
        }

        Matrix mtx = Matrix(*this);
        Matrix result = Matrix::identity(size);

        for (size_t stepIdx = 0; stepIdx < size; stepIdx++) {
            // Choose row with maximal absolute value of stepIdx_th element
            size_t maxByAbsRowIdx = stepIdx;
            T maxByAbsValue = mtx.get(maxByAbsRowIdx, stepIdx);
            for (size_t maxByAbsRowCandidate = stepIdx + 1; maxByAbsRowCandidate < size; maxByAbsRowCandidate++) {
                T candidateValue = std::abs(mtx.get(maxByAbsRowCandidate, stepIdx));
                if (candidateValue > std::abs(maxByAbsValue)) {
                    maxByAbsRowIdx = maxByAbsRowCandidate;
                    maxByAbsValue = mtx.get(maxByAbsRowIdx, stepIdx);
                }
            }

            // Swap current step row and max abs row
            mtx.swapRows(stepIdx, maxByAbsRowIdx);
            result.swapRows(stepIdx, maxByAbsRowIdx);

            // Divide max abs row by first stepIdx element
            mtx.multiplyRow(stepIdx, 1 / maxByAbsValue);
            result.multiplyRow(stepIdx, 1 / maxByAbsValue);

            // Nullify all elements under [stepIdx][stepIdx]
            for (size_t lowerRowIdx = 0; lowerRowIdx < size; lowerRowIdx++) {
                if (lowerRowIdx == stepIdx || std::abs(mtx.get(lowerRowIdx, stepIdx)) == 0) {
                    continue;
                }
                double coefficient = -mtx.get(lowerRowIdx, stepIdx);
                mtx.addRow(lowerRowIdx, stepIdx, coefficient);
                result.addRow(lowerRowIdx, stepIdx, coefficient);
            }
        }

        return result;
    }

};


// Input and Output
template<typename T>
std::istream &operator>>(std::istream &in, Matrix<T> &mtx) {
    size_t rows, cols;
    in >> rows >> cols;
    mtx = Matrix<T>(rows, cols);
    for (size_t i = 0; i < mtx.getRows(); i++)
        for (size_t j = 0; j < mtx.getCols(); j++) {
            T value;
            in >> value;
            mtx.set(i, j, value);
        }
    return in;
}

template<typename T>
std::ostream &operator<<(std::ostream &out, const Matrix<T> &mtx) {
    const size_t rows = mtx.getRows();
    const size_t cols = mtx.getCols();

    out << "[" << std::endl;
    for (size_t rowIdx = 0; rowIdx < rows; rowIdx++) {
        out << "\t[";
        for (size_t colIdx = 0; colIdx < cols; colIdx++) {
            out << "\t" << mtx.get(rowIdx, colIdx);
            if (colIdx < cols - 1)
                out << ",";
        }
        out << "]";
        if (rowIdx < rows - 1)
            out << ",";
        out << std::endl;
    }
    out << "]" << std::endl;

    return out;
}

// LHS multiplication
template<typename T>
Matrix<T> operator*(const T &value, const Matrix<T> &mtx) {
    return mtx.elementwise(mtx, [value](T a, T b) { return value * a; });
}


#endif //INC_23_HOME_2_MATRIX_HPP
