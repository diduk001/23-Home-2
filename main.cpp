#include "Matrix.hpp"
#include <iostream>
#include <fstream>

int main() {
    // mtx1 = {{8, 4, 5}, {8, 5, 7}, {0, 0, 7}}
    std::string inputFilename;
    std::cout << "Matrix input filename > ";
    std::cin >> inputFilename;
    std::ifstream in(inputFilename);

    Matrix<double> mtx1;
    in >> mtx1;

    Matrix mtx2 = Matrix(mtx1);
    mtx2.swapRows(0, 1);

    Matrix mtx3 = Matrix(mtx1);
    mtx3.multiplyRow(1, 4);

    Matrix mtx4 = Matrix(mtx1);
    mtx4.addRow(2, 0, 3);
    Matrix mtx5 = !mtx1;

    Matrix mtx6 = Matrix(mtx1).transpose();

    std::string outputFilename;
    std::cout << "Matrix output filename > ";
    std::cin >> outputFilename;
    std::ofstream out(outputFilename);

    out << "Original matrix (A):" << '\n';
    out << mtx1 << std::endl;

    out << "Swapped rows 0 and 1" << std::endl;
    out << mtx2 << std::endl;

    out << "Multiplied row 1 by 4" << std::endl;
    out << mtx3 << std::endl;

    out << "Added row 0 multiplied by 3 to row 2" << std::endl;
    out << mtx4 << std::endl;

    out << "Original matrix determinant: " << mtx1.det() << std::endl;
    out << "Inverse original matrix (I)" << std::endl;
    out << mtx5 << std::endl;

    out << "A * I = " << std::endl;
    out << mtx1 * mtx5 << std::endl;

    out << "I * A = " << std::endl;
    out << mtx5 * mtx1 << std::endl;

    out << "A^T" << std::endl;
    out << mtx6 << std::endl;

    return 0;
}