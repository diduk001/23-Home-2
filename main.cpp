//
// Created by Stepan Didurenko on 20.04.2024.
//

#include "Matrix.hpp"
#include <iostream>

int main() {
    Matrix<int> m(3, 3);
    m.set(0, 0, 1);
    m.set(0, 1, 2);
    m.set(0, 2, 3);
    m.set(1, 0, 4);
    m.set(1, 1, 5);
    m.set(1, 2, 6);
    m.set(2, 0, 7);
    m.set(2, 1, 8);
    m.set(2, 2, 9);
    m = m * 3;
    std::cout << m;
    return 0;
}