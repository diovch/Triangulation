#ifndef BYTE_MATRIX_H
#define BYTE_MATRIX_H

#include <vector>
#include <stdexcept>
#include <cmath>
#include "R3Graph.h"

template <class T> class Matrix {
public:
    int m;      // Number of rows
    int n;      // Number of columns
    std::vector<T> a;

    Matrix(int r = 1, int c = 1):
        m(r),
        n(c),
        a(m*n)
    {
        for (int i = 0; i < m*n; ++i)
            a[i] = T(0);
    }

    int rows() const { return m; }
    int cols() const { return n; }
    std::vector<T>& dataArray() { return a; }
    const std::vector<T>& dataArray() const { return a; }

    void resize(int numRows, int numCols) {
        m = numRows;
        n = numCols;
        a.resize(m*n);
        zero();
    }

    void zero() {
        for (int i = 0; i < m*n; ++i)
            a[i] = T(0);
    }

    void fill(const T& v = T()) {
        for (int i = 0; i < m*n; ++i)
            a[i] = v;
    }

    T* operator[](int i) {
        return &(a[i*n]);
    }

    const T* operator[](int i) const {
        return &(a[i*n]);
    }

    T& at(int i, int j) {
        if (i < 0 || i > m || j < 0 || j > n)
            throw std::range_error("Matrix indices are out of range");
        return a[i*n + j];
    }

    const T& at(int i, int j) const {
        if (i < 0 || i > m || j < 0 || j > n)
            throw std::range_error("Matrix indices are out of range");
        return a[i*n + j];
    }
};

typedef class Matrix<unsigned char> ByteMatrix;
typedef class Matrix<short> ShortMatrix;
typedef class Matrix<int> IntMatrix;
typedef class Matrix<double> DoubleMatrix;

class R3Matrix: public DoubleMatrix {
public:
    R3Matrix():
        DoubleMatrix(3, 3)
    {}

    void setUnit() {
        zero();
        at(0, 0) = 1.;
        at(1, 1) = 1.;
        at(2, 2) = 1.;
    }

    static R3Matrix unit() {
        R3Matrix res;
        res.setUnit();
        return res;
    }

    static R3Matrix rotationZ(double alpha) {
        R3Matrix res;
        double sa = sin(alpha);
        double ca = cos(alpha);
        res[0][0] = ca; res[0][1] = (-sa);
        res[1][0] = sa; res[1][1] = ca;
        res[2][2] = 1.;
        return res;
    }

    static R3Matrix rotationX(double alpha) {
        R3Matrix res;
        double sa = sin(alpha);
        double ca = cos(alpha);
        res[0][0] = 1.;
        res[1][0] = ca; res[1][1] = (-sa);
        res[2][0] = sa; res[2][1] = ca;
        return res;
    }

    static R3Matrix rotationY(double alpha) {
        R3Matrix res;
        double sa = sin(alpha);
        double ca = cos(alpha);
        res[0][0] = ca; res[2][1] = sa;
        res[1][1] = 1.;
        res[2][0] = (-sa); res[2][2] = ca;
        return res;
    }

    R3Matrix operator*(const R3Matrix& b) const;

    double det() const;
    R3Matrix inverse() const;

    static R3Matrix rotation(const R3Graph::R3Vector& axis, double alpha);
};

R3Graph::R3Vector operator*(const R3Matrix& a, const R3Graph::R3Vector& v);

#endif
