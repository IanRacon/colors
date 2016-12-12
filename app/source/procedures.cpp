#include "procedures.h"
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <cmath>
#include <sstream>
#include <future>
#include <time.h>
#include <chrono>

namespace procedures
{
void copy(double *lhs, const std::vector<double> &rhs)
{
    for (int i = 0; i < rhs.size(); ++i)
        lhs[i] = rhs[i];
}
void productFast(const CSR &csrMatrix, const double *vector, double *result)
{
    for (int i = 0; i < csrMatrix.rows(); ++i)
    {
        result[i] = 0;
        for (int j = csrMatrix.rowStartIndices[i]; j < csrMatrix.rowStartIndices[i + 1] - 1; ++j)
            result[i] += csrMatrix.allElements[j] * vector[csrMatrix.byColumnIndices[j]];
    }
}
std::vector<double> product(const CSR &csrMatrix, const std::vector<double> &myVector)
{
    std::vector<double> result(myVector.size());
    //result.reserve(myVector.size());
    // for (int i = 0; i < csrMatrix.rows(); ++i)
    //     result.push_back(0.0);
    if (csrMatrix.allElements.size() == 0)
        return result;
    if (csrMatrix.cols() != myVector.size())
        throw std::logic_error("Size of matrix and vector not the same");

    for (int i = 0; i < csrMatrix.rows(); ++i)
        for (int j = csrMatrix.rowStartIndices[i]; j < csrMatrix.rowStartIndices[i + 1]; ++j)
            result[i] += csrMatrix.allElements[j] * myVector[csrMatrix.byColumnIndices[j]];
    return result;
}
void substractFast(const double *lhs, const double *rhs, double *result, int size)
{
    for (int i = 0; i < size; ++i)
        result[i] = lhs[i] - rhs[i];
}
std::vector<double> substract(const std::vector<double> &minuend, const std::vector<double> &subtrahend)
{
    std::vector<double> result(minuend.size());
    for (int i = 0; i < minuend.size(); ++i)
        result[i] = minuend[i] - subtrahend[i];
    return result;
}
void multiplyFast(const double multiplier, const double *vector, double *result, int size)
{
    for (int i = 0; i < size; ++i)
        result[i] = multiplier * vector[i];
}
std::vector<double> multiply(double multiplier, const std::vector<double> &vector)
{
    std::vector<double> result(vector.size());
    for (int i = 0; i < vector.size(); ++i)
        result[i] = multiplier * vector[i];
    return result;
}
void addFast(const double *lhs, const double *rhs, double *result, int size)
{
    for (int i = 0; i < size; ++i)
        result[i] = lhs[i] + rhs[i];
}
std::vector<double> add(const std::vector<double> &lhs, const std::vector<double> &rhs)
{
    std::vector<double> result(lhs.size());
    for (int i = 0; i < lhs.size(); ++i)
        result[i] = lhs[i] + rhs[i];
    return result;
}
std::string printVector(const std::vector<double> &vector)
{
    std::stringstream ss;
    for (auto el : vector)
        ss << el << ", ";
    return ss.str();
}
double sum2DArray(const double *const *const array, double rows, double cols)
{
    double sum = 0.0;
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            sum += array[i][j];
    return sum;
}
double dot(const std::vector<double> &lhs, const std::vector<double> &rhs)
{
    return std::inner_product(lhs.begin(), lhs.end(), rhs.begin(), 0.0);
}
double norm(const std::vector<double> &vec)
{
    return sqrt(dot(vec, vec));
}
std::vector<double> conjugateGradient(const CSR &A,
                                      const std::vector<double> &b,
                                      const std::vector<double> &x0)
{
    std::vector<double> r = substract(b, product(A, x0));
    std::vector<double> p = r;
    std::vector<double> x_1 = x0;
    std::vector<double> x_11, r1;
    double convergence = 1000;
    double alpha, beta, rr;
    const double limit = 1e-9;
    int i = 0;
    std::vector<double> Ap;
    while (convergence > limit)
    {
        rr = dot(r, r);
        Ap = product(A, p);
        alpha = rr / dot(Ap, p);
        x_11 = add(x_1, multiply(alpha, p));
        r1 = substract(r, multiply(alpha, Ap));
        beta = dot(r1, r1) / rr;
        p = add(r1, multiply(beta, p));

        convergence = sqrt(rr);
        p = p;
        r = r1;
        x_1 = x_11;
        i++;
    }
    return x_1;
}
std::vector<double> bicgstab(const CSR &A, const std::vector<double> &b, const std::vector<double> &x0)
{
    std::vector<double> r = substract(b, product(A, x0));
    std::vector<double> p = r;
    std::vector<double> r0(r.size());
    for (int i = 0; i < r.size(); ++i)
        r0[i] = (double)i / (double)r.size();
    std::vector<double> x = x0;
    std::vector<double> s, Ap, As;
    double error = 1, alpha = 1, omega = 1, beta = 1, rr0 = 1;
    const double limit = 1e-16;
    int i = 0;
    error = 1;
    while (error > limit && i < 2000)
    {
        Ap = product(A, p);
        rr0 = dot(r, r0);
        alpha = rr0 / dot(Ap, r0);
        s = substract(r, multiply(alpha, Ap));
        As = product(A, s);
        omega = dot(As, s) / dot(As, As);
        x = add(add(x, multiply(alpha, p)), multiply(omega, s));
        r = substract(s, multiply(omega, As));
        beta = dot(r, r0) / rr0 * (alpha / omega);
        p = add(r, multiply(beta, substract(p, multiply(omega, Ap))));
        error = norm(r);
        ++i;
    }
    //std::cout << "Iteracji bicgstab : " << i << std::endl;
    return x;
}
std::vector<double> getInitialDensity(int centerX, int centerY, int rows, int cols, double moveStep, double sigma)
{
    std::vector<double> initialDensity(rows * cols);
    double x, y, r;
    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j)
        {
            x = moveStep * (i - centerY);
            y = moveStep * (j - centerX);
            r = x * x + y * y;
            initialDensity[i * cols + j] = exp(-r * pow(sigma, 2) / 2.0) / M_PI;
        }
    }
    return initialDensity;
}
}