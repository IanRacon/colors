#include "procedures.h"
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <cmath>
#include <sstream>

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
    std::vector<double> result;
    result.reserve(myVector.size());
    for (int i = 0; i < csrMatrix.rows(); ++i)
        result.push_back(0.0);
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
    std::vector<double> result;
    result.reserve(minuend.size());
    for (int i = 0; i < minuend.size(); ++i)
        result.push_back(minuend[i] - subtrahend[i]);
    return result;
}
void multiplyFast(const double multiplier, const double *vector, double *result, int size)
{
    for (int i = 0; i < size; ++i)
        result[i] = multiplier * vector[i];
}
std::vector<double> multiply(double multiplier, const std::vector<double> &vector)
{
    std::vector<double> result;
    result.reserve(vector.size());
    for (auto el : vector)
        result.push_back(multiplier * el);
    return result;
}
void addFast(const double *lhs, const double *rhs, double *result, int size)
{
    for (int i = 0; i < size; ++i)
        result[i] = lhs[i] + rhs[i];
}
std::vector<double> add(const std::vector<double> &lhs, const std::vector<double> &rhs)
{
    std::vector<double> result;
    result.reserve(lhs.size());
    for (int i = 0; i < lhs.size(); ++i)
        result.push_back(lhs[i] + rhs[i]);
    return result;
}
std::string printVector(const std::vector<double> &vector)
{
    std::stringstream ss;
    for (auto el : vector)
        ss << el << ", ";
    return ss.str();
}
std::vector<double> conjugateGradient(const CSR &csrMatrix,
                                      const std::vector<double> &b,
                                      const std::vector<double> &x0)
{
    std::vector<double> rj = substract(b, product(csrMatrix, x0));
    std::vector<double> pj = rj;
    std::vector<double> xj = x0;
    std::vector<double> xj1;
    std::vector<double> rj1;
    std::vector<double> pj1;
    double convergence = 1000;
    double alphaj = 0;
    double betaj = 0;
    const double limit = 1e-3;
    double rjrj = 0;
    std::vector<double> Apj;
    double alphajApj = 0;
    while (convergence > limit)
    {
        rjrj = std::inner_product(rj.begin(), rj.end(), rj.begin(), 0.0);
        Apj = product(csrMatrix, pj);
        alphaj = rjrj / std::inner_product(Apj.begin(), Apj.end(), pj.begin(), 0.0);
        xj1 = add(xj, multiply(alphaj, pj));
        rj1 = substract(rj, multiply(alphaj, Apj));
        betaj = std::inner_product(rj1.begin(), rj1.end(), rj1.begin(), 0.0) / rjrj;
        pj1 = add(rj1, multiply(betaj, pj));

        convergence = sqrt(rjrj);
        //std::cout << "Convergence: " << convergence << std::endl;
        // std::cout << "pj: " << printVector(pj) << ", pj1: " << printVector(pj1) << std::endl;
        // std::cout << "rj: " << printVector(rj) << ", rj1: " << printVector(rj1) << std::endl;
        // std::cout << "xj: " << printVector(xj) << ", xj1: " << printVector(xj1) << std::endl;

        pj = pj1;
        rj = rj1;
        xj = xj1;
    }
    return xj;
}
void conjugateGradientFast(const CSR &csrMatrix,
                           const std::vector<double> &barg,
                           const std::vector<double> &x0,
                           double *result)
{
    const int size = csrMatrix.rows();
    double *temp = new double[size];
    double *rj = new double[size];
    double *x0n = new double[size];
    double *b = new double[size];
    copy(b, barg);
    copy(x0n, x0);
    productFast(csrMatrix, x0n, temp);
    substractFast(b, temp, rj, size);
    double *pj = new double[size];
    pj = rj;
    double *xj = new double[size];
    copy(xj, x0);
    double *xj1 = new double[size];
    double *rj1 = new double[size];
    double *pj1 = new double[size];
    double convergence = 1000;
    double alphaj = 0;
    double betaj = 0;
    const double limit = 1e-3;
    double rjrj = 0;
    double *Apj = new double[size];
    double alphajApj = 0;
    while (convergence > limit)
    {
        rjrj = std::inner_product(rj, rj + size, rj, 0.0);
        productFast(csrMatrix, pj, Apj);
        alphaj = rjrj / std::inner_product(Apj, Apj + size, pj, 0.0);
        multiplyFast(alphaj, pj, temp, size);
        addFast(xj, temp, xj1, size);
        multiplyFast(alphaj, Apj, temp, size);
        substractFast(rj, temp, rj1, size);
        betaj = std::inner_product(rj1, rj1 + size, rj1, 0.0) / rjrj;
        multiplyFast(betaj, pj, temp, size);
        addFast(rj1, temp, pj1, size);

        convergence = sqrt(std::inner_product(rj, rj + size, rj, 0.0));
        //std::cout << "Convergence: " << convergence << std::endl;
        // std::cout << "pj: " << printVector(pj) << ", pj1: " << printVector(pj1) << std::endl;
        // std::cout << "rj: " << printVector(rj) << ", rj1: " << printVector(rj1) << std::endl;
        // std::cout << "xj: " << printVector(xj) << ", xj1: " << printVector(xj1) << std::endl;

        pj = pj1;
        rj = rj1;
        xj = xj1;
    }
}
}