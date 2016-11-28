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
    std::vector<double> rj = substract(b, product(A, x0));
    std::vector<double> pj = rj;
    std::vector<double> x_1 = x0;
    std::vector<double> x_11;
    std::vector<double> rj1;
    double convergence = 1000;
    double alphaj = 0;
    double betaj = 0;
    const double limit = 1e-9;
    double rjrj = 0;
    int counter = 0;
    std::vector<double> Apj;
    while (convergence > limit)
    {
        rjrj = std::inner_product(rj.begin(), rj.end(), rj.begin(), 0.0);
        Apj = product(A, pj);
        alphaj = rjrj / std::inner_product(Apj.begin(), Apj.end(), pj.begin(), 0.0);
        x_11 = add(x_1, multiply(alphaj, pj));
        rj1 = substract(rj, multiply(alphaj, Apj));
        betaj = std::inner_product(rj1.begin(), rj1.end(), rj1.begin(), 0.0) / rjrj;
        pj = add(rj1, multiply(betaj, pj));

        convergence = sqrt(rjrj);
        //std::cout << "Convergence: " << convergence << std::endl;
        // std::cout << "pj: " << printVector(pj) << ", pj: " << printVector(pj) << std::endl;
        // std::cout << "rj: " << printVector(rj) << ", rj1: " << printVector(rj1) << std::endl;
        // std::cout << "x_1: " << printVector(x_1) << ", x_11: " << printVector(x_11) << std::endl;

        pj = pj;
        rj = rj1;
        x_1 = x_11;
        counter++;
    }
    //std::cout << "Iteracji sprzężonego gradientu: " << counter << std::endl;
    return x_1;
}
std::vector<double> bicgstab(const CSR &A, const std::vector<double> &b, const std::vector<double> &x0)
{
    // for (int i = 0; i < b.size(); ++i)
    // {
    //     std::cout << b[i] << " ";
    // }
    std::vector<double> r_1 = substract(b, product(A, x0));
    //std::vector<double> pj = rj;
    std::vector<double> rj0p = r_1;
    std::vector<double> x_1 = x0;
    std::vector<double> x;
    std::vector<double> r;
    std::vector<double> s;
    std::vector<double> t;
    std::vector<double> p = r;
    std::vector<double> p_1(p.size());
    std::vector<double> h;
    std::vector<double> v;
    const double limit = 1e-9;
    double bnorm2 = norm(b);
    if (bnorm2 = 0.0)
        bnorm2 = 1.0;
    double error = norm(r) / bnorm2;
    if (error < limit)
    {
        std::cout << "error < limit " << std::endl;
        return x_1;
    }

    double wj = 0;
    std::vector<double> v_1(b.size());
    double rho = 1;
    double rho_1 = 1;
    double alpha = 1;
    double omega_1 = 1;
    double omega = 1;
    double beta = 0;
    std::vector<double> alphapj;

    double rjrj = 0;
    int counter = 0;
    std::vector<double> Apj;
    std::vector<double> Asj;
    while (error > limit)
    {
        rho = dot(r, rj0p);
        if (rho == 0.0)
        {
            std::cout << "rho = 0.0 " << std::endl;
            break;
        }
        beta = (rho / rho_1) * (alpha / omega_1);
        p = add(r_1, multiply(beta, substract(p_1, multiply(omega_1, v_1))));

        v = product(A, p);
        alpha = rho / dot(rj0p, v);
        x = add(x_1, multiply(alpha, p));
        if (norm(x) < error)
            return x;
        s = substract(r_1, multiply(alpha, v));
        t = product(A, s);
        omega = dot(t, s) / dot(t, t);
        x = add(h, multiply(omega, s));
        if (norm(x) < error)
            return x;
        r = substract(s, multiply(omega, t));

        omega_1 = omega;
        rho_1 = rho;
        r_1 = r;
        p_1 = p;
        x_1 = x;
        v_1 = v;
        counter++;
        // Apj = product(A, pj);
        // alphaj = pi / std::inner_product(Apj.begin(), Apj.end(), rj0p.begin(), 0.0);
        // sj = substract(rj, multiply(alphaj, Apj));
        // Asj = product(A, sj);
        // wj = std::inner_product(Asj.begin(), Asj.end(), sj.begin(), 0.0) /
        //      std::inner_product(Asj.begin(), Asj.end(), Asj.begin(), 0.0);
        // alphajpj = multiply(alphaj, pj);
        // x_11 = add(add(x_1, alphajpj), multiply(wj, sj));
        // rj1 = substract(sj, multiply(wj, Asj));
        // betaj = std::inner_product(rj1.begin(), rj1.end(), rj0p.begin(), 0.0) /
        //         pi * (alphaj / wj);
        // pj = add(rj1, multiply(betaj, substract(pj, multiply(wj, Apj))));

        // convergence = sqrt(rjrj);
        // std::cout << "Convergence: " << convergence << std::endl;
        // std::cout << "pj: " << printVector(pj) << ", pj: " << printVector(pj) << std::endl;
        // std::cout << "rj: " << printVector(rj) << ", rj1: " << printVector(rj1) << std::endl;
        // std::cout << "x_1: " << printVector(x_1) << ", x_11: " << printVector(x_11) << std::endl;

        // pi = pi1;
        // pj = pj;
        // rj = rj1;
        // x_1 = x_11;
        // counter++;
    }
    std::cout << "Iteracji bicgstab : " << counter << std::endl;
    return x_1;
}
}