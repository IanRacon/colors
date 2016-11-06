#pragma once
#include "CSR.h"
#include <string>

namespace procedures
{
std::string printVector(const std::vector<double> &vector);
std::vector<double> product(const CSR &csrMatrix, const std::vector<double> &vector);
std::vector<double> multiply(double multiplier, const std::vector<double> &vector);
std::vector<double> conjugateGradient(const CSR &csrMatrix, const std::vector<double> &b, const std::vector<double> &x0);
std::vector<double> substract(const std::vector<double> &lhs, const std::vector<double> &rhs);
std::vector<double> add(const std::vector<double> &lhs, const std::vector<double> &rhs);
void copy(double *lhs, const std::vector<double> &rhs);
// std::string printVector(const std::vector<double> &vector);
void productFast(const CSR &csrMatrix, const double *vector, double *result);
void multiplyFast(const double multiplier, const double *vector, double *result);
void conjugateGradientFast(const CSR &csrMatrix, const std::vector<double> &b, const std::vector<double> &x0, double *result);
void substractFast(const double *lhs, const double *rhs, double *result, int size);
void addFast(const double *lhs, const double *rhs, double *result, int size);
double sum2DArray(const double *const *const array, double rows, double cols);
}