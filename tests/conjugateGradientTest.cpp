#include <memory>
#include "gtest/gtest.h"
#include "procedures.h"
#include "CSR.h"

namespace procedures
{
TEST(procedureConjugateGradient, conjugateGradient_shouldReturnRightVector_for3x3)
{
    int size = 3;
    int csrValues[9] = {2, 17, 3, 17, 8, 0, 3, 0, 15};
    std::vector<double> b = {2, 7, 38};
    std::vector<double> x0 = {0, 0, 0};
    std::vector<double> expected = {0.5896, -0.3780, 2.4154};

    CSR csrMatrix(size, size, 9);
    for (int i = 0; i < size; ++i)
        for (int j = 0; j < size; ++j)
            csrMatrix.setValue(i, j, csrValues[i * size + j]);
    csrMatrix.setEndIndicator();

    std::vector<double> result = conjugateGradient(csrMatrix, b, x0);
    for (int i = 0; i < size; ++i)
    {
        ASSERT_NEAR(expected.at(i), result.at(i), 0.0001);
    }
}
TEST(procedureConjugateGradient, conjugateGradient_shouldReturnRightVector_for6x6)
{
    int size = 6;
    std::vector<double> b = {1, 2, 3, 4, 5, 6};
    std::vector<double> x0 = {0, 0, 0, 0, 0, 0};
    std::vector<double> expected = {0.6667, 0.4667, 0.2667, 0.0667, -0.1333, -0.3333};
    CSR csrMatrix(size, size, 6 * 6);
    for (int i = 0; i < size; ++i)
        for (int j = 0; j < size; ++j)
            csrMatrix.setValue(i, j, (i + 1) + (j + 1));
    csrMatrix.setEndIndicator();

    std::vector<double> result = conjugateGradient(csrMatrix, b, x0);
    for (int i = 0; i < size; ++i)
    {
        ASSERT_NEAR(expected.at(i), result.at(i), 0.0001);
    }
}

TEST(procedureConjugateGradient, conjugateGradient_shouldReturnRightVector_for1000x1000)
{
    int size = 1000;
    std::vector<double> b;
    b.reserve(size);
    for (int i = 0; i < size; ++i)
        b.push_back(i);
    std::vector<double> x0;
    x0.reserve(size);
    for (int i = 0; i < size; ++i)
        x0.push_back(0);

    CSR csrMatrix(size, size, size * size);
    for (int i = 0; i < size; ++i)
        for (int j = 0; j < size; ++j)
            csrMatrix.setValue(i, j, (i + 1) + (j + 1));
    csrMatrix.setEndIndicator();
    double *result = new double[size];
    conjugateGradientFast(csrMatrix, b, x0, result);
    //std::vector<double> result = conjugateGradient(csrMatrix, b, x0);
}
}
