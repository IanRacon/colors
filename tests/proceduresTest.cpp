#include "gtest/gtest.h"
#include "procedures.h"
#include "CSR.h"

namespace procedures
{
CSR setCSRmatrixSingleValue(int size, double value)
{
    CSR csrMatrix(size, size, size*size);
    for(int i=0;i<size;++i)
        for(int j=0;j<size;++j)
            csrMatrix.setValue(i, j, value);
    csrMatrix.setEndIndicator();
    return csrMatrix;
}
std::vector<double> createSingleValuedVector(int size, double value)
{
    std::vector<double> vec;
    vec.reserve(size);
    for(int i=0;i<size;++i)
        vec.push_back(value);
    return vec;
}
TEST(procedureProduct, product_shouldReturnValidProductWithSameValues)
{
    CSR csrMatrix = setCSRmatrixSingleValue(3, 2.0);
    std::vector<double> givenVec = createSingleValuedVector(3, 2.0);
    std::vector<double> resultVec = createSingleValuedVector(3, 12.0);
    ASSERT_EQ(resultVec, product(csrMatrix, givenVec));
}
TEST(procedureProduct, product_shouldReturnValidProductWithZeroMatrix)
{
    CSR csrMatrix = setCSRmatrixSingleValue(3, 0.0);
    std::vector<double> givenVec = createSingleValuedVector(3, 2.0);
    std::vector<double> resultVec = createSingleValuedVector(3, 0.0);
    ASSERT_EQ(resultVec, product(csrMatrix, givenVec));
}
TEST(procedureProduct, product_shouldReturnValidProductWithIdentityMatrix)
{
    CSR csrMatrix(3, 3, 3);
    for(int i=0;i<3;++i)
            csrMatrix.setValue(i, i, 1.0);
    csrMatrix.setEndIndicator();
    std::vector<double> givenVec = createSingleValuedVector(3, 2.0);
    std::vector<double> resultVec = givenVec;
    ASSERT_EQ(resultVec, product(csrMatrix, givenVec));
}
TEST(procedureProduct, product_shouldReturnValidProductWithMinusIdentityMatrix)
{
    CSR csrMatrix(3, 3, 3);
    for(int i=0;i<3;++i)
            csrMatrix.setValue(i, i, -1.0);
    csrMatrix.setEndIndicator();
    std::vector<double> givenVec = createSingleValuedVector(3, 2.0);
    std::vector<double> resultVec = createSingleValuedVector(3, -2.0);
    ASSERT_EQ(resultVec, product(csrMatrix, givenVec));
}
TEST(procedureProduct, product_shouldReturnValidProductWithZeroVector)
{
    CSR csrMatrix(3, 3, 3);
    for(int i=0;i<3;++i)
            csrMatrix.setValue(i, i, 1.0);
    csrMatrix.setEndIndicator();
    std::vector<double> givenVec = createSingleValuedVector(3, 0.0);
    std::vector<double> resultVec = givenVec;
    ASSERT_EQ(resultVec, product(csrMatrix, givenVec));
}
TEST(procedureProduct, product_shouldReturnValidProductWithMatrixOneRowZero)
{
    CSR csrMatrix(3, 3, 3);
    csrMatrix.setValue(0, 0, 1.0);
    csrMatrix.setValue(1, 1, 0.0);
    csrMatrix.setValue(2, 2, 1.0);
    csrMatrix.setEndIndicator();
    std::vector<double> givenVec = createSingleValuedVector(3, 2.0);
    std::vector<double> resultVec = {2, 0, 2};
    ASSERT_EQ(resultVec, product(csrMatrix, givenVec));
}
TEST(procedureProduct, product_shouldReturnValidProductWithMatrixScrambledValues)
{
    int values[9] = {1, 2, 3, 1, 2, 3, 1, 2, 3};
    CSR csrMatrix(3, 3, 3);
    for(int i=0;i<3;++i)
        for(int j=0;j<3;++j)
            csrMatrix.setValue(i, j, values[i*3+j]);
    csrMatrix.setEndIndicator();
    std::vector<double> givenVec = {10, 20, 30};
    std::vector<double> resultVec = {140, 140, 140};
    ASSERT_EQ(resultVec, product(csrMatrix, givenVec));
}
TEST(procedureProduct, product_shouldReturnValidProductWithMatrixFirstColumnsEmpty)
{
    int values[9] = {0, 0, 3, 0, 0, 0, 0, 0, 3};
    CSR csrMatrix(3, 3, 3);
    for(int i=0;i<3;++i)
        for(int j=0;j<3;++j)
            csrMatrix.setValue(i, j, values[i*3+j]);
    csrMatrix.setEndIndicator();
    std::vector<double> givenVec = {4, 5, 6};
    std::vector<double> resultVec = {18, 0, 18};
    ASSERT_EQ(resultVec, product(csrMatrix, givenVec));
}
TEST(procedureProduct, product_shouldReturnValidProductWithMatrixLastColumnsEmpty)
{
    int values[9] = {3, 0, 0, 0, 0, 0, 3, 0, 0};
    CSR csrMatrix(3, 3, 3);
    for(int i=0;i<3;++i)
        for(int j=0;j<3;++j)
            csrMatrix.setValue(i, j, values[i*3+j]);
    csrMatrix.setEndIndicator();
    std::vector<double> givenVec = {4, 5, 6};
    std::vector<double> resultVec = {12, 0, 12};
    ASSERT_EQ(resultVec, product(csrMatrix, givenVec));
}
}