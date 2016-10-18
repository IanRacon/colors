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
TEST(procedureProduct, product_shouldReturnValidProduct)
{
    CSR csrMatrix = setCSRmatrixSingleValue(3, 2.0);
    std::vector<double> givenVec = createSingleValuedVector(3, 2.0);
    std::vector<double> resultVec = createSingleValuedVector(3, 12.0);
    ASSERT_EQ(resultVec, product(csrMatrix, givenVec));
}
}