#include "gtest/gtest.h"
#include "procedures.h"
#include "imnmath.hpp"

namespace procedures
{
TEST(procedureSum2Darray, sumZeroValued2DShouldResultZero)
{
    double **zero2D = imn<double>::matrix(4, 3);
    imn<double>::set_matrix(zero2D, 4, 3, 0.0);
    ASSERT_EQ(0.0, sum2DArray(zero2D, 4, 3));
}
TEST(procedureSum2Darray, sumOneValued2DShouldBeValid)
{
    double **zero2D = imn<double>::matrix(100, 100);
    imn<double>::set_matrix(zero2D, 100, 100, 1.0);
    ASSERT_EQ(10000.0, sum2DArray(zero2D, 100, 100));
}
}