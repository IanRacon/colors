#include "gtest/gtest.h"
#include "imnmath.hpp"
#include "function.h"

namespace function
{
TEST(functionfillVDistribTest, fillVDistrib_shouldFill0InCenterPoints)
{
    int size = 100;
    double speedFactor = 2.0;
    double range = 10.0;
    double **velocityXMatrix = fillVDistrib(50, 50, range, speedFactor, size, clockwiseX);
    ASSERT_EQ(0.0, velocityXMatrix[50][50]);
    imn<double>::free_matrix(velocityXMatrix, size);
}
TEST(functionfillVDistribTest, fillVDistrib_shouldNotFillBeyondRange)
{
    int size = 100;
    double speedFactor = 2.0;
    double range = 10.0;
    double **velocityMatrix = fillVDistrib(50, 50, range, speedFactor, size, clockwiseX);
    ASSERT_EQ(0.0, velocityMatrix[50][61]);
    imn<double>::free_matrix(velocityMatrix, size);
}
TEST(functionfillVDistribTest, fillVDistrib_shouldFillInRange)
{
    int size = 100;
    double speedFactor = 2.0;
    double range = 10.0;
    double **velocityMatrix = fillVDistrib(49, 49, range, speedFactor, size, clockwiseX);
    ASSERT_TRUE(velocityMatrix[59][49] > 0.0);
    imn<double>::free_matrix(velocityMatrix, size);
}
TEST(functionfillVDistribTest, fillVDistribX_shouldBeZeroOnXAxis)
{
    int size = 100;
    double speedFactor = 2.0;
    double range = 10.0;
    double **velocityMatrix = fillVDistrib(49, 49, range, speedFactor, size, clockwiseX);
    for (int i = 0; i < (int)range; ++i)
    {
        ASSERT_NEAR(0.0, velocityMatrix[49][49 + i], 1e-12);
        ASSERT_NEAR(0.0, velocityMatrix[49][49 - i], 1e-12);
    }
    imn<double>::free_matrix(velocityMatrix, size);
}
TEST(functionfillVDistribTest, fillVDistribY_shouldBeZeroOnYAxis)
{
    int size = 100;
    double speedFactor = 2.0;
    double range = 10.0;
    double **velocityMatrix = fillVDistrib(49, 49, range, speedFactor, size, clockwiseY);
    for (int i = 0; i < (int)range; ++i)
    {
        ASSERT_NEAR(0.0, velocityMatrix[49 + i][49], 1e-12);
        ASSERT_NEAR(0.0, velocityMatrix[49 - i][49], 1e-12);
    }
    imn<double>::free_matrix(velocityMatrix, size);
}
}