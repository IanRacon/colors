#include <memory>
#include "gtest/gtest.h"
#include "function.h"

namespace function
{
TEST(functionCalculateAngle, calculateAngle_shouldReturnProperAngle)
{
    double angle = calculateAngle(10, 10, sqrt(100 + 100));
    ASSERT_NEAR(angle, M_PI / 4, 0.001);

    angle = calculateAngle(-10, 10, sqrt(100 + 100));
    ASSERT_NEAR(angle, 3 * M_PI / 4, 0.001);

    angle = calculateAngle(-10, -10, sqrt(100 + 100));
    ASSERT_NEAR(angle, 5 * M_PI / 4, 0.001);

    angle = calculateAngle(10, -10, sqrt(100 + 100));
    ASSERT_NEAR(angle, 7 * M_PI / 4, 0.001);
}
TEST(functionCalculateAngle, calculateAngle_shouldReturnProperAngle_forAxisPoints)
{
    double angle = calculateAngle(1, 0, sqrt(1));
    ASSERT_NEAR(angle, 0, 0.001);

    angle = calculateAngle(-1, 0, sqrt(1));
    ASSERT_NEAR(angle, M_PI, 0.001);

    angle = calculateAngle(0, -1, sqrt(1));
    ASSERT_NEAR(angle, 3 * M_PI / 2, 0.001);

    angle = calculateAngle(0, 1, sqrt(1));
    ASSERT_NEAR(angle, M_PI / 2, 0.001);
}
TEST(functionCalculateAngle, calculateAngle_shouldReturnProperAngle_forCenterPoint)
{
    double angle = calculateAngle(0, 0, sqrt(0));
    ASSERT_NEAR(angle, 0, 0.001);
}
}