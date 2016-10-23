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
}