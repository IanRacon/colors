#include "gtest/gtest.h"
#include "numerical.h"

namespace numerical
{
TEST(calculateModulusTest, modulusWithShouldBeValid)
{
    ASSERT_EQ(2.5, calculateModulus(1, 1, 10));
}
TEST(calculateModulusTest, modulusWithZeroVelocityShouldBeZero)
{
    ASSERT_EQ(0.0, calculateModulus(1, 1, 0));
}
}