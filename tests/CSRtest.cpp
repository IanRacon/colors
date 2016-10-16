#include "gtest/gtest.h"
#include "CSR.h"

class CSRTest : public ::testing::Test
{
public:
    void SetUp()
    {

    }
};

TEST_F(CSRTest, test_should_success)
{
    ASSERT_EQ(true, true);
}
TEST(myTest, test)
{
    ASSERT_EQ(1, 2);
}