#include "gtest/gtest.h"
#include "numerical.h"
#include "imnmath.hpp"

namespace numerical
{
const int size = 2;
const double spinCenterX = 2;
const double spinCenterY = 2;
const double range = 1;
const double speedFactor = 2;
const double timeStep = 0.1;
const double moveStep = 0.1;

TEST(fillRoDistribTest, csrMatrixShouldHaveRightOrder)
{
    double **velocityX = imn<double>::matrix(size, size);
    double **velocityY = imn<double>::matrix(size, size);
    imn<double>::set_matrix(velocityX, size, size, 1.0);
    imn<double>::set_matrix(velocityY, size, size, 2.0);

    const double matrix[size * size * size * size] = {1, 0.5, 0.25, 0,
                                                      -0.5, 1, 0.5, 0.25,
                                                      -0.25, -0.5, 1, 0.5,
                                                      0, -0.25, -0.5, 1};
    CSR scatterMatrix = fillAlphaMatrix(velocityX, velocityY, size, timeStep, moveStep);
    for (int i = 0; i < size * size * size * size; ++i)
        ASSERT_EQ(matrix[i], scatterMatrix.getValue(i / (size * size), i % (size * size)));
}
}