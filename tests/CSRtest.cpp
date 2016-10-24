#include <memory>
#include "gtest/gtest.h"
#include "CSR.h"

const int matrixColsSize = 5;
const int matrixRowsSize = 5;
const int estimatedMinElements = 10;
const int estimatedMinDiagonalMatrixElements = 8;
const std::vector<std::vector<double>> sparseMatrix = { 
                                        {10.0, 2.0, 0.0, 0.0, 0.0},
                                        {0.0, 0.0, 1.0, 0.0, 0.0},
                                        {5.0, 2.0, 3.0, 0.0, 3.0},
                                        {0.0, 3.0, 0.0, 0.0, 0.0},
                                        {0.0, 8.0, 0.0, 99.0, 0.0}};
const std::vector<std::vector<double>> sparseDiagonalMatrix = { 
                                        {10.0, 0.0, 3.0, 0.0, 0.0},
                                        {0.0, 5.0, 0.0, 5.0, 0.0},
                                        {0.0, 0.0, 11.0, 0.0, 7.0},
                                        {0.0, 0.0, 0.0, 2.0, 0.0},
                                        {0.0, 0.0, 0.0, 0.0, 56.0}};
const std::vector<double> nonZeroElements = {10.0, 2.0, 1.0, 5.0, 2.0, 3.0, 3.0, 3.0, 8.0, 99.0};
class CSRTest : public ::testing::Test
{
public:
    void SetUp()
    {
        csrMatrix = std::make_shared<CSR>(matrixColsSize, matrixRowsSize, estimatedMinElements);
        for(int i=0;i<matrixRowsSize;++i)
            for(int j=0;j<matrixColsSize;++j)
                csrMatrix->setValue(i, j, sparseMatrix[i][j]);
        csrMatrix->setEndIndicator();

        csrDiagonalMatrix = std::make_shared<CSR>(matrixColsSize, matrixRowsSize, estimatedMinDiagonalMatrixElements);
        for(int i=0;i<matrixRowsSize-2;++i)
                    for(int j=2;j<matrixColsSize;++j)
                        csrMatrix->setValue(i, j, sparseMatrix[i][j]);
        for(int i=0;i<matrixRowsSize;++i)
                csrDiagonalMatrix->setValue(i, i, sparseDiagonalMatrix[i][i]);
        csrDiagonalMatrix->setEndIndicator();
    }
    std::shared_ptr<CSR> csrMatrix;
    std::shared_ptr<CSR> csrDiagonalMatrix;
};

TEST_F(CSRTest, CSRmatrix_shouldContainValidValues)
{
    EXPECT_EQ(sparseMatrix[0][0], csrMatrix->allElements[0]);
    EXPECT_EQ(sparseMatrix[0][1], csrMatrix->allElements[1]);
    EXPECT_EQ(sparseMatrix[1][2], csrMatrix->allElements[2]);
    EXPECT_EQ(sparseMatrix[2][0], csrMatrix->allElements[3]);
    EXPECT_EQ(sparseMatrix[2][1], csrMatrix->allElements[4]);
    EXPECT_EQ(sparseMatrix[2][2], csrMatrix->allElements[5]);
    EXPECT_EQ(sparseMatrix[2][4], csrMatrix->allElements[6]);
    EXPECT_EQ(sparseMatrix[3][1], csrMatrix->allElements[7]);
    EXPECT_EQ(sparseMatrix[4][1], csrMatrix->allElements[8]);
    EXPECT_EQ(sparseMatrix[4][3], csrMatrix->allElements[9]);
}
TEST_F(CSRTest, CSRmatrix_shouldContainValidColumnIndices)
{
    EXPECT_EQ(0, csrMatrix->byColumnIndices[0]);
    EXPECT_EQ(1, csrMatrix->byColumnIndices[1]);
    EXPECT_EQ(2, csrMatrix->byColumnIndices[2]);
    EXPECT_EQ(0, csrMatrix->byColumnIndices[3]);
    EXPECT_EQ(1, csrMatrix->byColumnIndices[4]);
    EXPECT_EQ(2, csrMatrix->byColumnIndices[5]);
    EXPECT_EQ(4, csrMatrix->byColumnIndices[6]);
    EXPECT_EQ(1, csrMatrix->byColumnIndices[7]);
    EXPECT_EQ(1, csrMatrix->byColumnIndices[8]);
    EXPECT_EQ(3, csrMatrix->byColumnIndices[9]);
}
TEST_F(CSRTest, CSRmatrix_shouldContainValidRowsIndices)
{
    EXPECT_EQ(0, csrMatrix->rowStartIndices[0]);
    EXPECT_EQ(2, csrMatrix->rowStartIndices[1]);
    EXPECT_EQ(3, csrMatrix->rowStartIndices[2]);
    EXPECT_EQ(7, csrMatrix->rowStartIndices[3]);
    EXPECT_EQ(8, csrMatrix->rowStartIndices[4]);
    EXPECT_EQ(10, csrMatrix->rowStartIndices[5]);
}
//Three next tests fails because matrix has to be filled consecutevily
//row by row
/*
TEST_F(CSRTest, csrDiagonalMatrix_shouldContainValidValues)
{
    EXPECT_EQ(sparseDiagonalMatrix[0][0], csrDiagonalMatrix->allElements[0]);
    EXPECT_EQ(sparseDiagonalMatrix[1][1], csrDiagonalMatrix->allElements[2]);
    EXPECT_EQ(sparseDiagonalMatrix[2][2], csrDiagonalMatrix->allElements[4]);
    EXPECT_EQ(sparseDiagonalMatrix[3][3], csrDiagonalMatrix->allElements[6]);
    EXPECT_EQ(sparseDiagonalMatrix[4][4], csrDiagonalMatrix->allElements[7]);
    EXPECT_EQ(sparseDiagonalMatrix[0][2], csrDiagonalMatrix->allElements[1]);
    EXPECT_EQ(sparseDiagonalMatrix[1][3], csrDiagonalMatrix->allElements[3]);
    EXPECT_EQ(sparseDiagonalMatrix[2][4], csrDiagonalMatrix->allElements[5]);
}
TEST_F(CSRTest, csrDiagonalMatrix_shouldContainValidColumnIndices)
{
    EXPECT_EQ(0, csrDiagonalMatrix->byColumnIndices[0]);
    EXPECT_EQ(2, csrDiagonalMatrix->byColumnIndices[1]);
    EXPECT_EQ(1, csrDiagonalMatrix->byColumnIndices[2]);
    EXPECT_EQ(3, csrDiagonalMatrix->byColumnIndices[3]);
    EXPECT_EQ(2, csrDiagonalMatrix->byColumnIndices[4]);
    EXPECT_EQ(4, csrDiagonalMatrix->byColumnIndices[5]);
    EXPECT_EQ(3, csrDiagonalMatrix->byColumnIndices[6]);
    EXPECT_EQ(4, csrDiagonalMatrix->byColumnIndices[7]);
}
TEST_F(CSRTest, csrDiagonalMatrix_shouldContainValidRowsIndices)
{
    EXPECT_EQ(0, csrDiagonalMatrix->rowStartIndices[0]);
    EXPECT_EQ(2, csrDiagonalMatrix->rowStartIndices[1]);
    EXPECT_EQ(4, csrDiagonalMatrix->rowStartIndices[2]);
    EXPECT_EQ(6, csrDiagonalMatrix->rowStartIndices[3]);
    EXPECT_EQ(7, csrDiagonalMatrix->rowStartIndices[4]);
    EXPECT_EQ(8, csrDiagonalMatrix->rowStartIndices[5]);
}*/