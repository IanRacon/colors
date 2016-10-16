#pragma once
#include <vector>

class CSR
{
public:
    CSR(int absoluteMatrixCols, int absoluteMatrixRows, int estimatedMaxElements=0);
    void setValue(int col, int row, double value);
private:
    int absoluteMatrixCols;
    int absoluteMatrixRows;

    std::vector<double> allElements;
    std::vector<int> byColumnIndices;
    std::vector<int> rowStartIndices;

    void markZeroFilledRows(int presentRow);
    void setValueIndexInRow(int presentRow);
};