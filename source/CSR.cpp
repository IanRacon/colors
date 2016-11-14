#include "CSR.h"
#include <iostream>
//#include "easylogging++.h"

CSR::CSR(int absoluteMatrixCols, int absoluteMatrixRows, int estimatedMinElements) : absoluteMatrixCols(absoluteMatrixCols),
                                                                                     absoluteMatrixRows(absoluteMatrixRows)
{
    allElements.reserve(estimatedMinElements + 1);
    byColumnIndices.reserve(estimatedMinElements + 1);
    rowStartIndices.reserve(absoluteMatrixRows + 1);
}
void CSR::setValue(int row, int col, double value)
{
    if (value != 0.0)
    {
        allElements.push_back(value);
        byColumnIndices.push_back(col);
        setValueIndexInRow(row);
    }
}
void CSR::setValueIndexInRow(int presentRow)
{
    if (presentRow > rowStartIndices.size())
    {
        markZeroFilledRows(presentRow - rowStartIndices.size() - 1);
        rowStartIndices.push_back(allElements.size() - 1);
    }
    // else if(presentRow==rowStartIndices.size()+1)
    if (presentRow == rowStartIndices.size())
        rowStartIndices.push_back(allElements.size() - 1);
}
void CSR::markZeroFilledRows(int numberOfZeroRows)
{
    for (int i = 0; i < numberOfZeroRows; ++i)
        rowStartIndices.push_back(allElements.size() - 1);
}
void CSR::setEndIndicator()
{
    rowStartIndices.push_back(allElements.size());
}
int CSR::rows() const
{
    return absoluteMatrixRows;
}
int CSR::cols() const
{
    return absoluteMatrixCols;
}
double CSR::getValue(int row, int col) const
{
    int counter = 0;
    for (int j = rowStartIndices[row]; j < rowStartIndices[row + 1]; ++j)
    {
        if (byColumnIndices[j] == col)
            return allElements[rowStartIndices[row] + counter];
        counter++;
    }
    return 0;
}
