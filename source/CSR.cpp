#include "CSR.h"
#include "easylogging++.h"

CSR::CSR(int absoluteMatrixCols, int absoluteMatrixRows, int estimatedMaxElements):
absoluteMatrixCols(absoluteMatrixCols), 
absoluteMatrixRows(absoluteMatrixRows)
{
    allElements.reserve(estimatedMaxElements);
    byColumnIndices.reserve(estimatedMaxElements);
    rowStartIndices.reserve(absoluteMatrixRows);
}
void CSR::setValue(int col, int row, double value)
{
    if(value != 0.0)
    {
        allElements.push_back(value); 
        byColumnIndices.push_back(col);
        setValueIndexInRow(row);
    }
}
void CSR::setValueIndexInRow(int presentRow)
{
    if(presentRow > rowStartIndices.size()+1){
        markZeroFilledRows(presentRow-rowStartIndices.size()+1);
        LOG(ERROR) << "All values zero in rows: " << rowStartIndices.size()+1 << ", to: " << presentRow;
    }
    else if(presentRow==rowStartIndices.size()+1)
        rowStartIndices.push_back(allElements.size()); 
}
void CSR::markZeroFilledRows(int numberOfZeroRows)
{
    for(int i=0;i<numberOfZeroRows;++i)
        rowStartIndices.push_back(0);
}
