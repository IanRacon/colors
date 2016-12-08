#pragma once
#include <vector>

class CSR
{
  public:
    CSR(int maxRows, int maxCols, int estimatedMinElements = 0);
    void setValue(int row, int col, double value);
    double getValue(int row, int col) const;
    void setEndIndicator();
    int rows() const;
    int cols() const;
    std::vector<double> allElements;
    std::vector<int> byColumnIndices;
    std::vector<int> rowStartIndices;
    int absoluteMatrixRows;

  private:
    int absoluteMatrixCols;
    void markZeroFilledRows(int presentRow);
    void setValueIndexInRow(int presentRow);
};