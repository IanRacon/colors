#include "procedures.h"
#include <iostream>
#include <numeric>

namespace procedures
{
    std::vector<double> product(const CSR &csrMatrix, const std::vector<double>& vector){
        std::vector<double> result;
        result.reserve(vector.size());
        int indexOfFirstElementInRow;
        int indexOfLastElementInRow;
        if(vector.size() != csrMatrix.rows())
            std::cout << "In procedure product: Size of given vector not the same as csrMatrix!\n";
        for(int i=0; i<=csrMatrix.rows();++i){
            result[i]=0.0;
            indexOfFirstElementInRow = csrMatrix.rowStartIndices[i];
            indexOfLastElementInRow = csrMatrix.rowStartIndices[i+1]-1;
            for(int j=indexOfFirstElementInRow;j<indexOfLastElementInRow;++j)
                result[i] += csrMatrix.allElements[j]*vector[csrMatrix.byColumnIndices[j]];        
        }
        return result;
    }
}