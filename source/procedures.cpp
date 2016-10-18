#include "procedures.h"
#include <iostream>
#include <numeric>

namespace procedures
{
    std::vector<double> product(const CSR &csrMatrix, const std::vector<double>& myVector){
        if(myVector.size() != csrMatrix.rows())
            std::cout << "In procedure product: Size of given vector not the same as csrMatrix!\n";
        std::vector<double> result;
        result.reserve(myVector.size());
        int indexOfFirstElementInRow;
        int indexOfLastElementInRow;
        for(int i=0; i<csrMatrix.rows();++i){
            result.push_back(0.0);
            indexOfFirstElementInRow = csrMatrix.rowStartIndices[i];
            indexOfLastElementInRow = csrMatrix.rowStartIndices[i+1]-1;
            for(int j=indexOfFirstElementInRow;j<indexOfLastElementInRow+1;++j)
                result.at(i) += csrMatrix.allElements.at(j)*myVector.at(csrMatrix.byColumnIndices.at(j));    
        }
        return result;
    }
}