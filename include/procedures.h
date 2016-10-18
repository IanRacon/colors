#pragma once
#include "CSR.h"

namespace procedures
{
    std::vector<double> product(const CSR &csrMatrix, const std::vector<double>& vector);
}