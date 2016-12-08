#include "utils.h"
//#include "easylogging++.h"

namespace utils
{
void fill2Darray(int x1, int y1, int x2, int y2, double **array, double val)
{
    for(int i=x1;i<x2;++i)
        for(int j=y1;j<y2;++j)
            array[i][j] = val;
}
void logArray(int x1, int y1, int x2, int y2, double **array){
	for(int i=x1;i<x2;++i)
		for(int j=y1;j<y2;++j)
            return;
			//LOG(INFO) << array[i][j]; 
}
}