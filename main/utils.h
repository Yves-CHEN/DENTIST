
#ifndef UTILS
#define UTILS

//#include "dataTypeDef.h"
#include <numeric>
#include <vector>
#include <algorithm>
//



extern "C" {

void findDup(int* r, int* len, int* ifDup)
{
    int dim = *len;
    for (int i =0; i < dim; i ++)
    {
        if(ifDup [i] != 0 ) continue;
        for (int j =i+1; j < dim ; j ++)
            ifDup[j] += r[i*dim + j];
    }


}
};





template <typename T>
std::vector<size_t> sort_indexes(const std::vector<T>& v, unsigned int theSize) 
{
// initialize original index locations
    std::vector<size_t> idx(theSize);
    iota(idx.begin(), idx.end(), 0);
    sort(idx.begin(), idx.end(),   [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});
    return idx;
}



std::vector<size_t>   generateSetOfNumbers(int SIZE, int seed)
{
    int tempNum;            // temp variable to hold random number
    std::vector <int> numbers ( SIZE, 0);
    srand(seed);  /// fixed seeding instead of srand(time(NULL));
    numbers[0] = rand() ;   // generate the first number in the array
    for (int index = 1; index < SIZE; index++)  // loop to place other numbers
    {
        do
        {
            tempNum = rand() ;
            for (int index2 = 0; index2 < SIZE; index2++)
                if (tempNum == numbers[index2]) tempNum = -1;
        } while (tempNum == -1);
        numbers[index] = tempNum;
    }
    std::vector<size_t>  aa = sort_indexes<int>(numbers, SIZE);
    return aa;

}
#endif

