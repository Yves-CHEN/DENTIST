#ifndef SORTEDHEAP
#define SORTEDHEAP
#include "dataTypeDef.h"

template <class TT>
class sortedHeap
{
public:
    int  vecSize;
    TT* vec ;
    int stackSize;
    std::list<int> mylist;
    inline void copyToArray(int* arr) { copy(mylist.begin(),mylist.end(),arr); };
    inline int size() {return mylist.size();};
    sortedHeap (int vecSize, int stackSize, TT* vec)
    :vecSize(vecSize),stackSize(stackSize),vec(vec) {
        for(uint i = 0; i < vecSize; i ++)
            insert(i);
    };
    void  insert(uint idx)
    {
        if(idx > vecSize){ std::cout << "error occured! The index exceeds the size of vec" << std::endl; }
        TT valueIndexed = vec[idx];
        std::list<int>::iterator  it;
        for (it=mylist.begin(); it!=mylist.end(); ++it)
        {
            if(valueIndexed > vec [*it]) break;
        }
        mylist.insert (it, idx);
        if(mylist.size() > stackSize) /// delete the tail if exceeding
        {
            //cout << "delete tail element"  << vec[*(--mylist.end())] << endl;
            mylist.pop_back();
        }
    };
    void  print()
    {
        std::cout << "size: " << mylist.size() << std::endl;
        std::list<int>::iterator  it;
        for (it=mylist.begin(); it!=mylist.end(); ++it)
        {
            std::cout << *it << " ";
        }
        std::cout << "----------- "  << std::endl;
    };
    bool findDelete(int idx)
    {
        std::list<int>::iterator  it;
        for (it=mylist.begin(); it!=mylist.end(); ++it)
        {
            if( idx == *it) break;
        }
        if(it!=mylist.end())
        {
            mylist.erase(it);
            return 0;
        }
        else
        {
        return 1;
        }
    }
};
#endif
