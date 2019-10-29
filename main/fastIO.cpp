#include "dataTypeDef.h"



// ****************************************************************************
/// This is for faster io (than read.table in R) of bim, frq, summary Stat.
//  
//

void readBim(char** prefix, char** chrID, char** rsID,  int* M)
{
    std::string prefixStr = std::string(*prefix);


    FILE* bimReader = NULL;
    if ((bimReader = fopen ( (prefixStr  + ".bim").c_str(), "r"))==NULL) {
        printf("\t[error] Cannot open file %s\n", (prefixStr  + ".bim").c_str());
    }

//    for(uint i = 0; i < M; i ++)
//        fscanf(bimReader,"%s %s", chrID,rsID);
    fclose(bimReader);
}


