#include "DENTIST.h"
#include "version.h"
//#include <iostream>
//#include <algorithm>
//#include <numeric> 
//#include <vector> 
//#include <string> 
//
//





int main(int argc, char* argv[]) 
{


cout << "*******************************************************************" << endl;
cout << "* DENTIST (Detecting Errors iN analyses of summary staTISTics)" << endl;
cout << "* Version "<< AutoVersion::FULLVERSION_STRING << endl;
cout << "* (C) 2018 Wenhan Chen, Zhihong Zhu and Jian Yang" << endl;
cout << "* The University of Queensland" << endl;
cout << "* MIT License" << endl;
cout << "*******************************************************************" << endl;





    string summmaryFile  = "/home/uqwche11/30days/simulation.UK10K//anaHeight/summaryData/height_UKB_FULLSAMPLE_Julia.txt";
    string bfileName     = "/home/uqwche11/30days/simulation.UK10K/Data/hrs.hm2/hrs_hm2_chr22";
    Options opt;
    opt.parseOptions (argc, argv);



    runQC(opt);
    

    return 0;
}

