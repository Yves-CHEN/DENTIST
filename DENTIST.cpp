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


    // system/platform check
    assert (sizeof(off_t) == 8);
    assert (sizeof(int64) == 8);



    string summmaryFile  = "/home/uqwche11/30days/simulation.UK10K//anaHeight/summaryData/height_UKB_FULLSAMPLE_Julia.txt";
    string bfileName     = "/home/uqwche11/30days/simulation.UK10K/Data/hrs.hm2/hrs_hm2_chr22";
    Options opt;
    opt.parseOptions (argc, argv);
 //   if(opt.chrID == "")
 //       stop("[error] Please specify chrID by --chrID, since DENTIST only run for each chromosome.\n");


    if(opt.bldLDFile != "" && opt.doCheck)
    {
        bool ifPrint = true;
        int dim = opt.readTo - opt.readFrom ;
        readLDFromFile_FromTo(opt.bldLDFile, dim , opt.readFrom, opt.readTo, ifPrint);
    }

    if( opt.doWrite)
    {
//        saveLD<int>(opt.bfileName, opt.outPrefix.c_str(), 
//                            opt.maxDist, opt.thread_num);
        saveLD<int>(opt);
    }

    if( opt.doFreq)
    {
        getFreq(opt);
    }

    if(opt.doImpute && opt.summmaryFile != "" && opt.bfileName != "")
        runSummaryImpute(opt);

    if(opt.doQC && opt.summmaryFile != "" )
        runQC(opt);

    return 0;
}

