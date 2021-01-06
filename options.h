#ifndef OPTIONS_H
#define OPTIONS_H
#include "headers.h"
class Options
{
public:

    string summmaryFile;
    string bfileName   ;
    string outPrefix   ;
    string chrID;
    int    maxDim      ;
    int    minDim      ;
    int    maxDist;
    int    thread_num  ;
    double dupThresh   ;
    

    double pValueThreshold;
    double deltaMAF;
    string targetSNP   ;
    int64_t targetBP;
    int64_t radius;
    string extractFile;
    bool   ignoreWarnings ;
    int    withNA;
    float   propPCtrunc;

    // for extracting LD from bld file.
    string  bldLDFile;
    int   windSize;
    uint    readFrom;
    uint    readTo;
    bool     doWrite;
    bool    doCheck;
    bool    doQC;
    bool    doFreq;
    bool    doImpute;
    //const char *flgs[] ;
    vector<string> flags;
    float mafThresh     ;
    int   debugMode;
    bool  loadLD;
    bool  gcControl;
    /// to be implemented
    //  double dupThresh ;  the thresh of LD r^2 between two SNPs to be considered duplicates.
    int nIterations; // the number interations to be performed. At least one.
    int bytePerUnit;
    //  int qcMethod;    1) use zscore diff 2) zscore_diff/se
    //  int distMethod;  1) by number of snps, 2) by bp  3)by morgan distance.
    //
    //
    //
    static inline bool not_in_flags(vector<string> &flags, string str) { 
        return find(flags.begin(),flags.end(),str) == flags.end(); 
    };
    static inline bool has_prefix(const string &str, const string &prefix)
    {
        return str.size() >= prefix.size() &&
                str.compare(0, prefix.size(), prefix) == 0;
    };


    static inline void FileExist(string filename)
    {
        ifstream ifile(filename.c_str());
        if(!ifile)
        {
            cout <<  "Error: can not open the file ["+filename+"] to read." << endl;
            throw("Error: can not open the file ["+filename+"] to read.");
        }
        ifile.close();
    };

    static inline int FileExist2(string filename)
    {
        ifstream ifile(filename.c_str());
        if(!ifile)
        {
            return -1;
        }
        return 0;
    };






    static inline void bool_FLAG_VALID_CK(string str, const char* flag)
    {
        if( !(flag==NULL || has_prefix(flag, "--"))  )
        {
            // fprintf (stderr, "Please verify the flag %s!: \n",
            //          str.c_str());
            exit (EXIT_FAILURE);
        }


    }
    static inline void FLAG_VALID_CK(string str, const char* flag)
    {
        if(flag==NULL || has_prefix(flag, "--"))
        {
            fprintf (stderr, "Please verify the flag %s!: \n",
                      str.c_str());
            exit (EXIT_FAILURE);
        }
    };



    void FLAGS_VALID_CK(int option_num, char* option_str[]);
    


    inline Options()
    {
        summmaryFile  = ""; 
        flags.push_back("--gwas-summary"); //summmaryFile
        bfileName     = "";
        flags.push_back("--bfile");  // bed file

        bldLDFile = "";
        flags.push_back("--bld");  // bed file

        outPrefix     = "";  // default output prefix
        flags.push_back("--out");
        chrID         = "";  // default output prefix
        flags.push_back("--chrID");

        thread_num    = 1;      // number of threads for QC
        flags.push_back("--thread-num");
        dupThresh     = 0.99;    // percentage of probes to be filtered
        flags.push_back("--dup-threshold");
        pValueThreshold   = 5.0369e-8;    // percentage of probes to be filtered
        flags.push_back("--p-value-threshold");
        deltaMAF       =  -1;
        flags.push_back("--delta-MAF");

        maxDim        = 30000;   // default max window size for imputation
        // flags.push_back("--min-wind");
        mafThresh     = -1;   // default -1 for no restrictions on maf.
        flags.push_back("--maf");
        extractFile = "";
        flags.push_back("--extract");
        targetSNP = "";
        flags.push_back("--target");
        targetBP = -1;
        flags.push_back("--target-bp");
        radius   = 10e6;
        flags.push_back("--radius");

        withNA = 0;
        flags.push_back("--with-NA-geno");
        maxDist = 2000000;           // 2Mb by default
        flags.push_back("--wind-dist");
        minDim        = 2000;   // default min window size for imputation
        flags.push_back("--wind");
        debugMode  = 0;
        flags.push_back("--debug");
        loadLD = false;
        flags.push_back("--load-LD");

        nIterations = 8;
        flags.push_back("--iteration-num");
        bytePerUnit = 2;

        flags.push_back("--LD-unit-in-byte");

        gcControl = false;
        flags.push_back("--control-inflation");
        propPCtrunc = 0.5;
        flags.push_back("--SVD-trunc-prop");

   

        ignoreWarnings = false;

        readFrom = 0;
        readTo   = 0;
        doWrite       = false;
        doCheck       = false;
        doQC          = true;
        doFreq        = false;
        doImpute      = false;
        flags.push_back("--check-LD");
        flags.push_back("--write-LD");
        flags.push_back("--freq");
        flags.push_back("--impute");

        ignoreWarnings = false;
        flags.push_back("--ignore-warnings");


    }


    void parseOptions(int nArgs, char* option_str[]);



};


void Options::FLAGS_VALID_CK(int option_num, char* option_str[])
{
    if(option_num<3)
    {
        cout<<"Flags include:"<<endl;
        int cur_mark=0;
        for(int i=0;i<flags.size();i++)
        {
            int tmp=i>>2;
            if(tmp>cur_mark)
            {
                cout<<endl;
                cur_mark=tmp;
            }
            cout<<flags[i]<<",";
        }
        cout<<endl;
        exit (EXIT_FAILURE);
    }
    for(int i=0;i<option_num;i++)
    {
        if(has_prefix(option_str[i],"--"))
            if(not_in_flags(flags, option_str[i]))
            {
                //fprintf (stderr, "%s: Invalid option\n", option_str[i]);
                exit (EXIT_FAILURE);
            }
    }

}

void Options::parseOptions(int nArgs, char* option_str[])
{
    FLAGS_VALID_CK(nArgs, option_str);
    for(int i = 0; i < nArgs; i ++)
    {
        if(strcmp(option_str[i], "--load-LD") == 0)
        {
            loadLD = true;
            if(i+1 < nArgs)
                Options::bool_FLAG_VALID_CK(string("--load-LD"), option_str[i+1]);
            cout<< option_str[i] << " " <<  " TRUE" <<endl;
        }
        if(strcmp(option_str[i], "--control-inflation") == 0)
        {
            this->gcControl  = true;
            if(i+1 < nArgs)
                Options::bool_FLAG_VALID_CK(string("--control-inflation"), option_str[i+1]);
            cout<< option_str[i] << " " <<  " TRUE" <<endl;
        }
        if(strcmp(option_str[i], "--iteration-num") == 0){
            nIterations = atoi(option_str[++i]);
            if(nIterations < 1 || nIterations > 10) 
            {
                fprintf (stderr, "Error: --iteration-num should be between 1 and 10,  [1 10].\n");
                exit (EXIT_FAILURE);
            }
            cout << option_str[i-1] << " " << nIterations << endl;
        }

        if(strcmp(option_str[i], "--LD-unit-in-byte") == 0){
            bytePerUnit = atoi(option_str[++i]);
            if(!(bytePerUnit  ==1 || bytePerUnit  == 2 || bytePerUnit  == 4 || bytePerUnit  == 8) )
            {
                fprintf (stderr, "Error: --LD-unit-in-byte should be in [1,2,4,8] for char, short, int, double \n");
                exit (EXIT_FAILURE);
            }
            cout << option_str[i-1] << " " <<  bytePerUnit << endl;
        }







        if(strcmp(option_str[i], "--bld") == 0)
        {
            bldLDFile = option_str[++i];
            Options::FLAG_VALID_CK(string("--bld"), bldLDFile.c_str());
            cout<< option_str[i-1] << " " <<  bldLDFile <<endl;
            FileExist(bldLDFile + ".bld");
        }

        if(strcmp(option_str[i],"--check-LD")==0)
        {
            this->doCheck = true;
            this->doQC  = false;
            readFrom  = atoi(option_str[++i]);
            readTo    = atoi(option_str[++i]);
            //cout << option_str[i-2] << "\t" << " from ith SNP " << readFrom << " to jth SNP "
            //        << readTo << endl;
            printf("%s\t from %d-th SNP to %d-th SNP.\n", option_str[i-2], readFrom, readTo);
            if(readTo - readFrom < 1) 
                stop( "Error: the dimention of matrix < 1 \n");
            FileExist(bldLDFile + ".bld");
        }
        if(strcmp(option_str[i],"--freq")==0)
        {
            this->doFreq  = true;
            this->doQC    = false;
            this->doCheck = false;
            if(i+1 < nArgs)
                Options::bool_FLAG_VALID_CK(string("--freq"), option_str[i+1]);
            cout<< option_str[i] << " " <<  " TRUE" <<endl;
        }

        if(strcmp(option_str[i],"--write-LD")==0)
        {
            this->doWrite  = true;
            this->doQC  = false;
            if(i+1 < nArgs)
                Options::bool_FLAG_VALID_CK(string("--with-LD"), option_str[i+1]);
            cout<< option_str[i] << " " <<  " TRUE" <<endl;
        }

        if(strcmp(option_str[i],"--impute")==0)
        {
            this->doQC  = false;
            this->doImpute  = true;
            if(i+1 < nArgs)
                Options::bool_FLAG_VALID_CK(string("--impute"), option_str[i+1]);
            cout<< option_str[i] << " " <<  " TRUE" <<endl;
        }



        if(strcmp(option_str[i],"--gwas-summary")==0)
        {
            summmaryFile = option_str[++i];
            Options::FLAG_VALID_CK(string("--gwas-summmary"), summmaryFile.c_str());
            cout<< option_str[i-1] << " " <<  summmaryFile <<endl;
            FileExist(summmaryFile);
        }


        if(strcmp(option_str[i], "--bfile") == 0)
        {
            bfileName = option_str[++i];
            Options::FLAG_VALID_CK(string("--bfile"), bfileName.c_str());
            cout<< option_str[i-1] << " " <<  bfileName <<endl;
            FileExist(bfileName + ".bed");
            FileExist(bfileName + ".bim");
            FileExist(bfileName + ".fam");
        }
        if(strcmp(option_str[i], "--out") == 0)
        {
            outPrefix = option_str[++i];
            Options::FLAG_VALID_CK(string("--out"), outPrefix.c_str());
            cout<< option_str[i-1] << " "<< outPrefix <<endl;
        }
        if(strcmp(option_str[i], "--chrID") == 0)
        {
            chrID  = option_str[++i];
            Options::FLAG_VALID_CK(string("--chrID"), chrID.c_str());
            cout<< option_str[i-1] << " "<<  chrID <<endl;
        }

        if(strcmp(option_str[i], "--thread-num") == 0){
            thread_num = atoi(option_str[++i]);
            cout << option_str[i-1] << " " << thread_num<< endl;
            if(thread_num < 0)
            {
                fprintf (stderr, "Error: --thread-num should be over 0.\n");
                exit (EXIT_FAILURE);
            }
        }

        if(strcmp(option_str[i], "--delta-MAF") == 0){
            deltaMAF = atof(option_str[++i]);
            cout << option_str[i-1] << " " << deltaMAF << endl;
            if(deltaMAF < 0 || deltaMAF > 1) 
            {
                fprintf (stderr, "Error: --delta-MAF should be between 0 and 1 that is [0 1].\n");
                exit (EXIT_FAILURE);
            }
        }



        if(strcmp(option_str[i], "--dup-threshold") == 0){
            dupThresh = atof(option_str[++i]);
            cout << option_str[i-1] << " " << dupThresh << endl;
            if(dupThresh  < 0.01 || dupThresh > 1) 
            {
                fprintf (stderr, "Error: --dup-threshold should be between 0.01 and 1 that is [0.01 1].\n");
                exit (EXIT_FAILURE);
            }
        }
        if(strcmp(option_str[i], "--maf") == 0)
        {
            mafThresh = atof ( option_str[++i] );
            Options::FLAG_VALID_CK(string("--maf"), option_str[i] );
            cout<< option_str[i-1] << " " <<  mafThresh <<endl;
        }


        if(strcmp(option_str[i], "--wind-dist") == 0)
        {
            maxDist = atof ( option_str[++i] );
            Options::FLAG_VALID_CK(string("--wind-dist"), option_str[i] );
            cout<< option_str[i-1] << " " <<  maxDist <<endl;
        }

        if(strcmp(option_str[i], "--wind") == 0)
        {
            maxDim = atof ( option_str[++i] );
            Options::FLAG_VALID_CK(string("--wind"), option_str[i] );
            cout<< option_str[i-1] << " " <<  maxDim <<endl;
        }
        
        if(strcmp(option_str[i], "--target") == 0)
        {
            targetSNP = ( option_str[++i] );
            Options::FLAG_VALID_CK(string("--target"), targetSNP.c_str());
            cout<< option_str[i-1] << " " <<  targetSNP <<endl;
        }
        if(strcmp(option_str[i], "--target-bp") == 0)
        {
            std::istringstream ss(option_str[++i]);
            if (!(ss >> targetBP))
                stop("Fail to convert --target <value> into int64\n");
            Options::FLAG_VALID_CK(string("--target-bp"), option_str[i] );
            cout<< option_str[i-1] << " " <<  targetBP <<endl;
        }
        if(strcmp(option_str[i], "--radius") == 0)
        {
            std::istringstream ss(option_str[++i]);
            if (!(ss >> radius))
                stop("Fail to convert --target <value> into int64\n");
            Options::FLAG_VALID_CK(string("--radius"), option_str[i] );
            cout<< option_str[i-1] << " " <<  radius <<endl;
        }



        if(strcmp(option_str[i], "--extract") == 0)
        {
            extractFile = ( option_str[++i] );
            Options::FLAG_VALID_CK(string("--extract"), extractFile.c_str());
            cout<< option_str[i-1] << " " <<  extractFile <<endl;
        }

        if(strcmp(option_str[i], "--ignore-warnings") == 0)
        {
            ignoreWarnings = true;
            if(i+1 < nArgs)
            {
                Options::bool_FLAG_VALID_CK(string("--ignore-warnings"), option_str[i+1]);
            }
            cout<< option_str[i] << " " <<  " TRUE" <<endl;
        }
        if(strcmp(option_str[i], "--with-NA-geno") == 0)
        {
            withNA  = 1;
            if(i+1 < nArgs)
            {
                Options::bool_FLAG_VALID_CK(string("--with-NA-geno"), option_str[i+1]);
            }
            cout<< option_str[i] << " " <<  " TRUE" <<endl;
        }

        if(strcmp(option_str[i], "--debug") == 0)
        {
            debugMode = 1;
            if(i+1 < nArgs)
            {
                Options::bool_FLAG_VALID_CK(string("--debug"), option_str[i+1]);
            }
            cout<< option_str[i] << " " <<  " TRUE" <<endl;
        }


        if(strcmp(option_str[i], "--SVD-trunc-prop") == 0)
        {
            propPCtrunc = atof(option_str[++i]);
            Options::FLAG_VALID_CK(string("--SVD-trunc-prop"), option_str[i] );
            cout<< option_str[i-1] << " " <<  propPCtrunc <<endl;
            if(propPCtrunc <= 0 || propPCtrunc > 1)
                stop("[errors] --SVD-trunc-prop expects a number between (0,1].\n");
        }

    }
    if(outPrefix == "")
    {
        outPrefix = "out";
        cout<< "--out" << " " << outPrefix <<endl;

    }

}

#endif// OPTIONS_H

