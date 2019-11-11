#ifndef OPTIONS_H
#define OPTIONS_H
#include "headers.h"
class Options
{
public:

    string summmaryFile;
    string bfileName   ;
    string outPrefix   ;
    int    maxDim      ;
    int    minDim      ;
    int    maxDist;
    int    thread_num  ;
    double Degree_QC   ;
    double lambda      ;
    string targetSNP   ;
    string extractFile;
    bool   ignoreWarnings ;
    int    withNA;

    // for extracting LD from bld file.
    string  bldLDFile;
    uint    readFrom;
    uint    readTo;
    //const char *flgs[] ;
    vector<string> flags;
    float mafThresh     ;
    int   debugMode;
    /// to be implemented
    //  double dupThresh ;  the thresh of LD r^2 between two SNPs to be considered duplicates.
    //  int nIterations; the number interations to be performed. At least one.
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
        outPrefix     = "out";  // default output prefix
        flags.push_back("--out");
        thread_num    = 1;      // number of threads for QC
        flags.push_back("--thread-num");
        Degree_QC     = 0.3;    // percentage of probes to be filtered
        // flags.push_back("--degree-of-QC");
        lambda        = 0.1;    // the lambda for Ridge regression 
        // flags.push_back("--lambda");
        maxDim        = 25000;   // default max window size for imputation
        // flags.push_back("--min-wind");
        mafThresh     = -1;   // default -1 for no restrictions on maf.
        flags.push_back("--maf");
        extractFile = "";
        flags.push_back("--extract");
        targetSNP = "";
        flags.push_back("--target");
        withNA = 0;
        flags.push_back("--with-NA-geno");
        maxDist = -1;
        flags.push_back("--wind-dist");
        minDim        = 2000;   // default min window size for imputation
        flags.push_back("--wind");
        debugMode  = 0;
        flags.push_back("--debug");
        ignoreWarnings = false;
        flags.push_back("--ignore-warnings");

        ignoreWarnings = false;

        bldLDFile ="";
        readFrom = 0;
        readTo   = 0;
        flags.push_back("--read-LD");



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
        if(strcmp(option_str[i],"--read-LD")==0)
        {
            bldLDFile = option_str[++i];
            Options::FLAG_VALID_CK(string("--read-LD"), bldLDFile.c_str());
            cout<< option_str[i-1] << " " <<  bldLDFile <<endl;
            readFrom  = atoi(option_str[++i]);
            readTo    = atoi(option_str[++i]);
            cout << "\t" << " from ith SNP " << readFrom << " to jth SNP "
                    << readTo << endl;
            if(readTo - readFrom < 1) 
            {
                fprintf (stderr, "Error: the dimention of matrix < 1 \n");
                exit (EXIT_FAILURE);
            }
            FileExist(bldLDFile + ".bld");
            FileExist(bldLDFile + ".ridx");

            return ;
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
            // CommFunc::FileExist(oproblstName);
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

        if(strcmp(option_str[i], "--degree-of-QC") == 0){
            Degree_QC = atof(option_str[++i]);
            cout << option_str[i-1] << " " << Degree_QC  << endl;
            if(Degree_QC  < 0 || Degree_QC >= 1) 
            {
                fprintf (stderr, "Error: --degree-of-QC should be between 0 and 1.\n");
                exit (EXIT_FAILURE);
            }
        }
        if(strcmp(option_str[i], "--lambda") == 0){
            lambda = atof(option_str[++i]);
            cout << option_str[i-1] << " " << lambda << endl;
            if(Degree_QC  < 0 || Degree_QC >= 0.1) 
            {
                fprintf (stderr, "Error: --lambda should be between 0 and 0.1.\n");
                exit (EXIT_FAILURE);
            }
        }
        if(strcmp(option_str[i], "--maf") == 0)
        {
            mafThresh = atof ( option_str[++i] );
            Options::FLAG_VALID_CK(string("--maf"), bfileName.c_str());
            cout<< option_str[i-1] << " " <<  mafThresh <<endl;
        }


        if(strcmp(option_str[i], "--wind-dist") == 0)
        {
            maxDist = atof ( option_str[++i] );
            Options::FLAG_VALID_CK(string("--wind-dist"), bfileName.c_str());
            cout<< option_str[i-1] << " " <<  maxDist <<endl;
        }

        if(strcmp(option_str[i], "--wind") == 0)
        {
            maxDim = atof ( option_str[++i] );
            Options::FLAG_VALID_CK(string("--wind"), bfileName.c_str());
            cout<< option_str[i-1] << " " <<  maxDim <<endl;
        }

        if(strcmp(option_str[i], "--target") == 0)
        {
            targetSNP = ( option_str[++i] );
            Options::FLAG_VALID_CK(string("--target"), targetSNP.c_str());
            cout<< option_str[i-1] << " " <<  targetSNP <<endl;
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



            

    }

}

#endif// OPTIONS_H

