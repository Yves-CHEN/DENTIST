# Overview
DENTIST (Detecting Errors iN analyses of summary staTISTics) is a QC tool for summary-data-based analyses, given that such analyses based on summary statistics from association studies are susceptible to errors in summary statistics either originated from genotyping/imputation errors or other artifacts and  data heterogeneity as out-of-sample LD is usually used to approximate the in-sample LD. Our paper has summarized the effects of errors and data heterogeneity on commonly used analyses ( [GCTA COJO](https://cnsgenomics.com/software/gcta/#COJO),   [SMR HEIDI test](https://cnsgenomics.com/software/smr/#SMR&HEIDIanalysis),     [LD score regression]( https://github.com/bulik/ldsc) ), and how much improvement can be attained using DENTIST-based QC.

# Credits
The method is designed by Wenhan Chen, Zhihong Zhu and Jian Yang. The software is implemented and maintained by Wenhan Chen.  The idea of the LD consistency test is based on a previous study [PubMed ID: 24990607]. 

# Questions and Help Requests
If you meet any bugs or questions, please send an email to [Wenhan Chen](mailto:uqwche11@uq.edu.au) or [Jian Yang](mailto:jian.yang@uq.edu.au).

# Citations
In progress.

# Downloads
### Pre-compiled Executable Files 
The executable file below is compiled with "-static", and tested on 64-bit Linux distributions on the x86\_64 CPU platform.

Linux: [DENTIST  0.8.0.1](https://drive.google.com/open?id=1BYC3SDim9NBfjyYmWsXRGJcZVCUQh6_o)

# Basic usage
Note that the words starting with a $ sign are bash variables. The bash syntax can be found [here](https://linuxconfig.org/bash-scripting-tutorial-for-beginners).

To Run DENTIST with essential parameters,
> DENTIST --gwas-summary $summary --bfile $ref --out $filePrefix

To specify the number of cpus,
> DENTIST --gwas-summary $summary --bfile $ref --out $filePrefix   --thread-num $ncpus

To run at a target region specified by the a rsID (e.g. rs1234),
>DENTIST --gwas-summary $summary --bfile $ref --out $filePrefix --thread-num $ncpus --target $rsID


# Input and output

> \-\-bfile \<PLINK bed file prefix\>

Reads individual-level genotype data in PLINK bed format, e.g. test.fam, test.bim and test.bed.

> \-\-gwas-summary \<summary statistics in GCTA-COJO format\>

Reads GWAS summary data in in GCTA-COJO format.  Format of the GCTA-COJO file is like this,
```
SNP A1 A2 freq beta se p N
rs131538 A G 0.05 0.007 0.02 0.7 6000
rs140378 C G 0.05 0.007 0.02 0.7 6000
...
```
> \-\-out \<the output file prefix\>

Specifies the prefix of output file. The output of the QC test statistics is in the following format,
```
rsID        chisq   -log10(P)  ifDup
rs131538    0.012   0.91       0
rs140378    14.4    2          0
...
```
For each variants,  the first column specifies the rsID, followed by the QC test statistics used in DENTIST (which is a $\chi^2$ test), the p-value of the test statistics in -log10 (*) and a indicator of whether the variant is in strong correlation (|r| > 0.99) with any other variant, 0 for none is found and 1 for at least one is found.

# Parameters
> \-\-target \<rsID\>

Is trailed by an rsID to specify a region of 20Mb of interest centered at position specified by rsID. The rsID should present in the Plink bed file. A warning is reported if the target rsID is not found and DENTIST will run across the chromosome rather than the specified regions by --target. Notably, the locating of this rsID in the bfile is performed before --maf flag.

> \-\-maf \<minimal MAF threshold\>

Specifies the MAF threshold, so that markers from the bfile with a MAF smaller than this value would be precluded from the analysis. This exclusion is performed after --target.

> \-\-wind-dist \<window size\>

Is trailed by sliding window size measured in the number of base pairs. The default value is 2000000 base pairs. This is the default option unless overruled by --wind.

> \-\-thread-num \<ncpus\>

Specifies the number of threads for parallel computing, given that DENTIST is powered by OpenMP. The default value is 1.

> \-\-num-iterations \<number of iterations\>

Specifies the number of iterations for z-score consistency test (see Method). A large value increases the cost of computation, but a smaller value can reduce the accuracy and increase the false discovery rate. We experimented with this parameter and set it to a default value of 10, which presents the best trade-off.


> \-\-delta-MAF \<threshold\>

Specifies the threshold for the acceptable differences in MAFs between a reference sample and a summary data. This filter is not applied by default. The usually used thresholds include 0.1 and 0.2.
