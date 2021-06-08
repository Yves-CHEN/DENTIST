
# Overview
DENTIST (Detecting Errors iN analyses of summary staTISTics) is a quality control (QC) tool for summary-level data from genome-wide association studies (GWASs). It leverages the difference between the observed GWAS test-statistic of a variant and its predicted value (using the neighbouring variants and linkage equilibrium (LD) data from a reference panel) to remove problematic variants. It can detect genotyping/imputation errors in either the original GWAS or the LD reference samples, allelic errors (i.e., the effect alleles of the variants are mislabelled) in the GWAS summary data, as well as heterogeneity between the GWAS and LD reference samples. As shown in [our paper](#Citations), DENTIST can significantly improve the performance of the summary-data based conditional and jointly association analysis (COJO; Yang et al. 2012 Nat Genet) especially for rare variants. It can also improve the performance of LD score regression (LDSC; Builk-Sulivant et al. 2015 Nat Genet) and SMR-HEIDI (Zhu et al. 2015 Nat Genet). We believe that DENTIST can in principle work for all analyses that use GWAS summary data, although this has not been tested extensively. Apart from its application to summary data based analyses, DENTIST can also  be used in complementary to the conventional QCs in GWAS with individual-level data.

# Credits
The method is developed by Wenhan Chen, Zhihong Zhu and [Jian Yang](https://publons.com/researcher/2848531/jian-yang/) in the Yang Lab at The University of Queensland. The software is programmed and maintained by Wenhan Chen.

# Questions and Help Requests
If you meet any bugs or questions, please send an email to [Wenhan Chen](mailto:uqwche11@uq.edu.au) or [Jian Yang](mailto:jian.yang.qt@gmail.com).

# Citations
Chen, W., Wu, Y., Zheng, Z., Qi, T., Visscher, P. M., Zhu, Z., & Yang, J. (2020). Improved analyses of GWAS summary statistics by reducing data heterogeneity and errors. bioRxiv.


# Install
For compiling yourself, I have not made a user-friendly MAKEFILE yet. For now, you can type
```
git clone https://github.com/Yves-CHEN/DENTIST
cd DENTIST
make
```
and address the dependencies yourself on MKL, BOOST and EIGEN Libraries.


# Downloads
### Pre-compiled Executable Files
The executable file below is compiled with "-static" and tested on 64-bit Linux distributions on the x86_64 CPU platform,
[DENTIST_1.1.0.0](https://www.dropbox.com/s/1mtskir8qzqsmee/DENTIST.1.1.0.0.gz?dl=0)
[DENTIST  0.9.2.1](https://www.dropbox.com/s/37bc35azxwbzdos/DENTIST.0.9.2.1.gz?dl=0).

To download, you can run
> wget -O DENTIST_1.1.0.0.gz https://www.dropbox.com/s/1mtskir8qzqsmee/DENTIST.1.1.0.0.gz?dl=0


# Basic usage
To run DENTIST with essential parameters,
> DENTIST --gwas-summary summary_data --bfile ref --out prefix

To specify the number of CPUs,
> DENTIST --gwas-summary summary_data --bfile ref --out prefix --thread-num 10

To run DENTIST at a targeted region specified by the rsID of variant (e.g. rs1234),
>DENTIST --gwas-summary summary_data --bfile ref --out prefix --thread-num 10 --target rs1234


# Input and output

> \-\-bfile \<PLINK bed file prefix\>

Reads individual-level genotype data in PLINK bed format, e.g. test.fam, test.bim and test.bed.

> \-\-gwas-summary \<summary statistics in GCTA-COJO format\>

Reads GWAS summary data in in GCTA-COJO format (see below for an example).
```
SNP A1 A2 freq beta se p N
rs131538 A G 0.05 0.007 0.02 0.7 6000
rs140378 C G 0.05 0.007 0.02 0.7 6000
...
```
> \-\-out \<the output file prefix\>

Specifies the prefix of the output file. The output of DENTIST is in the following format.
```
rsID        chisq   -log10(P)  ifDup
rs131538    0.012   0.91       0
rs140378    14.4    2          0
...
```
For each variant, the first column shows the rsID, followed by the DENTIST test statistic (follows a <img src="https://render.githubusercontent.com/render/math?math=\chi^2"> distribution with 1 degree of freedom under the null), and the corresponding <img src="https://render.githubusercontent.com/render/math?math=-log_{10}(p-value)">. The last column is an indicator of whether the variant is in strong LD (|r| > 0.99) with any other variants, 0 for none and 1 for at least one.


# Parameters
| **Flag**  |**Description**                          |
|--------------------|------------|
| |**INPUT and OUTPUT files**  |
|-\-gwas-summary \<STR\>&nbsp;|Reads GWAS summary data in in GCTA-COJO format, e.g.  <br> `SNP A1 A2 freq beta se p N` <br> `rs131538 A G 0.05 0.007 0.02 0.7 6000` <br>`rs140378 C G 0.05 0.007 0.02 0.7 6000` |
|-\-bfile \<STR\>  &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;|Reads individual-level genotype data in PLINK bed format, e.g. test.fam, test.bim and test.bed|
|-\-out \<STR\>|Specifies the prefix of the output file. The output of DENTIST is in the following format, <br>`rsID        chisq   -log10(P)  ifDup` <br> `rs131538    0.012   0.91       0` <br> `rs140378    14.4    2          0` <br> `...` |
|  |**Select/filter a chromsome/region/given SNPs**  |
|-\-chrID     \<STR\>| Specifies the chromosome ID. If not provided, DENTIST will guess the chromosome ID based on if there is only one chromosome ID provided in the inputs; otherwise an error is reported.|
|-\-target    \<STR\>|  Specifies a region of 20 Mb centred at a position specified by rsID. The rsID should be present in the PLINK bim file. A warning will be reported if the target rsID is not found, and DENTIST will run across the chromosome rather than the specified region.|
|-\-target-bp \<INT\>|   has the same function as -\-target, except for that the center is defined by a bp position not rsID. |
|-\-radius    \<INT\>| Specifies the radius of a region to run DENTIST in bp unit. Together with --target or --target-bp to specify a region of interest. |
|-\-maf       \<INT\>|    Specifies a MAF threshold so that variants with MAFs smaller than this value will be precluded from the analysis. Note that this exclusion is performed after --target.|
|-\-extract   \<STR\>| Specifies a filename containing SNPs to be included for analyses. This file should contain variant identifiers (usually RS ID) in one column.|
|-\-delta-MAF  \<FLOAT\>    | Specifies a threshold for variants with MAF differences between the GWAS summary and LD reference data set. This filter is not applied by default. The commonly used thresholds include 0.1 and 0.2.|
|| **Tuning DENTIST model parameters (refer to our [paper](#Citations))** |
|-\-wind-dist    \<INT\>|  Specifies the size of a sliding window (in bp units). The default value is 2000000 bp.|
|-\-wind         \<INT\>|   has the same function as -\-wind-dist except for that this is in number of variants in a region. |
|-\-iteration-num \<INT\>|Specifies the number of iterations in the DENTIST analysis. A too large value will increase the computational costs and a too small value will increase the false discovery rate. We have experimented with this parameter and set a default value of 10 for a good trade-off.|
|-\-dup-threshold     \<FLOAT\>|    Specifies the correlation threshold for variants considered as duplicates for each other. The default value is 0.99, which corresponds to Pearson’s correlation (r2) of 0.99^2|
|-\-p-value-threshold \<FLOAT\>| Specifies the GWAS P-value threshold to group variants into null variant and significant variants. The default value is 0.01.|
|-\-with-NA-geno             | Specifies that there are NA genotypes in the input PLINK BED file. |
|-\-SVD-trunc-prop    \<FLOAT\>|    Specifies the degree SVD trunction preformed in DENTIST. It is a value between 0 and 1 specifying the proportion of components being retained from SVD truncation. The default value is 0.5. |
||**SPEED**|
|-\-thread-num \<INT\>|  Specifies the number of threads for parallel computing using OpenMP. The default value is 1.|
||**Others**|
| -\-debug             | turns on verbose output for debugging. |
|-\-write-LD                  | is yet to be implemented. |
|-\-LD-unit-in-byte \<INT\>     | is yet to be implemented. |
|-\-load-LD                   | is yet to be implemented. |
|-\-bld \<STR\>               | is yet to be implemented.|
|-\-check-LD | is yet to be implemented.  |
|-\-freq     | is yet to be implemented.  |
|-\-impute   | is yet to be implemented.  |

