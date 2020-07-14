# Overview
DENTIST (Detecting Errors iN analyses of summary staTISTics) is a quality control (QC) tool for summary-level data from genome-wide association studies (GWASs). It leverages the difference between the observed GWAS test-statistic of a variant and its predicted value (using the neighbouring variants and linkage equilibrium (LD) data from a reference panel) to remove problematic variants. It can detect genotyping/imputation errors in either the original GWAS or the LD reference samples, allelic errors (i.e., the effect alleles of the variants are mislabelled) in the GWAS summary data, as well as heterogeneity between the GWAS and LD reference samples. As shown in [our paper](#Citations), DENTIST can significantly improve the performance of the summary-data based conditional and jointly association analysis (COJO; Yang et al. 2012 Nat Genet) especially for rare variants. It can also improve the performance of LD score regression (LDSC; Builk-Sulivant et al. 2015 Nat Genet) and SMR-HEIDI (Zhu et al. 2015 Nat Genet). We believe that DENTIST can in principle work for all analyses that use GWAS summary data, although this has not been tested extensively. Apart from its application to summary data based analyses, DENTIST can also  be used in complementary to the conventional QCs in GWAS with individual-level data.

# Credits
The method is developed by Wenhan Chen, Zhihong Zhu and [Jian Yang](https://publons.com/researcher/2848531/jian-yang/) in the Yang Lab at The University of Queensland. The software is programmed and maintained by Wenhan Chen.

# Questions and Help Requests
If you meet any bugs or questions, please send an email to [Wenhan Chen](mailto:uqwche11@uq.edu.au) or [Jian Yang](mailto:jian.yang.qt@gmail.com).

# Citations
Chen et al. (2020) Improved analyses of GWAS summary statistics by eliminating heterogeneity and errors in data. bioRxiv.

# Downloads
### Pre-compiled Executable Files 
The executable file below is compiled with "-static" and tested on 64-bit Linux distributions on the x86_64 CPU platform.
Linux: [DENTIST  0.9.0.1](https://drive.google.com/file/d/1oUHc5HOTbTETLtq1ZtJMerKMpV4dSSC2/view?usp=sharing)

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
> \-\-target \<rsID\>

Specifies a region of 20 Mb centred at a position specified by rsID. The rsID should be present in the PLINK bim file. A warning will be reported if the target rsID is not found, and DENTIST will run across the chromosome rather than the specified region.

> \-\-maf \<minimal MAF threshold\>

Specifies a MAF threshold so that variants with MAFs smaller than this value will be precluded from the analysis. Note that this exclusion is performed after --target.

> \-\-wind-dist \<window size\>

Specifies the size of a sliding window (in bp units). The default value is 2000000 bp.
> \-\-thread-num \<ncpus\>

Specifies the number of threads for parallel computing using OpenMP. The default value is 1.

> \-\-num-iterations \<number of iterations\>

Specifies the number of iterations in the DENTIST analysis (see the Methods section of [our paper](#Citations)). A too large value will increase the computational costs and a too small value will increase the false discovery rate. We have experimented with this parameter and set a default value of 10 for a good trade-off.

> \-\-delta-MAF \<threshold\>

Specifies a threshold for variants with MAF differences between the GWAS summary and LD reference data set. This filter is not applied by default. The commonly used thresholds include 0.1 and 0.2.


