# QTL-seq User Guide
#### version 2.0.7

## Table of contents
- [What is QTL-seq?](#What-is-QTL-seq)
- [Installation](#Installation)
  + [Dependencies](#Dependencies)
  + [Installation using bioconda](#Installation-using-bioconda)
  + [Manual Installation](#Manual-Installation)
- [Usage](#Usage)
  + [Example 1 : run QTL-seq from FASTQ without trimming](#Example-1--run-QTL-seq-from-FASTQ-without-trimming)
  + [Example 2 : run QTL-seq from FASTQ with trimming](#Example-2--run-QTL-seq-from-FASTQ-with-trimming)
  + [Example 3 : run QTL-seq from BAM](#Example-3--run-QTL-seq-from-BAM)
  + [Example 4 : run QTL-seq from multiple FASTQs and BAMs](#Example-4--run-QTL-seq-from-multiple-FASTQs-and-BAMs)
  + [Example 5 : run QTL-plot from VCF](#Example-5--run-QTL-plot-from-VCF)
- [Outputs](#Outputs)

## What is QTL-seq?
<img src="https://210a94ef-a-fc7f2be1-s-sites.googlegroups.com/a/ibrc.or.jp/genome-e/home/bioinformatics-team/mutmap/QTL-seq_LOGO.jpg?attachauth=ANoY7coJoWWEVItHoU2_q2dDtQkdBtsHOiv228_rPBP8ZvISV9vvS7z-YRmxG5_JLCMhwPrMg92SBMU-ULQd5F1X2Tioz5Vvv_gLecE3WSMlrBYu5y69PHEuhTWWDM8hRruMJWCrVIdjeQoDL26KWbZUFVEZOHHXX5zL-s2B0UqH3zKDxurCfpmrg_gbE7y_8D9gvaGEAYe73HOR1Jl7WjpdeYijWeqanQzUnLwWzMnITpKxLzXD7fD5ebjBLI5wNZo2j7UTfRIv&attredirects=0" width=200>

Bulked segregant analysis, as implemented in  QTL-seq (Takagi et al., 2013), is a powerful and efficient method to identify agronomically important loci in crop plants. QTL-seq was adapted from MutMap to identify quantitative trait loci. It utilizes sequences pooled from two segregating progeny populations with extreme opposite traits (e.g. resistant vs susceptible) and a single whole-genome resequencing of either of the parental cultivars. While the original QTL-seq algorithm did not assume a highly heterozygous genome, a “modified QTL-seq” has been developed to handle this situation using high resolution mapping ([Itoh et al., 2019](https://doi.org/10.1007/s00122-019-03396-z)).

#### Citation
- Hiroki Takagi, Akira Abe, Kentaro Yoshida, Shunichi Kosugi, Satoshi Natsume, Chikako Mitsuoka, Aiko Uemura, Hiroe Utsushi, Muluneh Tamiru, Shohei Takuno, Hideki Innan, Liliana M. Cano, Sophien Kamoun, Ryohei Terauchi (2013).  QTL-seq: rapid mapping of quantitative trait loci in rice by whole genome resequencing of DNA from two bulked populations. Plant journal 74:174-183. [[URL]](https://doi.org/10.1111/tpj.12105)
- Yu Sugihara, Lester Young, Hiroki Yaegashi, Satoshi Natsume, Daniel J. Shea, Hiroki Takagi, Helen Booker, Ryohei Terauchi, Akira Abe (in preparation). High performance pipeline for MutMap and QTL-seq.

## Installation
### Dependencies
#### Software required
- [BWA](http://bio-bwa.sourceforge.net/)
- [SAMtools](http://samtools.sourceforge.net/)
- [BCFtools](http://samtools.github.io/bcftools/)
- [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
- [SnpEff](http://snpeff.sourceforge.net/) (optional)

The above software is included if bioconda is used to install QTL-seq. If you manually install QTL-seq the applications listed above will need to be installed and available

#### Python libraries
- matplotlib
- numpy
- pandas
- seaborn (optional)

### Installation using bioconda
You can install QTL-seq using [bioconda](https://bioconda.github.io/index.html).
```
$ conda install -c bioconda qtlseq
```

### Manual Installation
If you got a error during installation, you can install QTL-seq, using git.
```
$ git clone https://github.com/YuSugihara/QTL-seq.git
$ cd QTL-seq
$ pip install -e .
```
Manual installation of QTL-seq requires other software dependencies (listed in Software (#Software)). We highly recommend you  install SnpEff and Trimmomatic using bioconda.
```
$ conda install -c bioconda snpeff
$ conda install -c bioconda triimomatic
```
After installation, please check whether SnpEff and Trimmomatic work, for example, using the commands below.
```
$ snpEff --help
$ trimmomatic --help
```

## Usage
### File requirements
To run QTL-seq the following data files are required:
- a FASTA file containing your reference sequence e.g. myspecies_reference.fasta
- FASTQ files containing forward and reverse reads from individuals with one extreme phenotype
  e.g. resistant_bulk.1.fastq.gz and resistant_bulk.2.fastq.gz
- FASTQ files containing forward and reverse reads from individuals with the other extreme phenotype
  e.g. susceptible_bulk.1.fastq.gz and susceptible_bulk.2.fastq.gz
  
Alternatively, .bam files containing unaligned reads can be used instead of the FASTQ files
  e.g. resistant_reads.bam and susceptible_reads.bam
  
The QTL-seq pipeline can optionally identify causative SNPs and INDELs using snpEff (called with the "-e" option). This option requires a snpEFF file for your species to be avaialble to snpEff

### QTL-seq syntax

```
$ qtlseq -h

usage: qtlseq -r <FASTA> -p <BAM|FASTQ> -b1 <BAM|FASTQ>
              -b2 <BAM|FASTQ> -n1 <INT> -n2 <INT> -o <OUT_DIR>
              [-F <INT>] [-T] [-e <DATABASE>]

QTL-seq version 2.0.7

optional arguments:
  -h, --help         show this help message and exit
  -r , --ref         Reference fasta. [required]
  -p , --parent      fastq or bam of parent. If you specify
                     fastq, please separate forward and reverse 
                     pairs by commma, e.g. -p fastq1,fastq2. 
                     You can use this option multiple times. [required?]
  -b1 , --bulk1      fastq or bam of bulk 1. If you specify
                     fastq, please separate forward and reverse 
                     pairs by commma, e.g. -b1 fastq1,fastq2. 
                     You can use this option multiple times [required]
  -b2 , --bulk2      fastq or bam of bulk 2. If you specify
                     fastq, please separate forward and reverse 
                     pairs by commma, e.g. -b2 fastq1,fastq2.
                     You can use this option multiple times [required]
  -n1 , --N-bulk1    Number of individuals in bulk 1. [required]
  -n2 , --N-bulk2    Number of individuals in bulk 2. [required]
  -o , --out         Output directory. Specified name must not
                     exist and pipeline will quit if it is present. [required]
  -F , --filial      Filial generation. This parameter must be
                     greater than 1. [2]
  -t , --threads     Number of threads. If you specify a number
                     below one then QTL-seq will use as many threads 
                     as possible. [2]
  -w , --window      Window size (kb). [2000]
  -s , --step        Step size (kb). [100]
  -D , --max-depth   Maximum depth (number of reads) per variant (SNP or INDEL) 
                     to be used. [250]
  -d , --min-depth   Minimum depth of variants which will be used. [8]
  -N , --N-rep       Number of replicates for simulation to make
                     null distribution. [5000]
  -T, --trim         Trim fastq using trimmomatic.
  -a , --adapter     Path to FASTA file containing adapter sequences. Specifying 
                     -a is optional and is only used if trimming the fastq files
                     is required (that is, using "-T"). []
  --trim-params      Parameters for trimmomatic. Input parameters
                     must be separated by comma in the following
                     order: phred, ILLUMINACLIP, LEADING, TRAILING,
                     SLIDINGWINDOW, MINLEN. If you want to remove
                     Illumina adapters, please specify a FASTA file containing
                     the adapter sequences using "--adapter". The specified FASTA
                     file will be inserted into <ADAPTER_FASTA>. If adapter filtering
                     isn't specifed, adapter trimming will be skipped.
                     [33,<ADAPTER_FASTA>:2:30:10,20,20,4:15,75]
  -e , --snpEff      Predict causative variants using SnpEff. Please
                     check available databases in SnpEff.
  --mem              Maximum memory per thread when bam sorted;
                     suffix K/M/G recognized. [1G]
  -q , --min-MQ      Minimum mapping quality in mpileup. [40]
  -Q , --min-BQ      Minimum base quality in mpileup. [18]
  -C , --adjust-MQ   "adjust-MQ" in mpileup. Default parameter
                     is suited for BWA. [50]
  -v, --version      show program's version number and exit
```

QTL-seq can use FASTQ (without or with trimming) and/or BAM files for the parental, bulk1 and bulk2 reads. Once you run QTL-seq, QTL-seq automatically completes the subprocesses. That is, QTL-seq will trim (optional using "-T"), align readfiles (using bwa), filter and call variants (samtools and bcftools), calculate snp-index and delta(snp-index) values and create plots showing these values.
If you want to run QTL-seq using a VCF file previously created using the QTL-seq pipeline, please use QTL-plot (see example 5). 

+ [Example 1 : run QTL-seq from FASTQ without trimming](#Example-1--run-QTL-seq-from-FASTQ-without-trimming)
+ [Example 2 : run QTL-seq from FASTQ with trimming](#Example-2--run-QTL-seq-from-FASTQ-with-trimming)
+ [Example 3 : run QTL-seq from BAM](#Example-3--run-QTL-seq-from-BAM)
+ [Example 4 : run QTL-seq from multiple FASTQs and BAMs](#Example-4--run-QTL-seq-from-multiple-FASTQs-and-BAMs)
+ [Example 5 : run QTL-plot from VCF](#Example-5--run-QTL-plot-from-VCF)



### Example 1 : run QTL-seq using FASTQ files without trimming
```
$ qtlseq -r reference.fasta \
         -p parent.1.fastq,parent.2.fastq \
         -b1 bulk_1.1.fastq,bulk_1.2.fastq \
         -b2 bulk_2.1.fastq,bulk_2.2.fastq \
         -n1 20 \
         -n2 20 \
         -o example_dir
```

`-r` : reference fasta

`-p` : FASTQs of parent. Please supply pair-end reads separated by a comma. FASTQs can be gzipped.

`-b1` : FASTQs of bulk 1. Please supply pair-end reads separated by a comma. FASTQs can be gzipped.

`-b2` : FASTQs of bulk 2. Please supply pair-end reads separated by a comma. FASTQs can be gzipped.

`-n1` : number of individuals in bulk 1.

`-n2` : number of individuals in bulk 2.

`-o` : name of output directory. Specified name cannot exist.

### Example 2 : run QTL-seq using FASTQ files with trimming
```
$ qtlseq -r reference.fasta \
         -p parent.1.fastq,parent.2.fastq \
         -b1 bulk_1.1.fastq,bulk_1.2.fastq \
         -b2 bulk_2.1.fastq,bulk_2.2.fastq \
         -n1 20 \
         -n2 20 \
         -o example_dir \
         -T
```

`-r` : reference fasta

`-p` : FASTQs of parent. Please supply pair-end reads separated by a comma. FASTQs can be gzipped.

`-b1` : FASTQs of bulk 1. Please supply pair-end reads separated by a comma. FASTQs can be gzipped.

`-b2` : FASTQs of bulk 1. Please supply pair-end reads separated by a comma. FASTQs can be gzipped.

`-n1` : number of individuals in bulk 1.

`-n2` : number of individuals in bulk 2.

`-o` : name of output directory. Specified name cannot exist.

`-T` : trim reads using trimmomatic.

### Example 3 : run QTL-seq using BAM files
```
$ qtlseq -r reference.fasta \
         -p parent.bam \
         -b1 bulk_1.bam \
         -b2 bulk_2.bam \
         -n1 20 \
         -n2 20 \
         -o example_dir
```

`-r` : reference fasta

`-p` : BAM of parent.

`-b1` : BAM of bulk 1.

`-b2` : BAM of bulk 2.

`-n1` : number of individuals in bulk 1.

`-n2` : number of individuals in bulk 2.

`-o` : name of output directory. Specified name cannot exist.

### Example 4 : run QTL-seq using multiple FASTQs and BAMs
```
$ qtlseq -r reference.fasta \
         -p parent_1.1.fastq,parent_1.2.fastq \
         -p parent_1.bam \
         -b1 bulk_11.1.fastq,bulk_11.2.fastq \
         -b1 bulk_12.bam \
         -b1 bulk_13.bam \
         -b2 bulk_21.1.fastq,bulk_21.2.fastq \
         -b2 bulk_22.bam \
         -b2 bulk_23.bam \
         -n1 20 \
         -n2 20 \
         -o example_dir
```

QTL-seq can automatically merge multiple FASTQs and BAMs. Of course, you can manually merge FASTQs or BAMs using `cat` or `samtools merge` prior to using them in QTL-seq. If you specify `-p` multiple times, please make sure that the files are from only a single individual. On the other hand, `-b1` and `-b2` should include more than one individual because these are bulked samples. QTL-seq will assume that FASTQ files are present if a comma is present between the filenames. The pipeline will assume a BAM file is present if a comma is not present.

### Example 5 : run QTL-plot from VCF
```
$ qtlplot -h

usage: qtlplot -v <VCF> -n1 <INT> -n2 <INT> -o <OUT_DIR>
               [-F <INT>] [-t <INT>] [-w <INT>] [-s <INT>] [-D <INT>]
               [-d <INT>] [-N <INT>] [-m <FLOAT>] [-S <INT>] [-e <DATABASE>]
               [--igv] [--indel]

QTL-plot version 2.0.7

optional arguments:
  -h, --help            show this help message and exit
  -v , --vcf            VCF which contains parent, bulk1 and bulk2
                        in this order.
  -n1 , --N-bulk1       Number of individuals in bulk 1.
  -n2 , --N-bulk2       Number of individuals in bulk 2.
  -o , --out            Output directory. Specified name can exist.
  -F , --filial         Filial generation. This parameter must be
                        more than 1. [2]
  -t , --threads        Number of threads. If you specify a number
                        below one, then QTL-plot will use as many threads
                        as possible. [2]
  -w , --window         Window size (kb). [2000]
  -s , --step           Step size (kb). [100]
  -D , --max-depth      Maximum depth (number of reads) per variant (SNP or INDEL) 
                        to be used. [250]
  -d , --min-depth      Minimum depth of variants which will be used. [8]
  -N , --N-rep          Number of replicates for simulation to make
                        null distribution. [5000]
  -m , --min-SNPindex   Cutoff of minimum SNP-index for clear results. [0.3]
  -S , --strand-bias    Filter spurious homozygous genotypes in cultivar using
                        strand bias. If ADF (or ADR) is higher than this
                        cutoff when ADR (or ADF) is 0, then the SNP will be
                        filtered out. If you want to supress this filtering,
                        please set this cutoff to 0. [7]
  -e , --snpEff         Predict causative variants using SnpEff. Please
                        check available databases in SnpEff.
  --igv                 Output IGV format file to check results on IGV.
  --indel               Plot SNP-index with INDEL. Default is to plot SNPs only
  --fig-width           Width of chromosome plots. [7.5]
  --fig-height          Height of chromosome plots. [4.0]
  --white-space         White space between plots. (This option
                        only affect the vertical separation.) [0.6]
  --version             show program's version number and exit
```
QTL-plot is included in QTL-seq. QTL-plot will automatically run after QTL-seq has generated a VCF file using the default parameters. If you want to generate plots using altered parameters, use the VCF file generated by QTL-seq `(OUT_DIR/30_vcf/QTL-seq.vcf.gz)` in QTL-plot. An example is shown below.

```
$ qtlplot -v OUT_DIR/30_vcf/QTL-seq.vcf.gz \
          -o ANOTHER_DIR_NAME \
          -n1 20 \
          -n2 20 \
          -w 2000 \
          -s 100
```

#### QTL-plot using your own VCF file
It is possible to to use QTL-plot using your own VCF file. In this case, please make sure that:
1. Your VCF includes the AD field.
2. Variants were called using mpileup and three BAM alignments, parent, bulk1 and bulk2, in this order.

If you get an error, please try to run QTL-seq from FASTQ or BAM before asking in issues.

## Outputs
An example of the QTL-seq output written to `OUT_DIR` is shown below.
```
├── 10_ref
│  ├── IRGSP-1.0_genome.fasta
│  ├── IRGSP-1.0_genome.fasta.amb
│  ├── IRGSP-1.0_genome.fasta.ann
│  ├── IRGSP-1.0_genome.fasta.bwt
│  ├── IRGSP-1.0_genome.fasta.fai
│  ├── IRGSP-1.0_genome.fasta.pac
│  └── IRGSP-1.0_genome.fasta.sa
├── 20_bam
│  ├── bulk1.filt.bam
│  ├── bulk1.filt.bam.bai
│  ├── bulk2.filt.bam
│  ├── bulk2.filt.bam.bai
│  ├── parent.filt.bam
│  └── parent.filt.bam.bai
├── 30_vcf
│  ├── qtlseq.vcf.gz
│  └── qtlseq.vcf.gz.tbi
├── 40_qtlseq
│  ├── bulk1_SNPindex.png
│  ├── bulk2_SNPindex.png
│  ├── delta_SNPindex.png
│  ├── sliding_window.tsv
│  └── snp_index.tsv
└── log
   ├── bcftools.log
   ├── bgzip.log
   ├── bwa.log
   ├── samtools.log
   └── tabix.log
```
- If you run QTL-seq with trimming ("-T"), the trimmed FASTQ files will be placed in the `00_fastq` directory.
- The final results for the pipeline are in `40_QTL-seq`. Column headings for the results are:
  + `snp_index.tsv` : columns in this order.
    - **CHROM** : chromosome name
    - **POSI** : position in chromosome
    - **VARIANT** : SNP or INDEL
    - **DEPTH 1** : depth of bulk 1
    - **DEPTH 2** : depth of bulk 2
    - **p99** : 99% confidence interval for simulated delta SNP-index (absolute value)
    - **p95** : 95% confidence interval for simulated delta SNP-index (absolute value)
    - **SNP-index 1** : real SNP-index of variants in bulk 1
    - **SNP-index 2** : real SNP-index of variants in bulk 2
    - **DELTA SNP-index** : real delta SNP-index (bulk2 - bulk1)
  + `sliding_window.tsv` : columns in this order.
    - **CHROM** : chromosome name
    - **POSI** : central position of window
    - **MEAN p99** : mean 99% confidence interval for window
    - **MEAN p95** : mean 95% confidence interval for window
    - **MEAN SNP-index 1** : mean SNP-index of bulk 1 variants in window (absolute value)
    - **MEAN SNP-index 2** : mean SNP-index of bulk 2 variants in window (absolute value)
    - **MEAN DELTA SNP-index** : mean delta SNP-index
  + `QTL-seq_plot.png` : example plot shown below
    - **<span style="color: blue; ">BLUE dots</span>** : delta(SNP-index) values and location for SNPs
    - **<span style="color: red; ">RED line</span>** : mean SNP-index
    - **<span style="color: orange; ">ORANGE line</span>** : mean 99% confidence interval
    - **<span style="color: green; ">GREEN line</span>** : mean 95% confidence interval

<img src="https://user-images.githubusercontent.com/34593586/72580948-00e60500-3921-11ea-850f-dcf9a8d75e74.png" width=600>
