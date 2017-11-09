# CircRNAFisher
A systematic computational approach for de novo circular RNA identification

===============

Table of Contents
-----------------
    Prerequisites
    Installation
	Usage
	Commands&Parameter
-----------------
### Prerequisites

Perl and R are required to be installed. 
Packages "Math::CDF" and "Statistics::R" for perl are required to be installed.You can install them by:
perl -MCPAN -e shell
install Math::CDF
install Statistics::R

To take advantage of the built-in support for the Bowtie 2 alignment program and Samtools, you must have [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2) and Samtools installed.

Preparing Reference Sequences
You have to download sequences of reference genome from [UCSC Genome Bioinformatics](http://hgdownload.soe.ucsc.edu/downloads.html) and sequences of reference transcripts from [UCSC Table Browser]( http://genome.ucsc.edu/cgi-bin/hgTables?command=start).

Put all downloaded fasta files of reference genome into a folder, one file for each chromosome, filename prefix is its corresponding chromosome, end with ".fa", for instance chr1.fa. This directory is for argument "-refseq <seq> " in our pipeline.

Bowtie2 indexes
Bowtie2 indexes of whole genome and reference transcripts are required. You can download pre-built bowtie2 indexes for whole genome from http://bowtie-bio.sourceforge.net/bowtie2 and buid indexes of reference transcripts by running:
```
bowtie2-build -f reference_transcripts_seq.fa bt2-idx-trans.
```
reference_transcripts_seq.fa is the reference transcripts in fasta format,which you already have from the previous step.
These indexes are for arguments "-xg <bt2-idx>" and "-xt <bt2-idx-trans>" in our pipeline

Annotated exons
A file containing positional information for annotated exons is required.
Can be downloaded from [UCSC table browser](http://genome.ucsc.edu/cgi-bin/hgTables?command=start)(agree with the reference transcripts you downloaded from the previous step),by choosing "group:Genes and Gene Predictions","track:UCSC Genes" or "RefSeq Genes" or "Encode Genes" and "table:knownGene","output format: selected fields from primary and related tables" , then fill in a output file name in column "Output file". Choose "plain text" then enter "get output" to the next page,select "name,chrom,strand,txStart,txEnd,cdsStart,cdsEnd,exonCount,exonStarts,exonEnds,and alignID" then enter "get output". It will generate a file like this:
```
#name	chrom	strand	txStart	txEnd	cdsStart	cdsEnd	exonCount	exonStarts	exonEnds	alignID
uc001aaa.3	chr1	+	11873	14409	11873	11873	3	11873,12612,13220,	12227,12721,14409,	uc001aaa.3
uc010nxr.1	chr1	+	11873	14409	11873	11873	3	11873,12645,13220,	12227,12697,14409,	uc010nxr.1
...
```
This file is for argument "-ep <knownexons>" in our pipeline

Alternatively, you can download all pre-built indexes,reference sequences and all other files needed from our web site.
(on cluster:bowties indexes: /home/deepduoduo/projects/index/ /datamount/datas/deepduoduo/circRNA/mypipeline/mk_transcript_align/; 
hg19 sequnces /share/apps/hg19_chr. Their index_prefixes are "hg19" for hg19 genome and "UCSC_trans" for UCSC transcripts.
knownexons:/datamount/datas/deepduoduo/circRNA/hg19/UCSC_gene_exons.txt)

-----------------
### Installation 

* CircRNAFisher is writen mainly by perl script, so you can run it directly without compiling. To run CircRNAFither, just make a new folder and copy all these codes and files into it, set execute permission for everyone for all the files contained in our pipeline(by running command: chmod a+x *), then run this pipeline within this folder. 
-----------------
### Usage

* You can run CircRNAFisher step by step or put all the following commands into one file after all parameters have been set then run it together.
* If you run CircRNAFisher step by step, you can continue from any step you stoped.
-----------------
### Commands & Parameters

1 Generating unmapped reads.
```
./find_unmapped_reads.pl [options] -1 <read1> -2 <read2> -xg <bt2-idx> -xt <bt2-idx-trans> 

NAME
    find_unmapped_reads.pl

SYNOPSIS
    ./find_unmapped_reads.pl [options] -1 <read1> -2 <read2> -xg <bt2-idx>
    -xt <bt2-idx-trans>

ARGUMENTS
    -1 <read1>
        file with #1 mates, paired with file in <read2>, must be in fastq
        format, and same length or will cause issues

    -2 <read2>
        file with #2 mates, paired with file in <read1>, must be in fastq
        format, and same length or will cause issues

    -xg <bt2-idx>
        bowtie2 index filename prefix for genome alignment (minus trailing
        .X.bt2).

    -xt <bt2-idx-trans>
        bowtie2 index filename prefix for transcripts alignment (minus
        trailling .X.bt2).

OPTIONS
    -bowtie2 <path>
        The path to the Bowtie executables.(Default: the path to Bowtie
        executables is assumed to be in the user's PATH environment
        variable)

    -samtools <path>
        The path to the Samtools executables.(the path to Bowtie executables
        is assumed to be in the user's PATH environment variable)

    -p/--threads <int>
        number of alignment threads to launch (Default: 8)

    -h|--help
        Show help information.

OUTPUT
    file './output/unmapped_sample.sam'.
```

2 Generating candidate circRNA splicing sites.
```
./find_candidate_circRNAs.pl [options] -xg <bt2-idx> -refseq <seq>
NAME
    find_candidate_circRNAs.pl

SYNOPSIS
    ./find_candidate_circRNAs.pl [options] -xg <bt2-idx> -refseq <seq>

ARGUMENTS
    -xg <bt2-idx>
        bowtie2 index filename prefix for genome alignment (minus trailing
        .X.bt2).

    -refseq <seq>
        a folder contains fasta files of reference genome, one file for each
        chromosome, filename prefix is its corresponding chromosome, end
        with ".fa", for instance chr1.fa

OPTIONS
    -anchorl <int>
        Sequence length of an anchor (Default:20, better leave as default)

    -m <int>m
        maximum acceptable mismatch number of a read (Default:3)

    -anchormis <int>
        maximum acceptable mismatch number of anchors (Default:3)

    -em <int>
        maximum acceptable mismatch number during extension procedure in
        each direction (Default:3, please don't set it more than -m <int>)

    -circ_range <int>
        largest acceptable distance on genome between the two splice sites
        of circRNA (Default:2000000)

    -bowtie2 <path>
        The path to the Bowtie executables.(Default: the path to Bowtie
        executables is assumed to be in the user's PATH environment
        variable)

    -samtools <path>
        The path to the Samtools executables.(the path to Bowtie executables
        is assumed to be in the user's PATH environment variable)

    -p/--threads <int>
        number of alignment threads to launch (Default: 8)

    -h|--help
        Show help information.

OUTPUT
    file "./output/unique_spliceplot_withstand.txt".
```

3 Adjusting reads for junction-overlapping reads.
```
./circRNA_junction_read.pl -xg <bt2-idx> -refseq <seq> -ep <knownexons>
NAME
    circRNA_junction_read.pl

SYNOPSIS
    circRNA_junction_read.pl [options]-xg <bt2-idx> -refseq <seq> -ep
    <knownexons>

ARGUMENTS
    -xg <bt2-idx>
        bowtie2 index filename prefix for genome alignment (minus trailing
        .X.bt2).

    -refseq <folder>
        the folder contains fasta files of reference genome, one file for
        each chromosome, filename prefix is its corresponding chromosome,
        end with ".fa", for instance chr1.fa

    -ep <knownexons>
        File with annotated exon information.

OPTIONS
    -re-anchorl <int>
        Sequence length of an anchor for circRNA_junction_read procedure
        (Default:5, don't set it longer than the length of an anchor set in
        last step, recommand 5-10)

    -m <int>
        maximum acceptable mismatch number of a read (Default:3)

    -unique <int>
        whether only keep uniquely mapped reads. 1 for uniquely mapped
        reads, other values for all mapped reads (Default:1)

    -bowtie2 <path>
        The path to the Bowtie executables.(Default: the path to Bowtie
        executables is assumed to be in the user's PATH environment
        variable)

    -samtools <path>
        The path to the Samtools executables.(the path to Bowtie executables
        is assumed to be in the user's PATH environment variable)

    -p/--threads <int>
        number of alignment threads to launch (Default: 8)

    -h|--help
        Show help information.

OUTPUT
    file "re_unique_spliceplot_withstand.txt".
```

4 Calculating p-values and FDRs for candidate circRNAs.
```
./select_peak_circRNA.pl [options]
NAME
    select_peak_circRNA.pl

SYNOPSIS
    select_peak_circRNA.pl [options]

OPTIONS
    -fragment <int>
        fragmentsize of paired-end reads (Default:500)

    -unique <int>
        whether only keep uniquely mapped reads. 1 for uniquely mapped
        reads, other values for all mapped reads (Default:1)

    -h|--help
        Show help information.

OUTPUT
    file "./output/peak_pvalue_qvalue".
```

5 Generating final circRNAs which satisfying the FDR cur-offs.
```
./filter_by_fdr.pl [options]

NAME
    filter_by_fdr.pl

SYNOPSIS
    filter_by_fdr.pl [options]

OPTIONS
    -cut-offs <float>
        FDR cut-offs (Default:0.2)

    -h|--help
        Show help information.

OUTPUT
   file "./output/peak_fdr_passcutoff.txt".
```
-----------------
Authors
Acknowledgements
Licens
