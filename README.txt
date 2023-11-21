NGSEP - Next Generation Sequencing Experience Platform
Version 4.3.2 (19-07-2023)

===========================================================================

NGSEP provides an object model to enable different kinds of
analysis of DNA high throughput sequencing (HTS) data. The classic
use of NGSEP is a reference guided construction and downstream analysis of
large datasets of genomic variation. NGSEP performs accurate detection and
genotyping of Single Nucleotide Variants (SNVs), small and large indels, short
tandem repeats (STRs), inversions, and Copy Number Variants (CNVs). NGSEP also
provides utilities for downstream analysis of variation in VCF files, including
functional annotation of variants, filtering, format conversion, comparison,
clustering, imputation, introgression analysis and different kinds of
statistics. Version 4 includes new modules for read alignment and de-novo
analysis of short and long reads including calculations of k-mers, error
correction, de-novo analysis of Genotype-by-sequencing data and (coming soon)
de-novo assembly of long read whole genome sequencing (WGS) data.

--------------------
Building NGSEP
--------------------

NGSEP has been compiled and run successfully on the standard jdk version
11.0.8. To build the distribution library NGSEPcore.jar on a unix based
command line environment run the following commands in the directory where
NGSEPcore_4.3.2.tar.gz is located:

tar -xzvf NGSEPcore_4.3.2.tar.gz
cd NGSEPcore_4.3.2
make all

Note: Usage fields below do not include the version number. To remove the
version number, users can either copy the executable jar file:

cp NGSEPcore_4.3.2.jar NGSEPcore.jar

or just make a symbolic link:

ln -s NGSEPcore_4.3.2.jar NGSEPcore.jar

---------------
Asking for help
---------------

It is possible to obtain usage information for each module by typing:

java -jar NGSEPcore.jar <MODULE> --help

General information and the list of modules can be obtained by typing:

java -jar NGSEPcore.jar [ --help | --version | --citing ]

-------------------------------------------------------------------
-------------------------------------------------------------------
Group 1: Commands for de-novo and reference guided reads processing
-------------------------------------------------------------------
-------------------------------------------------------------------

--------------------
Demultiplexing reads
--------------------

Builds individual fastq files for different samples from fastq files of
complete sequencing lanes in which several samples were barcoded and sequenced.
Several lane files can be provided with the option -d or a single file can be
provided instead with the option -f (and -f2 for paired-end sequencing).
If neither the -d or the -f options are specified, the program tries to read
single sequencing reads from the standard input.

USAGE:

java -jar NGSEPcore.jar Demultiplex <OPTIONS>

OPTIONS:

        -i FILE		: Tab-delimited file with at least four columns by
			  default: flowcell, lane, barcode and sampleID. If
			  the -a option for dual barcode is activated, five
			  columns are expected: flowcell, lane, barcode1,
			  barcode2 and sampleID. The file must have a header
			  line. The same index file can be used to demultiplex
			  several FASTQ files (see option -d).
        -d FILE		: Tab-delimited file listing the lane FASTQ files to be
			  demultiplexed. Columns are: Flowcell, lane and fastq
			  file (which can be gzip compressed). A second fastq
			  file can be specified for pair-end sequencing. If the
			  reads sequenced for one lane are split in multiple
			  files, each file (or each pair of files) should be
			  included in a separate row. If this option is used,
			  the options -f, -f2, -c and -l are ignored.
        -o DIR		: Directory where the output fastq files will be saved.
			  Files will be gzip compressed by default.
        -f FILE		: File with raw reads in fastq format. It can be gzip
			  compressed.
        -f2 FILE	: File with raw reads in fastq format corresponding to
			  the second file for paired end reads. It can be gzip
			  compressed.
        -c STRING	: Id of the flowcell corresponding to the input fastq
			  file(s). Ignored if the -d option is specified but
			  required if -d option is not specified.
        -l STRING	: Id of the lane corresponding to the input fastq
			  file(s). Ignored if the -d option is specified but
			  required if the -d option is not specified.
        -t STRING	: Sequences to trim separated by comma. If any of the
			  given sequences is found within a read, the read will
			  be trimmed up to the start of the sequence.
        -u		: Output uncompressed files.
        -r INT		: Minimum read length to keep a read after trimming
			  adapter sequences. Default: 40.
        -a		: Activate demultiplexing with dual barcoding.


----------------------------------------
Obtaining k-mers spectrum from sequences
----------------------------------------

Extracts k-mers and generates a distribution of k-mer abundances from a file of
DNA sequences either in fastq or in fasta format (see -f option). Writes two
files, one with the k-mer distribution and a second file with the actual k-mers
and their counts. 

USAGE:

java -jar NGSEPcore.jar KmersExtractor <OPTIONS> <SEQUENCES_FILE>*

OPTIONS:

	-o FILE	: Prefix of the output files.
	-k INT		: K-mer length. Default: 15
	-m INT		: Minimum count to report a k-mer in the output file.
			  Default: 5
	-text		: Indicates that the sequences should be treated as
			  free text. By default it is assumed that the given
			  sequences are DNA and then only DNA k-mers are
			  counted. If this option is set, the -s option is also
			  activated to process the text only in the forward
			  direction.
	-s		: If set, only the forward strand would be used to
			  extract kmers. Mandatory for non-DNA sequences.
	-f INT		: Format of the input file(s). It can be 0 for fastq or
			  1 for fasta. Default: 0
	-c		: Ignore low complexity k-mers for counting and reporting.
	-t INT		: Number of threads. Default: 1


------------------------
Fixing sequencing errors
------------------------

Builds a k-mer abundance profile and use this profile to identify and correct
sequencing errors. For each predicted single nucleotide error, it looks for the
single change that would create k-mers within the normal distribution of
abundances. Using the option -e, this function can also receive a precalculated
table of k-mers, which could come from a larger number of reads or reads
sequenced using a different technology. For example, a k-mers profile based on
Illumina reads could be built using the KmersExtractor command, and then this
profile could be used to perform error correction on long reads.

USAGE:

java -jar NGSEPcore.jar ReadsFileErrorsCorrector <OPTIONS>

OPTIONS:

	-i FILE	: Input file with raw reads in fastq or fasta format. See
		  option -f for options on the file format. It can be gzip
		  compressed.
	-o FILE	: Output file with the corrected reads in fastq format
		  (gzip compressed).
	-e FILE	: Two column tab delimited file with k-mers and their
		  abundances.
	-k INT	: K-mer length. Default: 15
	-m INT	: Minimum k-mer count to consider a k-mer real. Default: 5
	-s	: If set, only the forward strand would be used to extract
		  kmers. Mandatory for non-DNA sequences.
	-f INT	: Format of the input file. It can be 0 for fastq or 1 for
		  fasta. Default: 0

----------------------------------------
Performing de-novo analysis of GBS reads
----------------------------------------

Performs de novo variants discovery from a genotype-by-sequencing (GBS) or
a double digestion RAD sequencing (ddRAD) experiment. Runs a clustering
algorithm based on quasi-exact matches to representative k-mers within the
first base pairs of each sequence. Then, it performs variants detection and
sample genotyping within each cluster using the same Bayesian model
implemented for the reference-guided analysis. By now it can only discover and
genotype Single Nucleotide Variants (SNVs).

USAGE:

java -jar NGSEPcore.jar DeNovoGBS <OPTIONS>

OPTIONS:

	-i FILE         : Directory with fastq files to be analyzed. Unless the
			  -d option is used, it processes as single reads all
			  fastq files within the given directory.
	-o FILE         : Prefix for the output VCF file with the discovered
			  variants and genotype calls as well as other output
			  files describing the behavior of this process.
	-d FILE         : Tab delimited text file listing the FASTQ files to be
			  processed for paired-end sequencing. It should have
			  three columns. sample id, first fastq file and second
			  fastq file. All files should be located within the
			  directory provided with the option -i.
	-k INT          : K-mer length. Default: 31
	-c INT          : Maximum number of read clusters to process. This
			  parameter controls the amount of memory spent by the
			  process. Default: 2000000
	-t INT          : Number of threads to process read clusters. Default: 1
	-maxBaseQS INT  : Maximum value allowed for a base quality score.
			  Larger values will be equalized to this value.
			  Default: 30
	-ignore5 INT	: Ignore this many base pairs from the 5' end of the
			  reads. Default: 0
	-ignore3 INT	: Ignore this many base pairs from the 3' end of the
			  reads. Default: 0
	-h DOUBLE       : Prior heterozygosity rate. Default: 0.001
	-minQuality INT : Minimum variant quality. In this command, this filter
			  applies to the QUAL column of the VCF, which is
			  calculated for each variant as the maximum of the
			  genotype qualities of samples with non-homozygous
			  reference genotype calls. See the command VCFFilter
			  to apply filters of quality and read depth on
			  individual genotype calls. Default: 40
	-ploidy INT     : Default ploidy of the samples. Default: 2

----------------------------------
Assembling genomes from long reads
----------------------------------

Builds a de-novo assembly from a set of long reads following an
overlap-layout-consensus (OLC) approach. It receives a fasta or fastq file with
raw PacBio HiFi or Nanopore reads and generates an assembly for the sample in
fasta format. It also generates a grap.gz file with the information of the
overlap graph. This graph can be provided as input in a second run skip the
graph construction step. Note: Although this functionality is already producing
competitive assemblies compared to other tools,  the following versions will
probably have improvements on contiguity and phasing of diploid assemblies.

USAGE:

java -jar NGSEPcore.jar Assembler <OPTIONS>

OPTIONS:

	-i FILE	: Input file. See option -f for options on the file
			  format. It can be gzip compressed.
	-o FILE	: Prefix of the output files.
	-g FILE	: File with a saved graph to perform layout and
			  consensus.
	-k INT		: K-mer length to identify overlaps Default: 15
	-f INT		: Format of the input file. It can be 0 for fastq or
			  1 for fasta. Default: 0
	-w INT		: Window length to calculate minimizers. Default: 30
	-m INT		: Minimum read length. Default: 5000
	-mspe DOUBLE	: Minimum proportion from the maximum score of the
			  edges of a sequence to keep an edge. Default: 0.3
	-ploidy INT	: Ploidy of the sample. Keep ploidy of 1 if the sample
			  is inbred, even if it is diploid or polyploid. This
			  option is still in progress and it has been tested
			  only in haploid and diploid samples. Default: 1
	-cml INT	: Maximum length of circular molecules Default: 0
	-cmof FILE	: Fasta file with known start sequences of circular
        		  molecules
	-ac STRING	: Algorithm used to build the consensus. It can be
			  Simple or Polishing. Default: Polishing
	-t INT		: Number of threads. Default: 1


------------------------------
Updating genomes from variants
------------------------------

Takes a VCF file with genotype information from one sample and the reference
genome used to build the VCF and generates a new genome in fasta format with
a ploidy consistent with the ploidy of the individual. For diploid or
polyploid assemblies, the VCF file must be properly phased. For non haploid
individuals, if the ploidy parameter is set to 1, this function performs
polishing of a haploid genome assembly.

USAGE:

java -jar NGSEPcore.jar IndividualGenomeBuilder <OPTIONS>

OPTIONS:

	-i FILE		: Fasta file with the original genome.
	-v FILE		: File in VCF format with the variants that will be
			  applied to the input genome.
	-ploidy INT	: Ploidy of the sample. To make polishing of a haploid
			  assembly for a non haploid individual, set this
			  parameter to 1. Default: 2
	-o FILE		: Output file in fasta format with the modified genome.

-------------------------------
Indexing genome reference files
-------------------------------

Creates a binary file containing an FM index for large sequences in fasta format
(usually a reference genome). This structure facilitates performing massive text
searches over the indexed sequence. This is a usual preparation step for
alignment of short reads.

USAGE:

java -jar NGSEPcore.jar GenomeIndexer <OPTIONS>

OPTIONS:

	-i FILE	: Input genome to index in fasta format.
		  It can be gzip compressed.
	-o FILE	: Output binary file with the FM index associated with the
		  input genome.

-----------------------------------
Aligning reads to reference genomes
-----------------------------------

Calculates a list of genomic regions for sites where the reads can be found in
a reference genome. It receives up to two files with raw reads in fastq format
and the reference genome. To map short reads to long genomes, a precalculated
FM index can also be provided with the option -d. See command GenomeIndexer for
construction of the FM index. It provides as output a file with alignments to
the reference genome in BAM format.

USAGE:

java -jar NGSEPcore.jar ReadsAligner <OPTIONS>

OPTIONS:

	-i FILE		: Input file with raw reads in fastq format. It can be
			  gzip compressed. Required if the second fastq file is
			  provided using the option -i2.
	-i2 FILE	: Input file with raw reads in fastq format
			  corresponding to the second file for paired end reads.
			  It can be gzip compressed.
	-o FILE		: Output file with the aligned reads in BAM format.
	-r GENOME	: Reference genome to align reads in FASTA format.
			  Required parameter. It can be gzip compressed.
	-d FILE		: FM-index of the reference genome to align short reads.
			  See GenomeIndexer for instructions to generate this
			  file. For large genomes it is more efficient to index
			  the reference once and provide the index with this
			  option.
	-s STRING	: Id of the sample. Default: Sample
	-p STRING	: Sequencing platform used to produce the reads.
			  Supported platforms include ILLUMINA, IONTORRENT,
			  PACBIO and ONT. Default: ILLUMINA
	-knownSTRs FILE	: Text file with location of known short tandem repeats
			  (STRs). It is a tab-delimited file with at least
			  three columns: Sequence name (chromosome), region
			  first base pair coordinate (1-based, inclusive) and
			  region last base pair coordinate (1-based, inclusive).
	-f INT		: Format of the input file. It can be 0 for fastq or 1
			  for fasta. Default: 0
	-k INT		: K-mer length. Default: 25
	-m INT		: Maximum alignments per read. Default: 3
	-minIL INT	: Minimum predicted insert length to consider an
			  alignment proper. Default: 0
	-maxIL INT	: Maximum predicted insert length to consider an
			  alignment proper. Default: 1000
	-w INT		: Window length to compute minimizers. Default: 20
	-t INT		: Number of threads used to align reads. Default: 1




-------------------------------------------------------
-------------------------------------------------------
Group 2: Commands for variants discovery and genotyping
-------------------------------------------------------
-------------------------------------------------------

----------------------------------------
Calculating base pair quality statistics
----------------------------------------

Takes one or more sets of alignments and a reference genome and counts the
number of mismatches with the reference for each read position from 5' to 3'
end. This report is useful to detect sequencing error biases. Requires one or
more alignment files in SAM, BAM or CRAM format, and the reference genome that
was used to produce the alignments. Writes to standard output unless the -o
option is used to specify an output file.

USAGE:

java -jar NGSEPcore.jar BasePairQualStats <OPTIONS> <ALIGNMENTS_FILE>*

OPTIONS:

	-o FILE		: Output file with the base pair quality statistics.
	-r FILE		: Fasta file with the reference genome.
	-minMQ INT	: Minimum mapping quality to call an alignment unique.
			  Default: 20
						  
The file(s) with alignments must be given in SAM, BAM or CRAM format and the
reference file in fasta format. The output is a text file with five columns:
- Position: 1- based from 5' to 3'
- Number of reads with a base call different than the reference (Considering
  all alignments)
- Number of reads with a base call different than the reference (Considering
  only reads with unique alignments)
- Number of total alignments counted with read length equal or larger than the
  position in the first column. The percentage of mismatches including all
  alignments is the ratio of column 2 divided by this column
- Number of uniquely aligned reads counted with read length equal or larger
  than the position in the first column. The percentage of mismatches for
  uniquely aligned reads is the ratio of column 3 divided by this column

-------------------------------
Calculating coverage statistics
-------------------------------

Calculates the number of base pairs that are covered by reads at each read
depth level from 1 to a maximum. Alignments must be in SAM, BAM or CRAM
format. Writes to standard output unless the -o option is used to specify an
output file.

USAGE:

java -jar NGSEPcore.jar CoverageStats <OPTIONS>

OPTIONS:

	-i FILE		: Input file with alignments to analyze.
	-o FILE		: Output file with the coverage distribution.
	-r GENOME	: Fasta file with the reference genome. Required for
			  CRAM files.
	-minMQ INT	: Minimum mapping quality to call an alignment unique.
			  Default: 20

The alignments file must be given in SAM or BAM format. The output is a text
file with three columns:
- Coverage
- Number of reference sites with this coverage (Considering all alignments)
- Number of reference sites with this coverage (Considering only reads with 
  unique alignments)

--------------------------------------
Calling variants over multiple samples
--------------------------------------

This module allows to call variants over a group of samples separated by files
or read group tags. This is now the recommended method to perform variants
detection on genotype-by-sequencing (GBS), RAD sequencing, whole exome
sequencing (WES), RNA-seq and low coverage (less than 10x) whole genome
sequencing (WGS) data. Although it can also be used on high coverage WGS data,
the classic sample-by-sample analysis (commands SingleSampleVariantsDetector,
MergeVariants and VCFMerge) is still recommended to identify structural
variants. This module requires one or more read alignment files in SAM, BAM or
CRAM format and the reference genome that was used to produce the alignments.
Alignments must be sorted by reference coordinates.

USAGE:

java -jar NGSEPcore.jar MultisampleVariantsDetector <OPTIONS> <ALIGNMENTS_FILE>*

OPTIONS:

	-r GENOME		: Fasta file with the reference genome.
	-o FILE			: Output VCF file with discovered variants and
				  genotype calls. Default: variants.vcf
	-ploidy INT             : Default ploidy of the samples. Default: 2
	-psp                    : Print id and ploidy of the sample in the VCF
				  header. The header generated with this option
				  is not a standard VCF header. However, it
				  helps NGSEP to keep track of the ploidy of
				  the samples through downstream analyses.
	-minMQ INT              : Minimum mapping quality to call an alignment
				  unique. Default: 20
	-knownVariants FILE     : VCF file with variants to be genotyped. Only
				  these variants will appear in the output VCF.
	-querySeq STRING	: Call variants just for this sequence.
	-first INT		: Call variants just from this position in the
				  given query sequence.
	-last INT		: Call variants just until this position in the
				  given query sequence.
	-ignoreLowerCaseRef	: Ignore sites where the reference allele is
				  lower case.
	-maxAlnsPerStartPos INT	: Maximum number of alignments allowed to start
				  at the same reference site. This parameter
				  helps to control false positives produced by
				  PCR amplification artifacts. In this command,
				  this filter is executed independently for
				  each read group. For GBS or RAD sequencing
				  data use a large value such as 100. Default: 5
	-p			: Process non unique primary alignments in the
				  pileup process. The default behavior is to
				  process alignments that are unique
				  (see option -minMQ).
	-s			: Consider secondary alignments in the pileup
				  process. Non-unique primary alignments will
				  also be considered in this mode.
	-ignore5 INT            : Ignore this many base pairs from the 5' end of
				  the reads. Default: 0
	-ignore3 INT            : Ignore this many base pairs from the 3' end of
				  the reads. Default: 0
	-h DOUBLE		: Prior heterozygosity rate. Default: 0.001
	-maxBaseQS INT          : Maximum value allowed for a base quality
				  score. Larger values will be equalized to
				  this value. This parameter allows to reduce
				  the effect of sequencing errors with high
				  base quality scores. Default: 30
	-knownSTRs FILE         : File with known short tandem repeats (STRs).
				  This is a text file with at least three
				  columns: chromosome, first position and last
				  position. Positions should be 1-based and
				  inclusive.
	-minQuality INT         : Minimum variant quality. In this command,
				  this filter applies to the QUAL column of the
				  VCF, which is calculated for each variant as
				  the maximum of the genotype qualities of
				  samples with non-homozygous reference
				  genotype calls. See the command VCFFilter to
				  apply filters of quality and read depth on
				  individual genotype calls. Default: 40
	-embeddedSNVs           : Flag to call SNVs within STRs. By default,
				  STRs are treated as a single locus and hence
				  no SNV will be called within an STR.

Alignments should be provided in SAM, BAM or CRAM format
(see http://samtools.github.io/hts-specs for details).
The reference genome should be provided in fasta format. It is assumed that
the sequence names in the alignments file correspond with the sequence names in
this reference assembly. For this module, alignment files must include RG headers
including the ID of each read group and the corresponding sample. A typical read
group header looks as follows

@RG	ID:<ReadGroupId>	SM:<SampleId>	PL:<Platform>

According to the specification, read groups must be unique across different
samples. Each alignment must include an RG tag indicating the id of the read
group of the aligned read. See the SAM format specification for details.
Check the documentation of your read aligner to make sure that alignment files
contain read group headers. This module uses read group headers to distribute the reads
that belong to the different samples.

DETAILS OF OUTPUT FILES: The output of this module is a VCF file
(see the standard format at http://samtools.github.io/hts-specs/VCFv4.2.pdf).
NGSEP uses standard output file formats such as VCF
and GFF to facilitate integration with other tools and visualization in
web genome browsers such as jbrowse, gbrowse and the UCSC genome browser, or
desktop browsers such as the Integrative Genomics Viewer (IGV). This allows
integrated visuallzation of read alignments, variants and functional elements.
Moreover, the NGSEP output files provide as optional fields custom
information on the variants and genotype calls. For each variant, NGSEP VCF
files include the following custom fields in the INFO column
(described also in the VCF header): 

TYPE (STRING)	: Type of variant for variants different than biallelic SNVs.
		  Possible types include MULTISNV, INDEL and STR (short tandem
		  repeat). Also, SNVs called within INDELS or STRs are tagged
		  with the EMBEDDED type.
CNV (INT)	: Number of samples with CNVs covering this variant. This will
		  not be generated by this command but by the classical per
		  sample analysis, and only if the read depth analysis is
		  executed
NS (INT)	: Number of samples genotyped.
MAF (DOUBLE)	: Minor allele frequency. Calculated by the VCFMerge and the
		  VCFFilter commands
AN (INT)	: Number of different alleles observed in called genotypes.
AFS (INT*)	: Allele counts over the population for all alleles, including
		  the reference. One number per allele.
		

Additionally, the Annotate command adds the INFO fields TA, TID, TGN, TCO and
TACH with the results of the functional annotation (See the Annotate command
below for details). NGSEP VCFs also include custom format fields with the
following information for each genotype call:

ADP (INT,INT,...)	: Number of base calls (depth) for alleles, including
			  the reference allele. The order of the counts is,
			  first the depth of the reference allele and then the
			  read depths of the alternative alleles in the order
			  listed in the ALT field. 
BSDP (INT,INT,INT,INT)	: For SNVs, number of base calls (depth) for the 4
			  possible nucleotides, sorted as A,C,G,T.
ACN (INT, INT, ...)	: Predicted copy number of each allele taking into
			  account the prediction of number of copies of the
			  region surrounding the variant. The order is the same
			  that in the ADP field

WARNING 1: Since v3.2.0, for RAD Sequencing or genotype-by-sequencing (GBS)
we now recommend to perform variants detection and genotyping using this module.
However, using the default value of the parameter to control for PCR duplicates
(maxAlnsPerStartPos) will yield very low sensitivity. We recommend to increase
the value of the parameter to about 100 to retain high sensitivity while
avoiding a severe penalty in memory usage. The default usage for RAD-Seq or GBS
samples becomes:

java -jar NGSEPcore.jar MultisampleVariantsDetector -maxAlnsPerStartPos 100 -r <REFERENCE> -o <OUTPUT_VCF> <BAM_FILE>*

WARNING 2: Unlike the behavior of the classical individual analysis per sample,
in this command the filter executed using the minQuality option applies to the
QUAL field of the VCF format, which corresponds to the probability of existence
of each variant (regardless of the individual genotype calls). In this module
the QUAL probability is calculated for each variant as the maximum of the
genotype qualities of samples with non-homozygous reference genotype calls. The
rationale for this calculation is that one variant should be real if it is
confidently called in at least one sample. Individual genotype calls are not
filtered by default and hence they could include some false positives. Please see
the VCFFilter command to perform custom filtering of genotype calls, either by
genotype quality (GQ format field) or by read depth (BSDP and ADP format fields).
Default values of other parameters are also set to maximize sensitivity. For
conservative variant detection including control for errors in base quality
scores and PCR amplification artifacts use:

java -jar NGSEPcore.jar MultisampleVariantsDetector -maxAlnsPerStartPos 2 -r <REFERENCE> -o <OUTPUT_VCF> <BAM_FILE>*

If the error rate towards the three prime end increases over 2% you can also
use the option -ignore3 to ignore errors at those read positions. If the
reference genome has lowercase characters for repetitive regions (usually
called softmasked), these regions can be directly filtered using the option
-ignoreLowerCaseRef. These regions can also be filtered at later stages of
the analysis using the VCFFilter command.

-----------------------------------------------------------------
Calling variants on individual samples with the variants detector
-----------------------------------------------------------------

This is the classic module of NGSEP to call SNVs, small indels and structural
variants from sequencing data of single individuals. Basic usage requires an
alignments file in SAM, BAM or CRAM format, the reference genome that
was used to produce the alignments, and a prefix for the output files.
Alignments must be sorted by reference coordinates.

USAGE: 

java -jar NGSEPcore.jar SingleSampleVariantsDetector <OPTIONS>

OPTIONS:

	-i FILE		: Input file with read alignments.
	-r FILE		: Fasta file with the reference genome.
	-o FILE		: Prefix for the output files.
	-sampleId STRING	: Id of the sample for the VCF file. If not set
				  it looks in the BAM file header for an SM
				  header tag. If this tag is not present, it
				  uses the default value. Default: Sample
	-ploidy INT		: Ploidy of the sample to be analyzed.
				  Default 2
	-psp			: Flag to print a header in the VCF file with
				  the id and the ploidy of the sample. The
				  header generated with this option is not a
				  standard VCF header. However, it helps NGSEP
				  to keep track of the ploidy of each sample
				  through downstream analyses.
	-minMQ INT		: Minimum mapping quality to call an alignment
				  unique. Default: 20
	-knownVariants FILE	: VCF file with variants to be genotyped. Only
				  these variants will appear in the output vcf
				  file. With this option homozygous calls to
				  the reference allele will be reported
	-querySeq STRING	: Call variants just for this sequence.
	-first INT		: Call variants just from this position in the
				  given query sequence
	-last INT		: Call variants just until this position in the
				  given query sequence
	-ignoreLowerCaseRef	: Ignore sites where the reference allele is
				  lower case. 
	-maxAlnsPerStartPos INT: Maximum number of alignments allowed to start
				  at the same reference site. This parameter
				  helps to control false positives produced by
				  PCR amplification artifacts. For GBS or RAD
				  sequencing data use a large value such as 100.
				  Default 5.
	-p			: Process non unique primary alignments in the
				  pileup process.  The default behavior is to
				  process alignments that are unique
				  (see option -minMQ) 
	-s			: Consider secondary alignments while calling
				  SNVs. Non-unique primary alignments will also
				  be considered in this mode.
	-ignore5 INT		: Ignore this many base pairs from the 5' end
				  of the reads. Default: 0
	-ignore3 INT		: Ignore this many base pairs from the 3' end
				  of the reads. Default: 0
	-h FLOAT		: Prior heterozygosity rate. Default: 0.001
	-maxBaseQS INT		: Maximum value allowed for a base quality
				  score. Larger values will be equalized to
				  this value. This parameter allows to reduce
				  the effect of sequencing errors with high
				  base quality scores. Default: 30
	-knownSTRs FILE	: File with known short tandem repeats (STRs).
				  This is a text file with at least three
				  columns, chromosome, first position and last
				  position. Positions should be 1-based and
				  inclusive.
	-minQuality INT	: Minimum genotype quality to accept a SNV call
				  Genotype quality is calculated as 1 minus the
				  posterior probability of the genotype given
				  the reads (in phred scale). Default: 0
	-embeddedSNVs		: Flag to call SNVs within STRs. By default,
				  STRs are treated as a single locus and hence
				  no SNV will be called within an STR.
	-csb			: Calculate a exact fisher test p-value for
				  strand bias between the reference and the
				  alternative allele
	-knownSVs FILE		: File with coordinates of known structural
				  variants in GFF format.
	-minSVQuality INT	: Minimum quality score (in PHRED scale) for
				  structural variants. Default: 20
	-runRep		: Turns on the procedure to find repetitive
				  regions based on reads with multiple
				  alignments.
	-runRD			: Turns on read depth (RD) analysis to identify
				  CNVs
	-genomeSize INT	: Total size of the genome to use during
				  detection of CNVs. This should be used when
				  the reference file only includes a part of
				  the genome (e.g. a chromosome or a partial
				  assembly).
	-binSize INT		: Size of the bins to analyze read depth.
				  Default: 100
	-algCNV STRING		: Comma-separated list of read depth algorithms
				  to run (e.g. CNVnator,EWT). Default: CNVnator
	-maxPCTOverlapCNVs INT	: Maximum percentage of overlap of a new CNV
				  with an input CNV to include it in the output
				  Default: 100 (No filter)
	-runRP			: Turns on read pair plus split-read analysis
				  (RP+SR) to identify large indels and
				  inversions.
	-maxLenDeletion INT	: Maximum length of deletions that the
				  RP analysis can identify. Default: 1000000
	-sizeSRSeed INT	: Size of the seed to look for split-read
				  alignments. Default: 8
	-ignoreProperPairFlag	: With this option, the proper pair flag will
				  not be taken into accout to decide if the ends
				  of each fragment are properly aligned. By
				  default, the distribution of insert length is
				  estimated only taking into account reads with
				  the proper pair flag turned on.
	-runOnlySVs		: Turns off SNV detection. In this mode, only
				  structural variation will be called
	-runLongReadSVs	: Runs the DBScan algorithm to identify structural
				  variants from alignments of long reads

Alignments should be provided in SAM, BAM or CRAM format
(see http://samtools.github.io/hts-specs for details).
The reference genome should be provided in fasta format. The output for SNVs,
small indels and STRs is a VCF file. These standard formats are used to
facilitate integration with other tools. See more details above in the
description of the MultisampleVariantsDetector command.
Structural variants are reported in a gff file (see the standard format at 
http://www.sequenceontology.org/gff3.shtml). This file can be used as a
parameter of the variants detector (option "-knownSVs") which is useful to
avoid recalculation of structural variants while genotyping known variants.
GFF files provided by NGSEP include the following INFO fields:

LENGTH (INT)	: Predicted length of the event. For insertions and deletions
		  identified with read pair analysis, this length is not the
		  reference span but the average of the lengths predicted by
		  each read pair having an alignment with a predicted insert
		  length significantly larger or shorter than the average
		  fragment length.
SOURCE (STRING)	: Algorithm that originated each variant call. Current values
		  include MultiAlns for repeats, Readpairs for read pair
		  analysis and CNVnator and EWT for read depth analysis.
NSF (INT)	: Number of fragments supporting the structural variation
		  event. For read depth algorithms is the (raw) number of reads
		  that can be aligned within the CNV. For read pair analysis
		  is the number of fragments (read pairs) that support the
		  indel or the inversion. For repeats is the number of reads
		  with multiple alignments
NC (DOUBLE)	: For CNVs called with the read depth algorithms this is the
		  estimated number of copies. It is kept as a real number
		  to allow users to filter by proximity to an integer value if
		  needed. 
HET (INT)	: For CNVs called with the read depth algorithms this is the
		  number of heterozygous genotype calls in the VCF file
		  enclosed within the CNV. Always zero if the option -noSNVS
		  is used.
NTADF (INT)	: For CNVs called with the read depth algorithms this is the
		  number of paired-end fragments showing an alignment pattern
		  consistent with a tandem duplication
NTRDF (INT)	: For CNVs called with the read depth algorithms this is the
		  number of paired-end fragments in which one read aligns
		  within the CNV and its pair aligns to another chromosome or
		  with a very long insert length. These fragments are useful
		  to classify the CNV as an interspersed (trans) duplication.
TGEN (STRING)	: For CNVs called with the read depth algorithms this is a
		  qualitative evaluation of the genotype call based on the
		  values of the fields NC, NTADF and NTRDF, and on the normal
		  ploidy of the sample. Possible values are DEL, TANDEMDUP
		  and TRANSDUP.
NSR (INT)	: Number of reads with split alignments supporting an insertion
		  or deletion. Events supported only by split-read analysis
		  have NSF=0 and NSR>0.
NUF (INT)	: For repeats identified from reads aligning to multiple
		  locations, this is the number of fragments with unique
		  alignments within the repeat. 

The default minimum genotype quality of the variants detector (0) will maximize
the number of called variants at the cost of generating some false positives in
samples with small coverage or high sequencing error rates. For conservative
variant calling from whole genome sequencing reads use:

java -jar NGSEPcore.jar SingleSampleVariantsDetector -maxAlnsPerStartPos 2 -minQuality 40 -r <REFERENCE> -i <INPUT_FILE> -o <OUTPUT_PREFIX>

If interested in structural variation, for Illumina reads you can add the
options to run read depth (RD) and read pair plus split read (RP+SR)
approaches to identify structural variation:

java -jar NGSEPcore.jar SingleSampleVariantsDetector -runRD -runRP -maxAlnsPerStartPos 2 -minQuality 40 -r <REFERENCE> -i <INPUT_FILE> -o <OUTPUT_PREFIX>

To detect structural variants from long reads, add the option -runLongReadSVs

java -jar NGSEPcore.jar SingleSampleVariantsDetector -runLongReadSVs -maxAlnsPerStartPos 2 -minQuality 40 -r <REFERENCE> -i <INPUT_FILE> -o <OUTPUT_PREFIX>

If the error rate towards the three prime end increases over 2% you can also
use the option -ignore3 to ignore errors at those read positions. If the
reference genome has lowercase characters for repetitive regions (usually
called softmasked), these regions can be directly filtered using the option
-ignoreLowerCaseRef. These regions can also be filtered at later stages of
the analysis using the VCFFilter command.

Since v 3.2.0, for RAD Sequencing or genotype-by-sequencing (GBS) we now
recommend the MultisampleVariantsDetector command described above. However,
the classic per-sample analysis pipeline using this command is still functional
with good quality. For both commands it is important to keep in mind that using
the default value of the parameter to control for PCR duplicates
(maxAlnsPerStartPos) will yield very low sensitivity. We recommend to increase
the value of the parameter to about 100 to retain high sensitivity while
avoiding a severe penalty in memory usage. Also, structural variants should not
be called using these data. The usage for conservative variant calling in
RAD-Seq or GBS samples becomes:

java -jar NGSEPcore.jar SingleSampleVariantsDetector -maxAlnsPerStartPos 100 -minQuality 40 -r <REFERENCE> -i <INPUT_FILE> -o <OUTPUT_PREFIX>

This module can also be used to discover variants at different allele frequencies
in pooled samples. For this use case, the ploidy parameter should be adjusted to
the number of haplotypes present within the pool. The prior heterozygosity rate
should also be increased according to the expected proportion of heterozygous
variants across the entire population. For targeted sequencing, the 
maxAlnsPerStartPos parameter should also be adjusted to make it larger than the
expected maximum read depth per site. For example, if 50 diploid individuals
were included in a pool and sequenced at 20x per individual, then this parameter
should be larger than 1000. This is a usage example to identify low frequency
variants from a pool of 48 individuals in a Tilling experiment:

java -jar NGSEPcore.jar SingleSampleVariantsDetector -maxAlnsPerStartPos 5000 -h 0.1 -ploidy 96 -r <REFERENCE> -i <INPUT_FILE> -o <OUTPUT_PREFIX>

-------------------------------------------
Molecular haplotyping of single individuals
-------------------------------------------

Performs molecular haplotyping on a single individual given a VCF and a set of
alignments in SAM, BAM or CRAM format. Although theoretically it can work with
Illumina reads, it is designed to work fine with long (PacBio) reads.
Alignments must be sorted by reference coordinates. Writes to standard output
unless the -o option is used to specify an output file.

USAGE:

java -jar NGSEPcore.jar SIH <OPTIONS>

OPTIONS:

	-i FILE		: Input VCF file with variants to phase.
	-b FILE		: Input file with read alignments.
	-o FILE		: Output VCF file with phased variants.
	-a STRING	: Algorithm for single individual haplotyping. It can
			  be Refhap or DGS. Default: DGS
	-minMQ INT	: Minimum mapping quality to call an alignment unique.
			  Default: 20
	-r GENOME	: Fasta file with the reference genome. Required for
			  CRAM files.

----------------------------------------
Merging variants from individual samples
----------------------------------------

The classical pipeline of NGSEP performs two steps to merge variants and
genotype calls from different samples into an integrated VCF file. After
independent variant discovery using the command SingleSampleVariantsDetector
for each sample, the next step is to generate a file including the whole set of
variants called in at least one of the samples. This can be done calling the
MergeVariants command, which has the following usage:

USAGE:

java -jar NGSEPcore.jar MergeVariants <OPTIONS> <VARIANTS_FILE>*

OPTIONS:

	-s FILE	: List of sequence names as they appear in the original
		  reference genome
	-o FILE	: Output VCF file with merged variants

The sequence names file is a text file which just has the ids of the sequences
in the reference. It is used by the program to determine the order of the
reference sequences. In unix systems this file can be obtained running the
following command on the fasta file with the reference genome:

awk '{if(substr($1,1,1)==">") print substr($1,2) }' <REFERENCE_FILE> > <SEQUENCE_NAMES_FILE>

If samtools is available. The fai index file provided by this tool can also be
used as a sequence names file. The fai index is generated with this command:

samtools faidx <REFERENCE_FILE>

The output file of the merge program is a vcf with the union of variants
reported by the input files but without any genotype information. 

The next step is to genotype for each sample the variants produced by
MergeVariants using the variants detector (See SingleSampleVariantsDetector).
For each sample, the command to execute at this stage (in conservative mode)
should look like this:

java -jar NGSEPcore.jar SingleSampleVariantsDetector -maxAlnsPerStartPos 2 -minQuality 40 -knownVariants <VARS_FILE> -r <REFERENCE> -i <INPUT_FILE> -o <OUTPUT_PREFIX>

where VARS_FILE is the output file obtained in the first step of the merging
process. At the end, this will produce a second set of vcf files which will
differ from the first set in the sense that they will include calls to the
reference allele. The last step is to join these new vcf files using the
VCFMerge command:

USAGE:

java -jar NGSEPcore.jar VCFMerge <OPTIONS> <GENOTYPED_VARIANTS_FILE>*

OPTIONS:

	-s FILE	: List of sequence names as they appear in the original
		  reference genome
	-o FILE	: Output VCF file with merged variants and genotype information

This command will write the final vcf file with the genotype calls for each
variant on each sample.

--------------------------------------
Tilling variants individual assignment
--------------------------------------

For tilling experiments, this module takes variants from pools and a pools
descriptor and calls individual variants. It receives a list of VCF files
generated either by the SingleSampleVariantsDetector or the
MultisampleVariantsDetector commands and a pools configuration file and
generates a VCF file with individual genotyping based on the variants
identified within the pools. 

USAGE:

java -jar NGSEPcore.jar TillingPoolsIndividualGenotyper <OPTIONS> <VCF_FILE>*

OPTIONS:

	-r GENOME	: Fasta file with the reference genome.
	-o FILE		: Output file with called variants in VCF format
	-d FILE		: File with the information of individuals assigned to
			  each pool
	-m INT		: Maximum number of pools in which a variant can
			  appear. Default 0 (no filter).
	-b		: Report only biallelic variants.

A pools configuration file must be provided with the option -d. It should be a
text file separated by semicolon and having one row for each individual. The
first entry on each row should be the individual id. The remaining entries
should be the ids of the different pools where the individual was included.
For example, if individual 20 was included in pools with ids 2, 10 and 14, the
line should look like this:

20;2;10;14

The sample ids within input pool VCF files must coincide with the pool
ids present in the pools configuration file. The VCF files can include
information for one or more poools.

----------------------------------------------------------
Obtaining relative allele counts from read alignment files
----------------------------------------------------------

Calculates a distribution of relative allele counts for sites showing base calls
for more than one nucleotide from read alignment files in SAM, BAM or CRAM
format. This analysis is useful to predict the ploidy of a sequenced sample.
Alignments must be sorted by reference coordinates.

USAGE:

java -jar NGSEPcore.jar RelativeAlleleCountsCalculator <OPTIONS>

OPTIONS:

	-i FILE		: Input file with read alignments.
	-o FILE		: Output file with statistics.
	-r GENOME	: Fasta file with the reference genome. Required for
			  CRAM files.
	-m INT		: Minimum read depth Default: 10
	-M INT		: Maximum read depth Default: 1000
	-q INT		: Minimum base quality score (Phred scale) Default: 20
	-frs FILE	: File with repeats (or any kind of genomic regions)
			  that should not be taken into account in the
			  analysis. The format of this file should contain
			  three columns: Sequence name (chromosome), first
			  position in the sequence, and last position in the
			  sequence. Both positions are assumed to be 1-based.
	-srs FILE	: File with genomic regions that should be taken into
			  account in the analysis. The format of this file
			  should contain three columns: Sequence name
			  (chromosome), first position in the sequence, and
			  last position in the sequence. Both positions are
			  assumed to be 1-based.
	-of FILE	: Separate output file with the complete counts
			  information for sites with more than one allele. The
			  file has the following data separated by tab:
			  sequence name, position, read depth, number of
			  different alleles, max depth, and second max depth
	-s		: Consider secondary alignments. By default, only
			  primary alignments are processed.

------------------------------------
Comparing read depth between samples
------------------------------------

This function compares the read depth of two samples to predict regions with
relative copy number variation (CNV) between a sample and a control. It takes
two alignment files and a reference genome, splits the genome into windows, and
for each window compares the read depth between the two samples. It outputs a
text file containing the list of windows of the genome in which the normalized
read depth ratio between the two samples is significantly different from 1.
Alignments can be provided in SAM, BAM or CRAM format and must be sorted by
reference coordinates.

USAGE:

java -jar NGSEPcore.jar CompareRD <OPTIONS>

OPTIONS:

	-i FILE		: Input file with alignments to a reference genome.
	-c FILE		: Input alignments file corresponding to the control
			  (wild type) sample.
	-o FILE		: File with genomic regions in which the two samples
			  have different read depth.
	-r FILE		: Fasta file with the reference genome.
	-w INT		: Window length to be used during the read depth
			  comparison. Default: 100
	-p DOUBLE	: Maximum p-value. Only the windows with a p-value
			  lower than that specified will be reported.
			  Default: 0.001
	-a		: Output an entry for every window in the genome.
	-gc		: Perform GC-correction of the read depth.
	-b		: Perform the Bonferroni correction for multiple 
			  testing.

The output text file contains the following columns:

1. Chromosome
2. Window start
3. Window end
4. Read depth sample 1
5. Read depth sample 2
6. Normalized read depth ratio
7. P-value

----------------------------------------------------------
----------------------------------------------------------
Group 3: Analysis of annotated gene models and transcripts
----------------------------------------------------------
----------------------------------------------------------

-----------------------------------
Evaluating transcriptome assemblies
-----------------------------------

Loads a transcriptome annotation in GFF3 format, logs format errors, provides
statistics on the assembled transcriptome, generates cDNA, CDS and protein
sequences.

USAGE:

java -jar NGSEPcore.jar TranscriptomeAnalyzer <OPTIONS>

OPTIONS:

	-i FILE	: Input GFF3 file with gene annotations. It can be gzip
			  compressed.
	-o FILE	: Prefix of the output files. It can be an absolute path
			  finished by the prefix
	-r FILE	: Fasta file with the reference genome. It can be gzip
			  compressed.

------------------------
Filtering transcriptomes
------------------------

Loads a transcriptome annotation in GFF3 format and generates a filtered file
by CDS length, presence of start and stop codons and intersection with other
regions. Writes to standard output unless the -o option is used to specify an
output file.

USAGE:

java -jar NGSEPcore.jar TranscriptomeFilter <OPTIONS>

OPTIONS:

	-i FILE	: Input GFF3 file with gene annotations. It can be gzip
			  compressed.
	-o FILE	: Output file with filtered genes. See option -f for
			  output format options.
	-r FILE	: Fasta file with the reference genome.
	-f INT		: Output format. 0: GFF3, 1: gene list, 2: gene
			  regions, 3: transcript list, 4: transcript regions.
			  Default: 0
	-c		: Output only complete transcripts (with start and stop
			  codons) in the output file.
	-l INT		: Minimum protein length for coding transcripts in the
			  output file. Default: 0
	-frs FILE	: File with genomic regions in which transcripts should
			  be filtered out. The format of this file should
			  contain three columns: Sequence name (chromosome),
			  first position in the sequence, and last position in
			  the sequence. Both positions are assumed to be 1-based.
	-srs FILE	: File with genomic regions in which transcripts should
			  be selected. The format of this file should contain
			  three columns: Sequence name (chromosome), first
			  position in the sequence, and last position in the
			  sequence. Both positions are assumed to be 1-based.
	-fgid FILE	: File with ids of genes that should be filtered out.
			  The first column should have the gene ids. Other
			  columns are ignored.
	-sgid FILE	: File with ids of genes that should be selected. The
			  first column should have the gene ids. Other columns
			  are ignored.


-----------------
Comparing genomes
-----------------

This module takes a list of assembled genomes in fasta format and their
corresponding transcriptomes in GFF3 format and runs whole genome comparisons.
It calculates orthogroups including orthologs and paralogs. It also identifies
synteny relationships between each pair of genomes. Finally, it calculates gene
presence/absence matrices and classifies gene families as core or accessory.

USAGE:

java -jar NGSEPcore.jar GenomesAligner <OPTIONS> <GENOME1> <TRANSCRIPTOME1> <GENOME2> <TRANSCRIPTOME2>

or

java -jar NGSEPcore.jar GenomesAligner -d <PATHTOFOLDER> -i <INPUTFILE>

OPTIONS:

	-d STRING	: Directory having the input genomes in fasta format
        		  and the genome annotations in gff3 format.
	-i STRING	: Input file with genome identifiers. These identifiers
        		  are used as prefixes for the fasta and gff3 files.
	-o STRING	: Prefix for output files. Default: genomesAlignment
	-k INT		: K-mer length to find orthologs. Default: 6
	-p INT		: Minimum percentage of k-mers to call orthologs. 
			  Default: 11
	-s		: Skip the MCL clustering phase and return unfiltered
			  orthogroups.
	-yh INT	: Minimum number of consistent homology units to call a synteny
			  block. Default: 6
	-yd INT	: Maximum distance (in basepairs) between homology units to
			  include them within the same synteny block. Default: 100000
	-f DOUBLE	: Minimum frequency to classify soft core gene families.
			  Default: 0.9
	-t INT		: Number of threads. Default: 1

The options -d and -i are useful to process large numbers of genome assemblies.
The file referred with the option -i should have one genome identifier for each
line:

Genome1
Genome2
Genome3
...
GenomeN

For each genome, the module will look into the directory referred with the
option -d for one fasta file and one gff3 file having as prefix the genome
identifier.
Possible suffixes for the fasta file include .fa, .fna, .fas and .fasta
and their gzip compressed extensions .fa.gz, .fna.gz, .fas.gz and .fasta.gz
Possible suffixes for the annotation file include .gff and .gff3 and their
gzip compressed extensions .gff.gz and .gff3.gz

		
The output is a series of text files having the ids and physical coordinates of
the paralogs within each genome and the orthologs between the two genomes.
The ortholog files, called <PREFIX>_orthologsG1.tsv and <PREFIX>_orthologsG2.tsv,
have the following format:

1. Id of the gene in the first genome
2. Chromosome of the gene in the first genome
3. Start of the gene in the first genome
4. End of the gene in the first genome
5. Number of paralogs of the gene in the first genome
6. Id of the second genome
7. Id of the ortholog in the second genome
8. Chromosome of the ortholog in the second genome
9. Start of the ortholog in the second genome
10. End of the ortholog in the second genome
11. Alignment type. It can be "L" if the gene has an ortholog in the
   second genome and it makes part of a synteny block. "U" if the gene has
   a unique ortholog but it does not make part of the syntheny block, and
   "M" if the gene has multiple orthologs in the second genome.

The files with the paralogs, called <PREFIX>_paralogsG1.tsv and
<PREFIX>_paralogsG2.tsv, have the same 10 first columns but columns 7 to 10
contain genes within the same genome as genes in column 1 to 4. 

Pangenome files are:
The file <PREFIX>_clusters.txt contains the clusters of homolog genes across
genomes that can be inferred from the pairwise homolog relationships.
The file <PREFIX>_paMatrix.txt contains the Presence/Absence matrix where each
row corresponds to a gene family and each column corresponds to a genome.
The file <PREFIX>_gfFreqs.txt contains the frequency of each gene family within
the genomes and its classification into exact/soft core/accesory genomes.

Finally, the files:

<PREFIX>_linearOrthologView.html
<PREFIX>_circularOrthologView.html and
<PREFIX>_circularParalogView.html

can be loaded in a web browser and provide an interactive view of the alignment
based on the d3 web development technology (https://d3js.org/).

---------------------------------------
Clustering orthologs from CDNA catalogs
---------------------------------------

This module takes cDNA transcriptomes in fasta format, infers orthology
relationships and cluster orthologs. This command is an alternative to the
genomes aligner to infer homologs nad build clusters without the synteny
analysis. Input cDNA files can be generated from annotated genome assemblies
using the command TranscriptomeAnalyzer.

USAGE:

java -jar NGSEPcore.jar CDNACatalogAligner <OPTIONS> <TRANSCRIPTOME>*

OPTIONS:

        -o STRING : Prefix of output files. Default: catalogsAlignment
        -k INT    : K-mer length to find orthologs. Default: 10
        -p INT    : Minimum percentage of k-mers to call orthologs Default: 50
        -s        : Skip the MCL clustering phase and returns unfiltered orthogroups.
        -y INT    : Type of sequences in the input file. 1 for CDNA, 2 for proteins. Default: 1
	-t INT		: Number of threads. Default: 1

This module produces two files as outputs. The first is a text file with
homology relationships. It has three columns separated by tab:

1. Id of the first gene
2. Id of the second gene
3. Homology score

The second file is also a tab delimited file with one line for each identified
cluster. Gene ids within each cluster are separated by tab.

---------------------------------
Identifying transposable elements
---------------------------------

Receives a genome assembly in fasta format and a file with known transposable
elements (TEs) and annotates regions in the assembly with TEs.

USAGE:

java -jar NGSEPcore.jar TransposonsFinder <OPTIONS>

OPTIONS:

	-i FILE	: Input genome to annotate in fasta format. It can be gzip
			  compressed.
	-o FILE	: Output file with annotations of transposable elements.
	-d FILE	: Database of transposable elements to annotate the genome.
	-m INT		: Minimum length (in basepairs) to call a transposable
			  element. Default: 200
	-r INT		: Number of search rounds to identify new TEs from
			  previously identified TEs. Default: 2
	-t INT		: Number of threads. Default: 1


------------------------------------
Masking regions in a genome assembly
------------------------------------

Receives a genome assembly in fasta format and a file of regions
(typically repeats) and masks the regions in the genome, either with lowercase
characters or with Ns.

USAGE:

java -jar NGSEPcore.jar GenomeAssemblyMask <OPTIONS>

OPTIONS:

	-i FILE	: Input genome to mask. It can be gzip compressed.
	-o FILE	: Output file with the masked genome
	-d FILE	: Genomic regions to mask. It must have at least three columns:
			  sequence name (chromosome), 1-based first position and
			  1-based last position.
	-h		: Mask with N characters. The default is to mask with lowercase
			  characters.


--------------------------------------------------------
--------------------------------------------------------
Group 4: Commands for Variants (VCF) downstream analysis
--------------------------------------------------------
--------------------------------------------------------

---------------------------------
Functional annotation of variants
---------------------------------

Generates a VCF file including the functional information related to each
variant. Requires a gff3 file with gene annotations, and the reference genome
in fasta format. Reads from standard input unless the -i option is used to
specify an input file. Writes to standard output unless the -o option is used
to specify an output file.

USAGE:

java -jar NGSEPcore.jar VCFAnnotate <OPTIONS>

OPTIONS:

	-i FILE : Input VCF file with variants to annotate. It can be gzip
		  compressed.
	-r FILE : Fasta file with the reference genome.
	-t FILE : Input GFF3 file with gene annotations.
	-o FILE	: Output VCF file with annotated variants.
	-u INT	: Maximum bp before the start of a transcript to classify a
		  variant as Upstream. Default: 1000
	-d INT	: Maximum bp after the end of a transcript to classify a
		  variant as Downstream. Default: 300
        -sd INT : Initial basepairs of an intron that should be considered as
		  splice donor. Default: 2
        -sa INT : Final basepairs of an intron that should be considered as
		  splice acceptor. Default: 2
        -si INT : Initial or final basepairs of an intron that should be
		  considered as part of the splice region. Default: 10
        -se INT : Initial or final basepairs of an exon that should be
		  considered as part of the splice region. Default: 2

Gene annotations related with the given genome should be provided in standard
GFF3 format. See http://www.sequenceontology.org/gff3.shtml for details.

Annotations in the output VCF file are included using the following custom
fields in the INFO column:

TA (STRING): 	Annotation based on a gene model. Annotation names are terms in
		the sequence ontology database (http://www.sequenceontology.org)
TID (STRING):	Id of the transcript related with the gene annotation in the TA
		field
TGN (STRING): 	Name of the gene related with the annotation in the TA field
TCO (FLOAT):	For variants in coding regions, position in the aminoacid
		sequence where the variant is located. The integer part is the
		1-based position of the mutated codon. The decimal part is the
		codon position.
TACH (String):	Description of the aminoacid change produced by a
		non-synonymous mutation. String encoded as reference aminoacid,
		position and mutated aminoacid

-------------------
Filtering VCF files
-------------------

This module implements different filters on VCF files with genotype
information and generates a VCF file with variants passing the
filtering criteria. The filtering order is as follows: first, it executes the
distance filter (-d option), then the filtering of samples and genotypes (-saf,
-fs, -q and -minRD options). Finally, it recalculates the number of samples
genotyped, the number of alleles called and the MAF to execute the remaining
filters. Since version 2.0.6, the default behavior does not perform any
filtering. 

USAGE:

java -jar NGSEPcore.jar VCFFilter <OPTIONS>

OPTIONS:

	-i FILE		: Input file in VCF format. It can be gzip compressed.
	-o FILE		: Output file in VCF format.
	-frs FILE	: File with genomic regions in which variants should be
			  filtered out. The format of this file should contain
			  at least three columns: Sequence name (chromosome),
			  first position in the sequence, and last position in
			  the sequence. Both positions are assumed to be
			  1-based.
	-srs FILE	: File with genomic regions in which variants should be
			  selected. The format of this file should contain at
			  least three columns: Sequence name (chromosome),
			  first position in the sequence, and last position in
			  the sequence. Both positions are assumed to be
			  1-based.
	-d INT		: Minimum distance between variants. Default: No filter
	-q INT		: Minimum genotyping quality score (GQ format field in
			  the VCF). Genotype calls with lower GQ become
			  undecided. Default: 0
        -minRD INT	: Minimum read depth (as reported in the DP genotype
			  format field) to keep a genotype call. Genotype calls
			  with less reads become undecided. Default: 0
	-s		: Keep only biallelic SNVs.
        -fi		: Filter out sites in which only one allele is
			  observed.
        -fir		: Filter sites in which only the reference allele is
			  observed.
        -fia		: Filter out sites in which only one alternative allele
			  is observed.
	-m INT		: Minimum number of samples genotyped to keep the
			  variant. Default: 0 (No filter)
	-minMAF FLOAT	: Minimum minor allele frequency over the samples in
			  the VCF file. Default: 0 (No filter)
	-maxMAF FLOAT	: Maximum minor allele frequency over the samples in
			  the VCF file. Default: 0.5 (No filter)
	-minOH FLOAT	: Minimum observed heterozygosity over the samples in
			  the VCF file. Default: 0 (No filter)
	-maxOH FLOAT	: Maximum observed heterozygosity over the samples in
			  the VCF file. Default: 1 (No filter)
	-g FILE		: File with the reference genome to calculate the
			  GC-Content of the region surrounding the variant.
	-minGC FLOAT	: Minimum percentage of GC of the 100bp region
			  surrounding the variant. Default: 40.0
	-maxGC FLOAT	: Maximum percentage of GC of the 100bp region
			  surrounding the variant. Default: 65.0
	-maxCNVs INT	: Maximum number of samples with copy number variation
			  in the region where the variant is located.
			  Default: No filter
	-gene STRING	: Id of the gene or the transcript related with the
			  variant.
	-a STRING	: Types of functional annotations related to the
			  variants. 
	-saf FILE	: File with the ids of the samples to be selected (or
			  removed, see -fs option). The file should have one
			  line per sample, being the first column the sample
			  id. Other columns in the file are ignored.
	-fs		: Flag to remove the samples provided with the -saf
			  option instead of selecting them. 

Names of functional annotations to use with the option -a should correspond to
standard sequence ontology terms (http://www.sequenceontology.org). More than
one annotation can be set as a comma-separated list. Common terms include
missense_variant,synonymous_variant,frameshift_variant,start_lost,stop_gained
among others.

----------------------------------
Convert VCF files to other formats
----------------------------------

Convert genotype calls in VCF format to other formats commonly used to perform
different kinds of downstream analysis.

USAGE:

java -jar NGSEPcore.jar VCFConverter <OPTIONS> <INPUT_FILE> <OUTPUT_PREFIX>

OPTIONS:

	-i FILE		: Input VCF file with variants and genotype data.
			  It can be gzip compressed.
	-o FILE		: Prefix of the output files.
	-darwin		: Generates the input files for DarWin
	-eigensoft	: Generates the input files for Eigensoft
	-emma		: Generates the input files for Emma.
	-fasta		: Generates a virtual multiple sequence alignment in
			  fasta format. It could be used to build distance
			  based dendograms.
        -fineStructure	: Generates the input files for FineStructure. The
			  option -s is required for this format.
	-flapjack	: Generates the input files for Flapjack.
	-genepop	: Generates the input format for GenePop.
	-GWASPoly	: Generates the input file for GWASPoly.
	-haploview	: Generates the input files for Haploview.
	-hapmap		: Generates the Hapmap format, which can be used in
			  programs such as Tassel.
	-joinMap	: Generates the input file to build genetic maps with
			  JoinMap. The options -p1 and -p2 are required for
			  this format.
	-matrix		: Generates a simple ACGT format which can be imported
			  to excel.
	-phase		: Generates the input file for PHASE. The option -s is
			  required for this format.
	-plink		: Generates the input files for Plink.
	-powerMarker	: Generates the input files for Powermarker.
	-rrBLUP		: Generates the input files for rrBLUP.
	-spagedi	: Generates the input files for Spagedi.
	-structure	: Generates the input format for structure.
	-treeMix	: Generates the input files for TreeMix. The option -p
			  is required for this format.
	-s STRING	: Name of the sequence (chromosome) for conversion to
			  PHASE.
	-p FILE		: File with population assignments for the samples.
			  This should be a two column text file with the sample
			  ids in the first column and the ids of the
			  populations in the second column. Required for
			  conversion to TreeMix.
	-p1 STRING	: Id of the first parent for conversion to JoinMap
	-p2 STRING	: Id of the second parent for conversion to JoinMap

WARNING: FASTA convertion does not use IUPAC codes, heterozygous SNPs are 
changed to N.
WARNING 2: Plink is only designed for humans, therefore it will only work for
22 sequences (chromosomes). If a sample exceeds this number, it is convenient 
to reduce the number of chromosomes and to remove all scaffolds.
WARNING 3: To generate dendograms in Tassel, it is better to use the HapMap 
format.

-------------------
Comparing VCF files
-------------------
Compares the genotype calls included in two different VCF files. Calculates
the number and percentage of homozygous and heterozygous differences between
every pair of samples. It requires the FASTA file with the reference genome
used to build the VCF files. If only the first input file is provided, this
module provides an internal comparison of the samples within the input file.
Writes to standard output unless the -o option is used to specify an output
file.

USAGE:

java -jar NGSEPcore.jar VCFComparator <OPTIONS>

OPTIONS: 

	-i FILE		: Input file in VCF format. It can be gzip compressed.
	-i2 FILE	: Input file in VCF format to compare with the first
			  file. It can be gzip compressed.
	-r FILE		: Fasta file with the reference genome.
	-o FILE		: Output file with the results of the comparison.
	-g DOUBLE	: Minimum percentage (0-100) of variants genotyped in
			  both samples. Default: 50.
	-d DOUBLE	: Maximum percentage (0-100) of differences between the
			  pair of samples. Default: 5.

Default values of optional parameters are set to facilitate the detection of
duplicated (or very similar) samples. To report the complete set of sample
pairs, use -g 0 -d 100. The output is a tab-delimited report with the following
fields:

1. Id sample VCF 1
2. Id sample VCF 2
3. Number of variants genotyped in sample 1
4. Number of variants genotyped in sample 2
5. Number of variants genotyped in both samples
6. Number of heterozygous differences
7. Percentage of heterozygous differences (sixth field / fifth field)
8. Number of homozygous differences
9. Percentage of homozygous differences (eighth field / fifth field)
10. Number of total differences
11. Percentage of total differences (tenth field / fifth field) 

------------------------------
Calculating summary statistics
------------------------------

Generate a report with the variants included in a VCF file for different
categories. It is specially useful when a complete population is being
processed and merged into a single annotated file. Reads from standard input
unless the -i option is used to specify an input file. Writes to standard
output unless the -o option is used to specify an output file.

USAGE:

java -jar NGSEPcore.jar VCFSummaryStats <OPTIONS>

OPTIONS:

	-i FILE	: Input file in VCF format. It can be gzip compressed.
	-o FILE	: Output file with statistics.
	-m INT	: Minimum number of samples genotyped to accurately calculate
		  the minor allele frequency. Default: 20

-----------------------------------------
Calculating diversity statistics per site
-----------------------------------------

Calculates basic diversity statistics for each variant in a VCF file. Reads
from standard input unless the -i option is used to specify an input file.
Writes to standard output unless the -o option is used to specify an output file.

USAGE:

java -jar NGSEPcore.jar VCFDiversityStats <OPTIONS>

OPTIONS:

	-i FILE	: Input file in VCF format. It can be gzip compressed.
	-o FILE	: Output file with statistics.
	-p FILE	: File with population assignments for the samples. A two
		  column text file with the sample ids in the first column and
		  the ids of the populations in the second column.

The output file contains the coordinates of each variant plus the following
statistics separated by semicolon:

1. Number of samples genotyped
2. Expected heterozygosity (under HWE)
3. Observed heterozygosity
4. F-statistic (1-OH/EH)
5. Minor allele frequency (MAF)
6. Chi-square value of departure from HWE
7. Uncorrected p-value of the Chi-square test for departure from HWE

If a file with population assignments is provided, this module will output one
column of statistics for the whole group and one column for each population.


----------------------------
Calculating variants density
----------------------------

Calculates the number of variants within a VCF file in non-overlapping windows
across the genome. Writes a text delimited file with four columns: sequence,
window first, window last and number of variants. Reads from standard input
unless the -i option is used to specify an input file. Writes to standard
output unless the -o option is used to specify an output file.

USAGE:

java -jar NGSEPcore.jar VCFVariantDensityCalculator <OPTIONS>

OPTIONS:

	-i FILE	: Input file in VCF format. It can be gzip compressed.
	-o FILE	: Output file with statistics.
	-r FILE	: Fasta file with the reference genome.
	-w INT  : Length of the window. Default: 100000

-------------------------------------------------------
Calculation of genetic distance matrices from VCF files
-------------------------------------------------------

Generates a distance matrix from a variants file in VCF format. The matrix is
calculated using the basic IBS (Identity by state) algorithm. However, four
options to infer the genotype call information are implemented. In particular,
users can choose predicted allele dosages of CNVs or direct estimations of
allele dosage per site per individual based on relative allele-specific read
counts. The latter option is useful to improve distance estimations in
polyploids. It writes the matrix of genetic distances in a generic format.
Reads from standard input unless the -i option is used to specify an input
file. Writes to standard output unless the -o option is used to specify an
output file.

USAGE:

java -jar NGSEPcore.jar VCFDistanceMatrixCalculator <OPTIONS>

OPTIONS:

	-i FILE	: Input file in VCF format. It can be gzip compressed.
	-o FILE	: Output file with the distance matrix. See -f for options on
		  the output format.
	-s INT	: Source of information in the VCF file to calculate distances.
		  0 for simple genotype calls (GT format field), 1 for allele
		  copy number (ACN format field), 2 for total copy number
		  (total of ACN format field), and 3 for raw allele depth (ADP
		  or BSDP format fields). Default: 0
	-f INT	: Matrix output format, 0 is full matrix, 1 lower-left matrix
		  and 2 is upper right matrix. Default: 0
	-p INT	: Default ploidy of the samples. Used if the distance source
		  (-s option) is the raw allele depths to recalculate allele
		  dosage based on these counts. Default: 2

--------------------------------------------------------
Building dendograms using the Neighbor-Joining algorithm
--------------------------------------------------------

Given a distance matrix file, this command builds a dendogram for graphical
display of genetic distances using the Neighbor Joining algorithm. The distance
matrix can be provided as an upper, lower or full matrix. The dendogram is
written in Newick format. Reads from standard input unless the -i option is
used to specify an input file. Writes to standard output unless the -o option
is used to specify an output file.

USAGE:

java -jar NGSEPcore.jar NeighborJoining <OPTIONS>

OPTIONS:

	-i FILE	: Input file with a distance matrix.
	-o FILE	: Output file with the dendrogam in Newick format.


-------------------------------------
Calculating allele sharing statistics
-------------------------------------

Calculates allele sharing diversity statistics, either through windows across
the genome or through the genes catalog of the species. This program calculates
the pairwise differences between every pair of samples in the VCF file and uses
that information to calculate diversity statistics such as the average number
of pairwise differences per Kbp, Fst and Tajima D. This functionality should
only be applied to VCFs containing populations of inbred samples. Each group
can either be one or more than populations wthin the populations file. Multiple
population names within one group should be separated by comma (without white
spaces). Reads from standard input unless the -i option is used to specify an
input file. Writes to standard output unless the -o option is used to specify
an output file.

USAGE:

java -jar NGSEPcore.jar VCFAlleleSharingStats <OPTIONS>

OPTIONS:


	-i FILE		: Input file in VCF format. It can be gzip compressed.
	-p FILE		: File with population assignments for the samples.
			  A two column text file with the sample ids in the
			  first column and the ids of the populations in the
			  second column.
	-o FILE		: Output file with statistics.
	-g1 STRING	: Comma-separated list of populations that should be
			  considered as group 1.
	-g2 STRING	: Comma-separated list of populations that should be
			  considered as group 2.
	-t FILE		: GFF3 file with the transcriptome. If this file is
			  provided, statistics will be provided by gene and not
			  by window.
	-n		: If set, introns will be included in the calculation
			  of pairwise differences. Only useful if the -t option
			  is set.
	-w INT		: Length of each genomic window to calculate pairwise
			  differences between samples. Default: 100000
	-s INT		: Step between windows to calculate pairwise
			  differences between samples. Default: 10000


The populations file is a tab-delimited text file with two columns:
sample id and population id. Writes a tab-delimited report with the following
fields:

1. Chromosome
2. Window start
3. Length of the region in Kbp
4. Total number of variants within the window
5. Segregating sites within the group 1
6. Segregating sites within the group 2
7. Segregating sites within the two groups
8. Diversity measured as average number of pairwise differences per Kbp within
   group 1
9. Diversity within group 2
10. Diversity between the two groups
11. Diversity within the two groups
12. Diversity across all samples in the two groups
13. Diversity across all samples in the file
14. Fst between the two groups measured as the difference between diversity
    between and within groups divided by the diversity between groups.
15. Tajima D within the group 1
16. Tajima D within the group 2 

If the -t option is used, the first two columns are replaced by the transcript
id and gene id respectively.

-------------------
Genotype imputation
-------------------

This module allows imputation of missing genotypes from unphased multilocus
SNP genotype data in a VCF. The current version is a reimplementation of the
Hidden Markov Model (HMM) implemented in the package fastPHASE
(http://stephenslab.uchicago.edu/software.html). This implementation allows to
process VCF files and produces its output also as a VCF. However, only
biallelic SNPs are imputed and included in the output VCF file. The current
version supports imputation of either highly homozygous or heterozygous
populations. Parental lines can be provided for both types of populations using
the -p option. The options -ip and -is tell the model that either the parental
accessions (-ip) or the entire population (-is) are inbred samples with low
heterozygosity. In the latter mode, the model will only produce homozygous
genotype calls.

USAGE:

java -jar NGSEPcore.jar ImputeVCF <OPTIONS>

OPTIONS: 

	-i FILE		: Input file in VCF format. It can be gzip compressed.
	-o FILE		: Prefix of the output files.
	-p STRING	: Comma-separated list of sample ids of the parents of
			  the breeding population. This should only be used for
			  bi-parental or multi-parental breeding populations.
	-k INT		: Maximum number of groups in which local haplotypes
			  will be clustered. See (PMID:16532393) for details
			  of the HMM implemented in the fastPHASE algorithm.
			  For bi-parental or multi-parental breeding
			  populations please set explicitly the number of
			  parents of the population even if the list of parents
			  is provided with the -p option. This allows to take
			  into account cases of populations in which some of
			  the parents are missing. Default: 8
	-w INT		: Window size as number of variants within each window.
			  Default: 5000
        -v INT		: Overlap as number of variants shared between neighbor
			  windows. Default: 50
	-c FLOAT	: Estimated average number of centiMorgans per Kbp on
			  euchromatic regions of the genome. This value is used
			  by the model to estimate initial transitions between
			  the states of the HMM. Typical values of this
			  parameter are 0.001 for human populations, 0.004 for
			  rice and 0.35 for yeast populations. Default: 0.001
	-t		: If set, transition probabilities in the HMM will NOT
			  be updated during the Baum-Welch training of the HMM. 
			  Not recommended unless the -c option is set to a
			  value allowing a reasonable initial estimation of the
			  transition probabilities.
        -ip		: Specifies that parents of the population are inbred.
        -is		: Specifies that the samples to impute are inbred.


This module outputs two files, the first is a VCF file including the imputed
genotypes for the datapoints having an undecided genotype call in the input
file. The second outputs for each SNP and each sample the index of the parent
that most likely originated the observed haplotype of the individual.

--------------------------------
Finding haplotype introgressions
--------------------------------

Runs a window-based analysis to identify the most common haplotype within each
of the populations described in the given populations file and then identify
common haplotypes of one population introgressed in samples of a different
population. Although it can be run on any VCF file, it is particularly designed
to work with populations of inbred samples. Reads from standard input unless the
-i option is used to specify an input file.

USAGE:

java -jar NGSEPcore.jar VCFIntrogressionAnalysis <OPTIONS>

OPTIONS:

	-i FILE		: Input file in VCF format. It can be gzip compressed.
	-p FILE		: File with population assignments for the samples.
			  A two column text file with the sample ids in the
			  first column and the ids of the populations in the
			  second column.
	-o FILE		: Prefix of the output files.
	-g FLOAT	: Minimum percentage of samples genotyped within a
			  population to identify the most common allele.
			  Default: 80
	-d FLOAT	: Minimum difference between reference allele
			  frequencies of at least two populations to consider a
			  variant discriminative. Default: 0.6
	-m FLOAT	: Maximum minor allele frequency within a population to
			  consider the major allele of a variant as
			  representative allele for such population.
			  Default: 0.4
	-w INT		: Window size as number of variants within each window.
			  Default: 50
        -v INT		: Overlap as number of variants shared between neighbor
			  windows. Default: 0
	-a INT		: Score given of a match between homozygous genotypes
			  comparing haplotypes Default: 1
	-t INT		: Score given of a mismatch between homozygous
			  genotypes comparing haplotypes Default: -1
	-s INT		: Minimum score to match an individual haplotype with a
			  population-derived haplotype Default: 30
	-c		: Outputs a VCF file with the biallelic variants that
			  showed segregation between at least one pair of
			  groups and hence were selected for the analysis.
        -u		: If set, reports introgression events for unassigned
			  haplotypes according to the minimum score defined by
			  the options -a -i and -s

By default, this function outputs three files:
1. <OUT_PREFIX>_assignments.txt: Table with population assignments for each
   genomic region (in rows) and each sample (in columns). If a sample does not
   have enough variants genotyped within a region, an "M" (missing) will appear
   in the corresponding population assignment. If a sample has a haplotype that
   does not match any population haplotype according to the minimum score
   defined by the options -a -i and -s, a "U" (Unassigned) will appear in the
   population assignment. If the difference between scores of the best and the
   second best population assignments is less than 10, the two populations and
   their scores will be reported.
2. <OUT_PREFIX>_introgressions.txt: Introgressions identified by the analysis.
   Includes the genomic region, the sample id, the sample population, the
   population where the haplotype is most common (e.g. introgression origin),
   the total number of variants analyzed within the region, the number of
   variants genotyped for the sample, and the score obtained for each
   population. This report aggregates in a single event assignments over
   consecutive regions reported in the assignments file. If the -u option is
   set, it also outputs events in which the haplotype does not match any
   known population haplotype, which could indicate an introgression from a
   population not included in the dataset.
3. <OUT_PREFIX>_assignmentStats.txt: Summary report with the number of
   variants that segregate between each pair of populations. For each
   population, the report also contains the number of variants that do not meet
   the minimum percentage of samples genotyped (according to the -g option) and
   the number variants with high MAF (heterozygous according to the -m option).
   Finally, it also reports for each sample the number of regions assigned to
   each population, and the number of regions non genotyped, unassigned and
   assigned to more than one population.

----------------------------------------
Mapping de-novo GBS variants to a genome
----------------------------------------

Given a VCF file with coordinates relative to a set of consensus sequences, and
a reference genome, aligns the consensus sequences and provides a new VCF file
with coordinates relative to the reference genome. This command is useful to
quickly translate variants identified with the DeNovoGBS command to an
assembled genome.

USAGE:

java -jar NGSEPcore.jar VCFRelativeCoordinatesTranslator <OPTIONS>

OPTIONS:

	-i FILE		: Input VCF file having variants with coordinates
			  relative to a given set of consensus sequences.
	-r GENOME	: Fasta file with the reference genome.
	-o FILE		: Prefix of the output files.
	-b FILE		: BAM file with alignments of the consensus sequences
			  to the given reference genome.
	-c FILE		: Fasta file with consensus sequences. Only used if the
			  -b option is not used.
	-d FILE		: FM-index file of the reference genome calculated with
			  the command GenomeIndexer. Only used if the consensus
			  sequences are provided in FASTA format (See option -c).

This command produces as main output the file <OUT_PREFIX>.vcf with the
translated variants. It also produces a file called <OUT_PREFIX>.info providing
statistics on number and percentage of aligned consensus sequences and
translated variants. Finally, if consensus sequences are provided in FASTA
format, it produces the file <OUT_PREFIX>_alns.bam with the alignments of the
given consensus sequences.

---------------------------------
---------------------------------
Group 5: Simulation and benchmark
---------------------------------
---------------------------------

----------------------------------------------
Simulating individuals from a reference genome
----------------------------------------------

This simulator takes a (haploid) genome assembly and simulates a single
individual including homozygous and heterozygous mutations (SNPs, indels and
mutated STRs) relative to the input assembly. It produces two files, a fasta
file with the simulated genome, and a phased VCF file with the simulated
variants.

USAGE:

java -jar NGSEPcore.jar SingleIndividualSimulator <OPTIONS>

OPTIONS:

	-i FILE		: Fasta file with the genome to simulate an individual.
	-o FILE		: Prefix of the output files.
	-s DOUBLE	: Proportion of reference basepairs with simulated SNV
			  events. Default: 0.001
	-n DOUBLE	: Proportion of reference basepairs with simulated
			  indel events. Default: 1.0E-4
	-f DOUBLE	: Fraction of input STRs for which a mutation will be
			  simulated. Default: 0.1
	-t FILE		: Path to a text file describing the known STRs in the
			  given genome.
	-u INT		: Zero-based index in the STR file where the unit
			  sequence is located. Default: 14
	-d STRING	: ID of the simulated sample. Appears in the VCF header
			  and as part of the name of the sequences in the
			  simulated genome. Default: Simulated
	-p INT		: Ploidy of the simulated sample. Default: 2

The file with known STRs should have at least four columns:

1. Sequence name (chromosome)
2. First basepair of the STR (1-based inclusive)
3. Last basepair of the STR (1-based inclusive)
4. STR unit sequence

The option -u allows to indicate the actual column where the unit sequence is
located. At this moment, the default value corresponds to the column where this
sequence is located in the (proceesed) output of tandem repeats finder (TRF).

----------------
Simulating reads
----------------

Generates single reads randomly distributed from a given reference genome.

USAGE:

java -jar NGSEPcore.jar SingleReadsSimulator <OPTIONS>

OPTIONS:

	-i FILE		: Fasta file with the genome to simulate reads.
	-o FILE		: Gzip compressed output file with simulated reads.
			  See option -f for options on the file format.
	-n INT		: Number of reads. Default: 30000
	-u INT		: Average read length. Default: 20000
	-s INT		: Standard deviation read length. Default: 10000
        -m INT		: Minimum read length. Default: 50
	-e DOUBLE	: Substitution error rate. Default: 0.01
	-d DOUBLE	: Indel error rate. Default: 0.01
	-f INT		: Output format. 0 for fastq, 1 for fasta Default: 0

------------------------------
Simulating TILLING experiments
------------------------------

Simulates a mutagenized population from selected regions on the given reference
genome. Distributes samples in pools and simulates reads from amplicons
assigned to each pool.

USAGE:

java -jar NGSEPcore.jar TillingPopulationSimulator <OPTIONS>

OPTIONS:

	-i FILE		: File with the description of the regions that will be
			  used as amplicons for the simulation.
	-g GENOME	: Fasta file with the genome to simulate reads.
	-o FILE		: Prefix of the output files
	-d INT		: Number of individuals to simulate. It should be less
			  or equal than the product of the three dimensions of
			  the pool design (parameters d1, d2 and d3).
			  Default: 288
	-n INT		: Number of fragments to sequence for each pool.
			  Default: 50000
	-m INT		: Number of mutations to generate. Default: 300
	-u INT		: Read length. Default: 200
	-e1 DOUBLE	: Minimum substitution error rate (at the 5' end).
			  Default: 0.001
	-e2 DOUBLE	: Maximum substitution error rate (at the 3' end).
			  Default: 0.01
	-d1 INT		: First dimension of the pool design. Default: 6
	-d2 INT		: Second dimension of the pool design. Default: 8
	-d3 INT		: Third dimension of the pool design. Default: 6

The following files are generated with the given prefix:

- A VCF file with the simulated mutations for each individual.
- A text file separated by semicolon with the pools assignment to each individual.
  This file can be loaded in the TilligPoolsIndividualGenotyper.
- Two fastq files for each pool with the simulated reads.

--------------------------
Benchmarking variant calls
--------------------------

Takes a VCF file with genotype information from one sample, the reference
genome used to build the VCF and a phased VCF file with gold standard calls and
calculates quality statistics comparing gold-standard with test calls. Writes
to standard output unless the -o option is used to specify an output file.

USAGE:

java -jar NGSEPcore.jar VCFGoldStandardComparator <OPTIONS>

OPTIONS:

	-i FILE	: Input test file in VCF format. It can be gzip compressed.
	-g FILE	: Gold standard file in VCF format. It can be gzip compressed.
	-r FILE	: Fasta file with the reference genome.
	-o FILE	: Output file with statistics.
	-c FILE	: File with coordinates of complex regions (such as STRs).
	-f FILE	: File with coordinates of regions in which the gold standard
		  can be trusted.
	-e	: Indicates that the gold standard VCF is genomic, which means
		  that confidence regions can be extracted from annotated
		  regions with homozygous reference genotypes.

The output is a tab delimited file with the following fields:
1. Minimum genotype quality score (GQ field)
2. Homozygous reference calls in homozygous reference regions
3. Heterozygous calls in homozygous reference regions
4. Homozygous alternative calls in homozygous reference regions
5. Homozygous reference calls in heterozygous regions
6. Heterozygous calls in heterozygous regions
7. Homozygous alternative calls in heterozygous regions
8. Homozygous reference calls in homozygous alternative regions
9. Heterozygous calls in homozygous alternative regions
10. Homozygous alternative calls in homozygous alternative regions
11. Non matched gold standard homozygous reference calls
12. Non matched gold standard heterozygous calls
13. Non matched gold standard homozygous alternative calls
14. Non matched test homozygous reference calls
15. Non matched test heterozygous calls
16. Non matched test homozygous alternative calls
17. Total gold standard homozygous reference calls
18. Total gold standard heterozygous calls
19. Total gold standard homozygous alternative calls
20. Total test homozygous reference calls
21. Total test heterozygous calls
22. Total test homozygous alternative calls
23. Recall heterozygous calls
24. False discoveries heterozygous calls
25. FPPM heterozygous calls
26. FDR heterozygous calls
27. Precision heterozygous calls
28. F1 heterozygous calls
29. Recall homozygous calls
30. False discoveries homozygous calls
31. FPPM homozygous calls
32. FDR homozygous calls
33. Precision homozygous calls
34. F1 homozygous calls

The current output also includes distributions of gold standard variants per
cluster, heterozygous test variants per cluster and genome span per cluster

------------------------------
Citing and supporting packages
------------------------------

The latest algorithms implemented in NGSEP 3 to improve accuracy for variants
detection and genotyping were recently published in bioinformatics:

Tello D, Gil J, Loaiza CD, Riascos JJ, Cardozo N, and Duitama J. (2019)
NGSEP3: accurate variant calling across species and sequencing protocols.
Bioinformatics 35(22): 4716–4723.
http://doi.org/10.1093/bioinformatics/btz275

Further details on the pipeline built for variants detection on
Genotype-By-Sequencing (GBS) data can be found at BMC Genomics:

Perea C, Hoz JFDL, Cruz DF, Lobaton JD, Izquierdo P, Quintero JC, Raatz B and Duitama J. (2016).
Bioinformatic analysis of genotype by sequencing (GBS) data with NGSEP.
BMC Genomics, 17:498.
http://doi.org/10.1186/s12864-016-2827-7

The first manuscript with the initial description of the main modules of NGSEP
is available at Nucleic Acids research:

Duitama J, Quintero JC, Cruz DF, Quintero C, Hubmann G, Foulquie-Moreno MR, Verstrepen KJ, Thevelein JM, and Tohme J. (2014). 
An integrated framework for discovery and genotyping of genomic variants from high-throughput sequencing experiments. 
Nucleic Acids Research. 42(6): e44. 
http://doi.org/10.1093/nar/gkt1381


Details of algorithms implemented in NGSEP for different functionalities can be
found in the following publications:

CNV detection (Read depth analysis):
Abyzov, A., Urban, A. E., Snyder, M., and Gerstein, M. (2011). 
CNVnator: an approach to discover, genotype, and characterize typical and atypical CNVs from family and population genome sequencing.
Genome research, 21(6), 974–84. 
http://doi.org/10.1101/gr.114876.110

Yoon S, Xuan Z, Makarov V, Ye K, Sebat J. (2009).
Sensitive and accurate detection of copy number variants using read depth of coverage. 
Genome Research 19(9):1586-1592.
http://doi.org/10.1101/gr.092981.109

Genotype imputation:
Scheet, P and Stephens, M. (2006).
A Fast and Flexible Statistical Model for Large-Scale Population Genotype Data: Applications to Inferring Missing Genotypes and Haplotypic Phase.
American Journal of Human Genetics 78: 629-644.
http://doi.org/10.1086/502802

Read Depth comparison:
Xie C and Tammi MT. (2009).
CNV-seq, a new method to detect copy number variation using high-throughput sequencing.
BMC Bioinformatics 10:80.
http://doi.org/10.1186/1471-2105-10-80

Since version 2.1.2, we implemented a new model to integrate paired-end
and split-read analysis for detection of large indels. A recent benchmark
experiment of this algorithm against other software tools using data from the
3000 rice genomes project is available at Genome Research:

Fuentes RR, Chebotarov D, Duitama J, Smith S, De la Hoz JF, Mohiyuddin M, et al. (2019).
Structural variants in 3000 rice genomes.
Genome Research 29: 870-880.
http://doi.org/10.1101/gr.241240.118

NGSEP is also supported by the following open source software packages:

Bowtie2: http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
Picard: http://picard.sourceforge.net/
Jsci: http://jsci.sourceforge.net/
XChart: http://xeiam.com/xchart/
Trimmomatic: http://www.usadellab.org/cms/?page=trimmomatic. We borrowed one
class from Trimmomatic 0.35 to allow correct reading of gzip files
