NGSEP - Next Generation Sequencing Experience Platform
Version 3.3.1 (08-03-2019)
===========================================================================

NGSEP provides an object model to enable different kinds of
analysis of DNA high throughput sequencing (HTS) data. The most important
use of NGSEP is the construction and downstream analysis of large datasets of
genomic variation. NGSEP performs accurate detection and genotyping of Single
Nucleotide Variants (SNVs), small and large indels, short tandem repeats (STRs),
inversions, and Copy Number Variants (CNVs). NGSEP also provides utilities for
downstream analysis of variation in VCF files, including functional annotation
of variants, filtering, format conversion, comparison, clustering, imputation,
introgression analysis and different kinds of statistics. Other functionalities
include calculations of k-mer distributions from fasta or fastq files,
demultiplexing of barcoded sequencing reads, and comparative analysis of read
depth distributions.

--------------------
Building NGSEP
--------------------

NGSEP has been compiled and run successfully on the standard jdk version
1.8.0. To build the distribution library NGSEPcore.jar on a unix based
command line environment run the following commands in the directory where
NGSEPcore_3.3.1.tar.gz is located:

tar -xzvf NGSEPcore_3.3.1.tar.gz
cd NGSEPcore_3.3.1
make all

Note: Usage fields below do not include the version number. To remove the
version number, users can either copy the executable jar file:

cp NGSEPcore_3.3.1.jar NGSEPcore.jar

or just make a symbolic link:

ln -s NGSEPcore_3.3.1.jar NGSEPcore.jar

---------------
Asking for help
---------------

It is possible to obtain usage information for each module by typing:

java -jar NGSEPcore.jar <MODULE> --help

General information and the list of modules can be obtained by typing:

java -jar NGSEPcore.jar [ --help | --version | --citing ]

--------------------------------------
Calling variants over multiple samples
--------------------------------------

This modules allows to call variants over a group of samples separated by files
or read group tags. This is now the recommended method to perform variants
detection on genotype-by-sequencing (GBS), RAD sequencing, whole exome
sequencing (WES), RNA-seq and low coverage (less than 10x) whole genome
sequencing (WGS) data. Although it can also be used on high coverage WGS data,
the classic sample-by-sample analysis (commands FindVariants, MergeVariants and
MergeVCF) is still recommended to identify structural variants. This module
requires one or more read alignment files in BAM format and the reference genome
that was used to produce the alignments.

USAGE:

java -jar NGSEPcore.jar MultisampleVariantsDetector <OPTIONS> <BAM_FILE>*

OPTIONS:

	-r GENOME		: Fasta file with the reference genome.
	-o FILE			: Output file. Default: variants.vcf
	-h DOUBLE		: Heterozygosity rate. Default: 0.001
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
				  this filter is executed independently for each
				  read group. For GBS or RAD sequencing data use
				  a large value such as 100. Default: 5
	-p			: Process non unique primary alignments in the
				  pileup process. The default behavior is to
				  process alignments that are unique
				  (see option -minMQ).
	-s			: Consider secondary alignments in the pileup
				  process. Non-unique primary alignments will
				  also be considered in this mode.
	-minQuality INT         : Minimum variant quality. In this command, this
				  filter applies to the QUAL column of the VCF,
				  which is calculated for each variant as the
				  maximum of the genotype qualities of samples
				  with non-homozygous reference genotype calls.
				  See the command FilterVCF to apply filters of
				  quality and read depth on individual genotype
				  calls. Default: 40
	-maxBaseQS INT          : Maximum value allowed for a base quality score.
				  Larger values will be equalized to this value. Default: 100
	-ignore5 INT            : Ignore this many base pairs from the 5' end of
				  the reads. Default: 0
	-ignore3 INT            : Ignore this many base pairs from the 3' end of
				  the reads. Default: 0
	-knownSTRs FILE         : File with known short tandem repeats (STRs).
				  This is a text file with at least three
				  columns: chromosome, first position and last
				  position. Positions should be 1-based and
				  inclusive.
	-knownVariants FILE     : VCF file with variants to be genotyped. Only
				  these variants will appear in the output VCF.
	-embeddedSNVs           : Flag to call SNVs within STRs. By default,
				  STRs are treated as a single locus and hence
				  no SNV will be called within an STR.
	-minMQ INT              : Minimum mapping quality to call an alignment
				  unique Default: 20
	-ploidy INT             : Default ploidy of the samples. Default: 2
	-psp                    : Print id and ploidy of the sample in the VCF
				  header. The header generated with this option
				  is not a standard VCF header. However, it helps
				  NGSEP to keep track of the ploidy of the samples
				  through downstream analyses

Alignments should be provided in BAM v1 format
(see http://samtools.github.io/hts-specs/SAMv1.pdf for details).
The reference genome should be provided in fasta format. It is assumed that
the sequence names in the alignments file correspond with the sequence names in
this reference assembly. For this module, BAM files must include RG headers
including the ID of each read group and the corresponding sample. A typical read
group header looks as follows

@RG	ID:<ReadGroupId>	SM:<SampleId>	PL:<Platform>

Each alignment must include an RG tag indicating the id of the read group of
the aligned read. See the BAM format specification for details.
Check the documentation of your read aligner to make sure that BAM files contain
read group headers. This module uses read group headers to distribute the reads
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
MAF (DOUBLE)	: Minor allele frequency. Calculated by the MergeVCF and the
		  FilterVCF commands
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

java -jar NGSEPcore.jar FindVariants -maxAlnsPerStartPos 100 -r <REFERENCE> <BAM_FILE>*

WARNING 2: Unlike the behavior of the classical individual analysis per sample,
in this command the filter executed using the minQuality option applies to the
QUAL field of the VCF format, which corresponds to the probability of existence
of each variant (regardless of the individual genotype calls). In this module
the QUAL probability is calculated for each variant as the maximum of the
genotype qualities of samples with non-homozygous reference genotype calls. The
rationale for this calculation is that one variant should be real if it is
confidently called in at least one sample. Individual genotype calls are not
filtered by default and hence they could include some false positives. Please see
the FilterVCF command to perform custom filtering of genotype calls, either by
genotype quality (GQ format field) or by read depth (BSDP and ADP format fields).
Default values of other parameters are also set to maximize sensitivity. For
conservative variant detection including control for errors in base quality
scores and PCR amplification artifacts use:

java -jar NGSEPcore.jar FindVariants -maxAlnsPerStartPos 2 -maxBaseQS 30 -r <REFERENCE> <BAM_FILE>*

If the error rate towards the three prime end increases over 2% you can also
use the option -ignore3 to ignore errors at those read positions. If the
reference genome has lowercase characters for repetitive regions (usually
called softmasked), these regions can be directly filtered using the option
-ignoreLowerCaseRef. These regions can also be filtered at later stages of
the analysis using the FilterVCF command.



-----------------------------------------------------------------
Calling variants on individual samples with the variants detector
-----------------------------------------------------------------

This is the classic module of NGSEP to call SNVs, small indels and structural
variants from sequencing data of single individuals. Basic usage requires an
alignments file in BAM format, the reference genome that was used to produce the
alignments, and a prefix for the output files.

USAGE: 

java -jar NGSEPcore.jar FindVariants <OPTIONS> <REFERENCE> <INPUT_FILE> <OUTPUT_PREFIX>

OPTIONS:	
	-h FLOAT		: Heterozygosity rate. Default: 0.001
	-querySeq STRING	: Call variants just for this sequence name 
	-first INT		: Call variants just from this position in the
				  given query sequence
	-last INT		: Call variants just until this position in the
				  given query sequence
	-ignoreLowerCaseRef	: Ignore sites where the reference allele is
				  lower case. 
	-maxAlnsPerStartPos INT	: Maximum number of alignments allowed to start
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
	-csb		: Calculate a exact fisher test p-value for strand bias between
				  the reference and the alternative allele
	-minQuality	INT	: Minimum genotype quality to accept a SNV call
				  Genotype quality is calculated as 1 minus the
				  posterior probability of the genotype given
				  the reads (in phred scale). Default: 0
	-maxBaseQS INT		: Maximum value allowed for a base quality
				  score. Larger values will be equalized to
				  this value. This parameter allows to reduce
				  the effect of sequencing errors with high
				  base quality scores. Default: 100
	-ignore5 INT		: Ignore this many base pairs from the 5' end
				  of the reads. Default: 0
	-ignore3 INT		: Ignore this many base pairs from the 3' end
				  of the reads. Default: 0
	-knownSTRs FILE		: File with known short tandem repeats (STRs)
				  This is a text file with at least three
				  columns, chromosome, first position and last
				  position. Positions should be 1-based and
				  inclusive.
	-knownVariants FILE	: VCF file with variants to be genotyped. Only
				  these variants will appear in the output vcf
				  file. With this option homozygous calls to
				  the reference allele will be reported
	-embeddedSNVs		: Flag to call SNVs within STRs. By default,
				  STRs are treated as a single locus and hence
				  no SNV will be called within an STR.
	-minSVQuality INT	: Minimum quality score (in PHRED scale) for
				  structural variants. Default: 20
	-genomeSize INT		: Total size of the genome to use during
				  detection of CNVs. This should be used when
				  the reference file only includes a part of
				  the genome (e.g. a chromosome or a partial
				  assembly)   
	-binSize INT		: Size of the bins to analyze read depth.
				  Default: 100
	-algCNV	STRING		: Comma-separated list of read depth algorithms
				  to run (e.g. CNVnator,EWT). Default: CNVnator
	-maxPCTOverlapCNVs INT	: Maximum percentage of overlap of a new CNV
				  with an input CNV to include it in the output
				  Default: 100 (No filter)
	-maxLenDeletion INT	: Maximum length of deletions that the read-pair
				  analysis can identify. Default: 1000000
	-sizeSRSeed INT		: Size of the seed to look for split-read
				  alignments. Default: 8
	-ignoreProperPairFlag	: With this option, the proper pair flag will
				  not be taken into accout to decide if the ends
				  of each fragment are properly aligned. By
				  default, the distribution of insert length is
				  estimated only taking into account reads with
				  the proper pair flag turned on
	-knownSVs FILE		: File with coordinates of known structural
				  variants in GFF format.
	-minMQ INT		: Minimum mapping quality to call an alignment unique.
				  Default: 20
	-sampleId STRING	: Id of the sample that will appear in the
				  output vcf file
	-ploidy INT		: Ploidy of the sample to be analyzed. 
				  Default 2
	-psp			: Flag to print a header in the VCF file with
				  the id and the ploidy of the sample. The
				  header generated with this option is not a
				  standard VCF header. However, it helps NGSEP
				  to keep track of the ploidy of each sample
				  through downstream analyses
	-runRep			: Turns on the procedure to find repetitive
				  regions based on reads with multiple alignments
	-runRD			: Turns on read depth (RD) analysis to identify
				  CNVs
	-noNewCNV		: Turns off finding new CNVs with the read depth
				  analysis. Input CNVs and repeats will still
				  be genotyped using the RD distribution
	-runRP			: Turns on read pair plus split-read analysis
				  (RP+SR) to identify large indels and inversions
	-noSNVS			: Turns off SNV detection. In this mode, only
				  structural variation will be called

Alignments should be provided in BAM v1 format
(see http://samtools.github.io/hts-specs/SAMv1.pdf for details).
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

WARNING 1: The default minimum genotype quality of the variants detector (0) 
will maximize the number of called variants at the cost of generating some
false positives in samples with small coverage or high sequencing error rates.
For conservative variant calling from whole genome sequencing reads use:

java -jar NGSEPcore.jar FindVariants -maxAlnsPerStartPos 2 -minQuality 40 -maxBaseQS 30 <REFERENCE> <INPUT_FILE> <OUTPUT_PREFIX>

If interested in structural variation, you can add the options to run read
depth (RD) and read pair plus split read (RP+SR) approaches to identify
structural variation:

java -jar NGSEPcore.jar FindVariants -runRD -runRP -maxAlnsPerStartPos 2 -minQuality 40 -maxBaseQS 30 <REFERENCE> <INPUT_FILE> <OUTPUT_PREFIX>

If the error rate towards the three prime end increases over 2% you can also
use the option -ignore3 to ignore errors at those read positions. If the
reference genome has lowercase characters for repetitive regions (usually
called softmasked), these regions can be directly filtered using the option
-ignoreLowerCaseRef. These regions can also be filtered at later stages of
the analysis using the FilterVCF command.

WARNING 2: Since v 3.2.0, for RAD Sequencing or genotype-by-sequencing (GBS)
we now recommend the MultisampleVariantsDetector command described above.
However, the classic per-sample analysis pipeline using this command is still
functional with good quality. For both commands it is important to keep in mind
that using the default value of the parameter to control for PCR duplicates
(maxAlnsPerStartPos) will yield very low sensitivity. We recommend to increase
the value of the parameter to about 100 to retain high sensitivity while
avoiding a severe penalty in memory usage. Also, structural variants should not
be called using these data. The usage for conservative variant calling in
RAD-Seq or GBS samples becomes:

java -jar NGSEPcore.jar FindVariants -maxAlnsPerStartPos 100 -minQuality 40 -maxBaseQS 30 <REFERENCE> <INPUT_FILE> <OUTPUT_PREFIX>

-----------------------------------
Calculating base quality statistics
-----------------------------------

This module takes one or more sets of alignments and a reference genome and
writes to standard output a report counting the number of mismatches with the 
reference for each read position from 5' to 3' end. This report is useful to 
detect sequencing error biases. The usage for this tool is the following:

USAGE:

java -jar NGSEPcore.jar QualStats <OPTIONS> <REFERENCE_FILE> <ALIGNMENTS_FILE>*

OPTIONS:
	-minMQ INT	: Minimum mapping quality to call an alignment unique.
			  Default: 20
						  
The file(s) with alignments must be given in SAM or BAM format and the
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

This module calculates the number of base pairs that are covered by reads at
each coverage level from 1 to a maximum. This statistic is useful to visualize
how uniform was the sequencing process over the genome. The usage is as follows

USAGE:

java -jar NGSEPcore.jar CoverageStats <OPTIONS> <ALIGNMENTS_FILE> <OUTPUT_FILE>

OPTIONS:
	-minMQ INT	: Minimum mapping quality to call an alignment unique.
			  Default: 20


The alignments file must be given in SAM or BAM format. The output is a text
file with three columns:
- Coverage
- Number of reference sites with this coverage (Considering all alignments)
- Number of reference sites with this coverage (Considering only reads with 
  unique alignments)

---------------------------------
Functional annotation of variants
---------------------------------

This module takes a VCF file produced by NGSEP, the reference genome in fasta
format, and a gff3 file with gene annotations related with the given genome
(see http://www.sequenceontology.org/gff3.shtml for details) and generates a
VCF file which includes the functional information related with each variant.
The usage is as follows

USAGE:

java -jar NGSEPcore.jar Annotate <OPTIONS> <VARIANTS_FILE> <TRANSCRIPTOME_MAP> <REFERENCE_FILE>

OPTIONS:
	-u INT	: Maximum bp before a gene to classify a variant as Upstream.
		  Default: 1000
	-d INT	: Maximum bp after a gene to classify a variant as Downstream.
		  Default: 300
        -sd INT : Initial basepairs of an intron that should be considered as
		  splice donor. Default: 2
        -sa INT : Final basepairs of an intron that should be considered as
		  splice acceptor. Default: 2
        -si INT : Initial or final basepairs of an intron that should be
		  considered as part of the splice region. Default: 10
        -se INT : Initial or final basepairs of an exon that should be
		  considered as part of the splice region. Default: 2

The vcf file with functional annotations is written in the standard output.
Annotations are included using the following custom fields in the INFO column:

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

----------------------------------------
Merging variants from individual samples
----------------------------------------

NGSEP can be used to merge variants from different samples into an
integrated VCF file. The pipeline for this purpose is as follows.

The first step is to generate a file including the whole set of variants called
in at least one of the samples. This can be done calling the MergeVariants
command as follows:

USAGE:

java -jar NGSEPcore.jar MergeVariants <SEQUENCE_NAMES_FILE> <OUTPUT_FILE> <VARIANTS_FILE>*


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

The second step is to genotype for each sample the variants produced at the
first step using the variants detector (See FindVariants command). For each
sample, the command to execute at this stage (in conservative mode) should look
like this:

java -jar NGSEPcore.jar FindVariants -maxAlnsPerStartPos 2 -minQuality 40 -maxBaseQS 30 -knownVariants <VARS_FILE> <REFERENCE> <INPUT_FILE> <OUTPUT_PREFIX>

where VARS_FILE is the output file obtained in the first step of the merging
process. At the end, this will produce a second set of vcf files which will
differ from the first set in the sense that they will include calls to the
reference allele. The third step is to join these new vcf files using the
following command:

java -jar NGSEPcore.jar MergeVCF <SEQUENCE_NAMES_FILE> <GENOTYPED_VARIANTS_FILE>*

This command will write to standard output the final vcf file with the genotype
calls for each variant on each sample.

-------------------
Filtering VCF files
-------------------

This module implements different filters on VCF files with genotype
information. It writes to standard output a VCF file with variants passing the
filtering criteria. Since version 2.0.6, the default behavior does not perform
any filtering. The filtering order is as follows: first, it executes the
distance filter (-d option), then the filtering of samples and genotypes (-saf,
-fs, -q and -minC options). Finally, it recalculates the number of samples 
genotyped, the number of alleles called and the MAF to execute the remaining 
filters.

USAGE:

java -jar NGSEPcore.jar FilterVCF <OPTIONS> <INPUT_FILE>

OPTIONS:
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
	-d INT		: Minimum distance between variants.
	-g FILE		: File with the reference genome to calculate the
			  GC-Content of the region surrounding the variant.
	-minGC FLOAT	: Minimum percentage of GC of the 100bp region
			  surrounding the variant.
	-maxGC FLOAT	: Maximum percentage of GC of the 100bp region
			  surrounding the variant.
	-q INT		: Minimum genotyping quality score (GQ field for each
			  genotype call in the vcf file).
	-s		: Keep only biallelic SNVs.
	-fi		: Flag to filter sites in which only one allele is
			  observed.
	-fir		: Flag to filter sites in which only the reference
			  allele is observed.
	-fia		: Flag to filter sites in which only one alternative
			  allele is observed.
	-minI INT	: Minimum number of individuals genotyped to keep the
			  variant.
	-minC INT	: Minimum coverage to keep a genotype call.
	-minMAF FLOAT	: Minimum minor allele frequency over the samples in
			  the VCF file.
	-maxMAF FLOAT	: Maximum minor allele frequency over the samples in
			  the VCF file.
	-minOH FLOAT	: Minimum observed heterozygosity over the samples in
			  the VCF file.
	-maxOH FLOAT	: Maximum observed heterozygosity over the samples in
			  the VCF file.
	-maxCNVs INT	: Maximum number of samples with copy number variation
			  in the region where the variant is located.
	-gene STRING	: Id of the gene or the transcript related with the
			  variant.
	-a STRING	: Types of functional annotations (Missense, Nonsense,
			  Synonymous, etc) related with the variant. More than
			  one annotation can be set as a comma-separated list
	-saf FILE	: File with the ids of the samples to be selected (or
			  filtered, see -fs option). The file should have one
			  line per sample, being the first column the sample
			  id. Other columns in the file are ignored.
	-fs		: Flag to filter the samples provided with the -saf
			  option instead of selecting them. 

----------------------------------
Convert VCF files to other formats
----------------------------------

This module allows to convert genotype calls in VCF format to other formats
commonly used to perform different kinds of analysis.

USAGE:

java -jar NGSEPcore.jar ConvertVCF <OPTIONS> <INPUT_FILE> <OUTPUT_PREFIX>

OPTIONS:
	-printStructure		: Prints input format for structure
	-printFasta		: Prints a virtual multiple sequence alignment
				  in fasta format. Useful to build phylogenetic
				  trees
	-printrrBLUP		: Prints the input files for rrBLUP
	-printMatrix		: Prints genotypes in a simple ACGT format
				  which can be imported to excel 
	-printHapmap		: Prints Hapmap format, which can be used in
				  programs such as Tassel
	-printSpagedi		: Prints the input files for Spagedi
	-printPlink		: Prints the input files for Plink
	-printHaploview		: Prints the input files for Haploview
	-printEmma		: Prints the input files for Emma
	-printPowerMarker	: Prints the input files for Powermarker
	-printEigensoft		: Prints the input files for Eigensoft
	-printFlapjack		: Prints the input files for Flapjack
	-printDarwin		: Prints the input files for DarWin
	-printTreeMix		: Prints the input files for TreeMix
	-printJoinMap		: Prints the input file to build genetic maps
				  with JoinMap
	-printPhase		: Prints the input file for PHASE
	-p1 STRING		: Id of the first parent for conversion to
				  JoinMap
	-p2 STRING		: Id of the second parent for conversion to
				  JoinMap
	-s STRING		: Name of the sequence (chromosome) for
				  conversion to PHASE 
	-p FILE			: File with population assignments for the
				  samples. This should be a two column text
				  file with the sample ids in the first column
				  and the ids of the populations in the second
				  column. Required for conversion to TreeMix

WARNING: FASTA convertion does not use IUPAC codes, heterozygous SNPs are 
changed to N.
WARNING 2: Plink is only designed for humans, therefore it will only work for
22 sequences (chromosomes). If a sample exceeds this number, it is convenient 
to reduce the number of chromosomes and to remove all scaffolds.
WARNING 3: To generate dendograms in Tassel, it is better to use the HapMap 
format.

------------------------------
Calculating summary statistics
------------------------------

This module writes to the standard output a report with the numbers of variants
included in a VCF file for different categories. Although it can be called for
any VCF file generated by the pipeline, this report is specially useful when a
complete population is being processed and merged into a single annotated file.

USAGE:

java -jar NGSEPcore.jar SummaryStats <OPTIONS> <INPUT_FILE> 

OPTIONS:
	-m INT			: Minimum number of samples genotyped to
				  accurately calculate the minor allele
				  frequency. Default: 20

-----------------------------------------
Calculating diversity statistics per site
-----------------------------------------

This module produces basic diversity statistics for each variant in a VCF file.
It receives a VCF file and an optional text file with population assignments for
each sample included in the VCF and writes to the standard output the
coordinates of each variant plus the following statistics separated by
semicolon:

1. Number of samples genotyped
2. Expected heterozygosity (under HWE)
3. Observed heterozygosity
4. F-statistic (1-OH/EH)
5. Minor allele frequency (MAF)
6. Chi-square value of departure from HWE
7. Uncorrected p-value of the Chi-square test for departure from HWE

If a file with population assignments is provided, this module will output one
column of statistics for the whole group and one column for each population.

USAGE:

java -jar NGSEPcore.jar DiversityStats <VCF_FILE> <POPULATIONS_FILE>


The populations file is a tab-delimited text file with two columns: sample id
and population id.

-------------------------------------------------------
Calculation of genetic distance matrices from VCF files
-------------------------------------------------------

Generates a distance matrix from a variants file in VCF format. The matrix is
calculated using the basic IBS (Identity by state) algorithm. However, four
options to infer the genotype call information are implemented. In particular,
users can choose predicted allele dosages of CNVs or direct estimations of
allele dosage per site per individual based on relative allele-specific read
counts. The latter option is useful to improve distance estimations in
polyploids. It writes to standard output the matrix of genetic distances in a
generic format.

USAGE:

java -jar NGSEPcore.jar VCFDistanceMatrixCalculator <OPTIONS> <VCF_FILE>

OPTIONS:

	-t INT	: Matrix output format, 0 is full matrix, 1 lower-left matrix
		  and 2 is upper right matrix. Default: 0
	-s INT	: Source of information in the VCF file to calculate distances.
		  0 for simple genotype calls (GT format field), 1 for allele
		  copy number (ACN format field), 2 for total copy number
		  (total of ACN format field), and 3 for raw allele depth (ADP
		  or BSDP format fields). Default: 0
	-p INT	: Default ploidy of the samples. Used if the distance source
		  (-s option) is the raw allele depths to recalculate allele
		  dosage based on these counts. Default: 2

--------------------------------------------------------
Building dendograms using the Neighbor-Joining algorithm
--------------------------------------------------------

Given a distance matrix file, this command builds a dendogram for graphical
display of genetic distances using the Neighbor Joining algorithm. The distance
matrix can be provided as an upper, lower or full matrix. The dendogram is
written to standard output in Newick format.

USAGE:

java -jar NGSEPcore.jar NeighborJoining <MATRIX_FILE>


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
spaces). 

USAGE:

java -jar NGSEPcore.jar AlleleSharingStats <VCF_FILE> <POPULATIONS_FILE> <GROUP_1> <GROUP_2>

OPTIONS:
	-t FILE		: GFF3 file with the transcriptome. If this file is
			  provided, statistics will be provided by gene and not
			  by window.
        -i		: If set, introns will be included in the calculation
			  of pairwise differences. Only useful if the -t option
			  is set.
        -w INT		: Length of each genomic window to calculate pairwise
			  differences between samples. Default: 100000
        -s INT		: Step between windows to calculate pairwise
			  differences between samples. Default: 10000



The populations file is a tab-delimited text file with two columns:
sample id and population id. Writes to the standard output a tab-delimited
report with the following fields:

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
id and gene id respectively



-------------------
Comparing VCF files
-------------------
This module allows to compare the genotype calls included in two different VCF
files, calculating the number and percentage of homozygous and heterozygous
differences between every pair of samples. This writes to standard output a
tab-delimited report with the following fields:

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

USAGE:

java -jar NGSEPcore.jar CompareVCF <OPTIONS> <REFERENCE_FILE> <FIRST_VCF_FILE> <SECOND_VCF_FILE>

OPTIONS: 

	-g FLOAT	: Minimum percentage (0-100) of variants genotyped in
			  both samples. Default: 50.
	-d FLOAT	: Maximum percentage (0-100) of differences between the
			  pair of samples. Default: 5.

The first required argument is the FASTA file with the reference genome used to
generate the VCF files. Default values of optional parameters are set to
facilitate the detection of duplicated (or very similar) samples. To report the
complete set of sample pairs, use -g 0 -d 100. 

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
expected heterozygosity. In the latter mode, the model will only produce
homozygous genotypes

USAGE:

java -jar NGSEPcore.jar ImputeVCF <OPTIONS> <VCF_FILE> <OUT_PREFIX>

OPTIONS: 

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
	-w INT		: Size of the window to process variants at the same
			  time. Default: 5000
        -o INT		: Overlap between windows. Default: 50
	-c FLOAT	: Estimated average number of centiMorgans per Kbp on
			  euchromatic regions of the genome. This value is used
			  by the model to estimate initial transitions between
			  the states of the HMM. Typical values of this
			  parameter are 0.001 for human populations, 0.004 for
			  rice and 0.35 for yeast populations. We expect to
			  implement an option to allow setting the estimated
			  recombination rate per site in future versions.
			  Default: 0.001
	-t		: If set, transition probabilities in the HMM will NOT
			  be updated during the Baum-Welch training of the HMM. 
			  Not recommended unless the -c option is set to a
			  value allowing a reasonable initial estimation of the
			  transition probabilities.
        -ip		: Specifies that parents of the population are inbreds
        -is		: Specifies that the samples to impute are inbreds


This module outputs two files, the first is a VCF file including the imputed
genotypes for the datapoints having an undecided genotype call in the input
file. The second outputs for each SNP and each sample the index of the parent
that most likely originated the observed haplotype of the individual.

--------------------------------
Finding haplotype introgressions
--------------------------------

This module runs a window-based analysis to identify the most common haplotype
within each of the populations described in the given populations file and then
identifies common haplotypes of one population introgressed in samples of a
different population. Although it can be run on any VCF file, it is
particularly designed to work with populations of inbred samples.

USAGE:

java -jar NGSEPcore.jar IntrogressionAnalysis <OPTIONS> <VCF_FILE> <POPULATIONS_FILE> <OUT_PREFIX>

OPTIONS:

	-p FLOAT	: Minimum percentage of samples genotyped within a
			  population to identify the most common allele.
			  Default: 80
	-d FLOAT	: Minimum difference between reference allele
			  frequencies of at least two populations to consider a
			  variant discriminative. Default: 0.6
	-m FLOAT	: Maximum minor allele frequency within a population to
			  consider the major allele of a variant as
			  representative allele for such population.
			  Default: 0.4
	-w INT		: Window size as number of variants within each window
			  Default: 50
	-o INT		: Overlap as number of variants shared between neighbor
			  windows Default: 0
	-a INT		: Score given of a match between homozygous genotypes
			  comparing haplotypes Default: 1
	-i INT		: Score given of a mismatch between homozygous
			  genotypes comparing haplotypes Default: -1
	-s INT		: Minimum score to match an individual haplotype with a
			  population-derived haplotype Default: 30
	-v		: Outputs a VCF file with the biallelic variants that
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
   assigned to more than one population

------------------------------
Building genomes from variants
------------------------------

This module takes a VCF file with genotype information from one sample and the
reference genome used to build the VCF and generates a new genome in fasta
format modified using the alternative alleles from variants called as
homozygous alternative within the individual. This can be useful to perform
polishing of new genome assemblies using Illumina data, or in general to
construct a haploid version of an individual genome.

USAGE:

java -jar NGSEPcore.jar VCFIndividualGenomeBuilder <VCF_FILE> <REFERENCE_GENOME> <OUT_GENOME>

-----------------
Comparing genomes
-----------------

This module takes two assembled genomes in fasta format and their corresponding
transcriptome gene annotations in GFF3 format and runs a whole genome
comparison taking unique genes as orthology units. It also calculate paralogs
within each genome. 

USAGE:

java -jar NGSEPcore.jar GenomesAligner <OPTIONS> <GENOME1> <TRANSCRIPTOME1> <GENOME2> <TRANSCRIPTOME2>

OPTIONS:

        -o STRING : Prefix of output files Default: genomesAlignment

The output is a series of text files having the ids and physical coordinates of
the paralogs within each genome and the orthologs between the two genomes.
The ortholog files, called <PREFIX>_uniqueG1.tsv and <PREFIX>_uniqueG2.tsv,
have the following format:

1. Id of the gene in the first genome
2. Chromosome of the gene in the first genome
3. Start of the gene in the first genome
4. End of the gene in the first genome
5. Id of the ortholog in the second genome
6. Chromosome of the ortholog in the second genome
7. Start of the ortholog in the second genome
8. End of the ortholog in the second genome
9. Alignment type. It can be "L" if the gene has a unique ortholog in the
   second genomeand it makes part of a synteny block. "U" if the gene has
   a unique ortholog but it does not make part of the syntheny block, and
   "M" if the gene has multiple orthologs in the second genome.

The files with the paralogs, called <PREFIX>_paralogsG1.tsv and
<PREFIX>_paralogsG2.tsv, have the same 8 first columns but columns 5 to 8
contain genes within the same genome as genes in column 1 to 4. Finally,
the file called <PREFIX>_linearView.html can be loaded in a web browser
and provides an interactive view of the alignment based on the d3 web
development technology (https://d3js.org/).

-------------------
Demultiplexing reads
-------------------

This option allows to build individual fastq files for different samples from
a single file containing the reads for a whole sequencing lane in which several
samples were barcoded and sequenced.

USAGE:

java -jar NGSEPcore.jar Demultiplex <OPTIONS> <INDEX_FILE> <FASTQ_FILE_1> (<FASTQ_FILE_2>) 

OPTIONS: 
	-o DIRECTORY	: Directory where the output fastq files will be saved
	-t STRING	: If this sequence is found within a read, the read
			  will be trimmed up to the start of this sequence.
			  useful to remove adapter contamination
	-u		: Output uncompressed files
        -a		: Activate demultiplexing with dual barcoding.
	-d FILE		: Tab-delimited file storing physical locations of the
			  files to be demultiplexed. Columns of the file should
			  be Flowcell, lane and fastq file (which can be gzip
			  compressed). A second fastq file can be specified if
			  the lane was sequenced in paired-end mode. If the
			  reads sequenced for one lane are split in multiple
			  files, each file (or each pair of files) should be
			  included in a separate row. If this option is used,
			  the options -f, -l and the input fastq file(s) are
			  ignored.
	-f STRING	: Id of the flowcell corresponding to the input fastq
			  file(s). Ignored if the -d option is specified but
			  required if the -d option is not specified.
	-l STRING	: Id of the lane corresponding to the input fastq
			  file(s). Ignored if the -d option is specified but
			  required if the -d option is not specified.

INDEX_FILE is a tab-delimited text file with four columns by default: flowcell,
lane, barcode and sampleID. If the -a option for dual barcode is activated,
five columns are expected: flowcell, lane, barcode1, barcode2 and sampleID. The
file must have a header line. The same index file can be used to demultiplex
several FASTQ files. Out FASTQ files will be gzip compressed by default.

------------------------------------
Comparing read depth between samples
------------------------------------

This function compares the read depth of two samples. It takes two alignment
files and a reference genome, splits the genome into windows, and for each
window compares the read depth between the two samples. It outputs a text file
containing the list of windows of the genome in which the normalized read depth 
ratio between the two samples is significantly different from 1. The text file
contains the following columns:

1. Chromosome
2. Window start
3. Window end
4. Read depth sample 1
5. Read depth sample 2
6. Normalized read depth ratio
7. P-value

USAGE:

java -jar NGSEPcore.jar CompareRD <OPTIONS> <ALIGNMENTS_FILE_1> <ALIGNMENTS_FILE_2> <REFERENCE> <OUT_PREFIX>

OPTIONS:

	-binSize INT	: Window size to be used during the read depth
			  comparison. Default: 100
	-p FLOAT	: Maximum p-value. Only the windows with a p-value
			  lower than that specified will be reported.
			  Default: 0.001
	-w		: Output an entry for every window in the genome
	-g		: Perform GC-correction of the read depth
	-b		: Perform the Bonferroni correction for multiple 
			  testing

----------------------------------------
Obtaining k-mers spectrum from sequences
----------------------------------------

Generate a distribution of k-mer abundances from a file of DNA sequences either
in fastq or in fasta format. Writes to standard output the number of k-mers
obtained at each specific read depth.

USAGE:

java -jar NGSEPcore.jar KmersCounter <OPTIONS> <SEQUENCES_FILE>

OPTIONS:

        -b     : Count k-mers from both strands.
        -k INT : K-mer length. Default: 21
        -fasta : Input is a fasta file.

----------------------------------------------------------
Obtaining relative allele counts from read alignment files
----------------------------------------------------------

Calculates a distribution of relative allele counts for sites showing base calls
for more than one nucleotide from read alignment files in BAM format. This
analysis is useful to predict the ploidy of a sequenced sample.

USAGE:

java -jar NGSEPcore.jar RelativeAlleleCounts <OPTIONS> <ALIGNMENTS_FILE>

OPTIONS:

        -m INT	: Minimum read depth Default: 10
        -M INT	: Maximum read depth Default: 1000
        -q INT	: Minimum base quality score (Phred scale) Default: 20
        -r FILE : File with repeats (or any kind of genomic regions) that
		  should not be taken into account in the analysis. The format
		  of this file should contain three columns: Sequence name
		  (chromosome), first position in the sequence, and last
		  position in the sequence. Both positions are assumed to be
		  1-based.
        -f FILE : File with genomic regions that should be taken into account
		  in the analysis. The format of this file should contain three
		  columns: Sequence name (chromosome), first position in the
		  sequence, and last position in the sequence. Both positions
		  are assumed to be 1-based.
        -s      : Consider secondary alignments. By default, only primary
		  alignments are processed

----------------------------------------------
Simulating individuals from a reference genome
----------------------------------------------

This simulator takes a (haploid) genome assembly and simulates a single
individual including homozygous and heterozygous mutations (SNPs, indels and
mutated STRs) relative to the input assembly. It produces two files, a fasta
file with the simulated genome, and a phased VCF file with the simulated
variants.

USAGE:

java -jar NGSEPcore.jar SingleIndividualSimulator <OPTIONS> <GENOME> <OUT_PREFIX>

OPTIONS:

	-s DOUBLE	: Proportion of reference basepairs with simulated SNV
			  events. Default: 0.001
	-i DOUBLE	: Proportion of reference basepairs with simulated
			  indel events. Default: 0.0001
	-f DOUBLE	: Fraction of input STRs for which a mutation will be
			  simulated. Default: 0.1
	-t FILE		: Path to a text file describing the known STRs in the
			  given genome
        -u INT		: Zero-based index in the STR file where the unit
			  sequence is located. Default: 14
	-d STRING	: ID of the simulated sample. Appears in the VCF header
			  and as part of the name of the sequences in the
			  simulated genome. Default: Simulated
	-p INT		: Ploidy of the simulated sample. Default: 2

------------------------------
Citing and supporting packages
------------------------------

A manuscript with the description of the main modules of NGSEP is available at
Nucleic Acids research:

Duitama J, Quintero JC, Cruz DF, Quintero C, Hubmann G, Foulquie-Moreno MR, Verstrepen KJ, Thevelein JM, and Tohme J. (2014). 
An integrated framework for discovery and genotyping of genomic variants from high-throughput sequencing experiments. 
Nucleic Acids Research. 42(6): e44. 
http://doi.org/10.1093/nar/gkt1381

A description of some of the latest modules and recent benchmarks with other
tools for variants detection on Genotype-By-Sequencing (GBS) data is available
at BMC Genomics:

Perea C, Hoz JFDL, Cruz DF, Lobaton JD, Izquierdo P, Quintero JC, Raatz B and Duitama J. (2016).
Bioinformatic analysis of genotype by sequencing (GBS) data with NGSEP.
BMC Genomics, 17:498.
http://doi.org/10.1186/s12864-016-2827-7

Details of variant detection algorithms implemented in NGSEP can be found in
the following publications:

SNV detection:
Duitama J, Srivastava PK, and Mandoiu II. (2012). 
Towards accurate detection and genotyping of expressed variants from whole transcriptome sequencing data. 
BMC Genomics, 13(Suppl 2), S6. 
http://doi.org/10.1186/1471-2164-13-S2-S6

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

NOTE: Since version 2.1.2, we implemented a new model to integrate paired-end and split-read analysis for detection of large indels.
Benchmarking with other tools is in progress.

NGSEP is also supported by the following open source software packages:

Bowtie2: http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
Picard: http://picard.sourceforge.net/
Jsci: http://jsci.sourceforge.net/
XChart: http://xeiam.com/xchart/
Trimmomatic: http://www.usadellab.org/cms/?page=trimmomatic. We borrowed one
class from Trimmomatic 0.35 to allow correct reading of gzip files
