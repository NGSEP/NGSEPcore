package ngsep.gbs;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.logging.Logger;

import ngsep.alignments.ReadAlignment;
import ngsep.alignments.ReadAlignmentPair;
import ngsep.alignments.ReadsAligner;
import ngsep.alignments.io.ReadAlignmentFileReader;
import ngsep.alignments.io.ReadAlignmentFileWriter;
import ngsep.genome.GenomicRegionComparator;
import ngsep.genome.ReferenceGenome;
import ngsep.genome.ReferenceGenomeFMIndex;
import ngsep.main.CommandsDescriptor;
import ngsep.main.OptionValuesDecoder;
import ngsep.main.ProgressNotifier;
import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.DNASequence;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.RawRead;
import ngsep.sequences.io.FastaFileReader;
import ngsep.variants.CalledGenomicVariant;
import ngsep.variants.CalledSNV;
import ngsep.variants.GenomicVariant;
import ngsep.variants.GenomicVariantAnnotation;
import ngsep.variants.GenomicVariantImpl;
import ngsep.variants.SNV;
import ngsep.variants.Sample;
import ngsep.vcf.VCFFileHeader;
import ngsep.vcf.VCFFileReader;
import ngsep.vcf.VCFFileWriter;
import ngsep.vcf.VCFRecord;

public class VCFRelativeCoordinatesTranslator {

	// Logging and progress
	private Logger log = Logger.getLogger(VCFRelativeCoordinatesTranslator.class.getName());
	private ProgressNotifier progressNotifier = null;
	
	// Parameters
	private String inputFile;
	private String outputPrefix;
	private ReferenceGenome genome;
	private String fmIndexFile;
	private String filenameConsensusFA;
	private String filenameAlignmentBAM;
	
	
	//Statistics
	int numTranslatedRecords = 0;
	int totalRecords = 0;
	int skippedRecords = 0;
	int nullReadName = 0;
	int TruePosNotAlign = 0;
	
	int consensusTot = 0;
	int pairedConsensus = 0;
	int singleConsensus = 0;
	int biallelic = 0;
	int triallelic = 0;
	
	
	// issues with translation
	int untranslated = 0;
	int notSNV = 0;
	int nonVariant = 0;
	int refNotInAlleles = 0;
	int notRefSeq = 0;
	int notDNA = 0;
	int translatedmapped = 0;
	int refSeqLess0 = 0;
	int trueCallsNull = 0;
	int recordWihoutAlign = 0;
	
	// issues with mapping
	int unmappedReadSingle = 0;
	int unmappedReadPaired = 0;
	int unmappedRead = 0;
	int singlemapfor = 0;
	int singlemaprev = 0;
	int unpaired = 0;
	int partial = 0;
	int noAlignSingle = 0;
	int noAlignPaired = 0;
	int oddAlign = 0;
	int unmappedCluster = 0;

	
	// Get and set methods
	public Logger getLog() {
		return log;
	}
	public void setLog(Logger log) {
		this.log = log;
	}
	
	public ProgressNotifier getProgressNotifier() {
		return progressNotifier;
	}
	public void setProgressNotifier(ProgressNotifier progressNotifier) { 
		this.progressNotifier = progressNotifier;
	}
	
	public String getInputFile() {
		return inputFile;
	}
	public void setInputFile(String inputFile) {
		this.inputFile = inputFile;
	}
	
	public ReferenceGenome getGenome() {
		return genome;
	}
	public void setGenome(ReferenceGenome genome) {
		this.genome = genome;
	}
	public void setGenome(String genomeFile) throws IOException {
		this.genome = OptionValuesDecoder.loadGenome(genomeFile, log);
	}
	
	public String getOutputPrefix() {
		return outputPrefix;
	}
	public void setOutputPrefix(String outputPrefix) {
		this.outputPrefix = outputPrefix;
	}
	
	public String getFmIndexFile() {
		return fmIndexFile;
	}
	public void setFmIndexFile(String fmIndexFile) {
		this.fmIndexFile = fmIndexFile;
	}
	
	public String getFilenameConsensusFA() {
		return filenameConsensusFA;
	}
	public void setFilenameConsensusFA(String filenameConsensusFA) {
		this.filenameConsensusFA = filenameConsensusFA;
	}
	
	public String getFilenameAlignmentBAM() {
		return filenameAlignmentBAM;
	}
	public void setFilenameAlignmentBAM(String filenameAlignmentBAM) {
		this.filenameAlignmentBAM = filenameAlignmentBAM;
	}
	
	public static void main(String[] args) throws Exception {
		VCFRelativeCoordinatesTranslator instance = new VCFRelativeCoordinatesTranslator();
		CommandsDescriptor.getInstance().loadOptions(instance, args);
		instance.run();
	}

	public void run() throws IOException {
		if(inputFile==null) throw new IOException("The VCF file with the variants to translate is a required parameter");
		if (genome == null) throw new IOException("The file with the reference genome is a required parameter");
		if(outputPrefix==null) throw new IOException("The prefix of the output files is a required parameter");
		
		logParameters();
		Map<String, ReadAlignment> alignments;
		if (filenameAlignmentBAM!=null) alignments = loadAlignmentsFromBAM(filenameAlignmentBAM);
		else if (filenameConsensusFA!=null) {
			log.info("Aligning input consensus sequences from "+filenameConsensusFA);
			alignments = alignConsensusSequences(filenameConsensusFA);
			saveAlignments(alignments,outputPrefix+"_alns.bam");
		} else {
			throw new IOException("Either a fasta file with the consensus sequences or a BAM fie with alignments is required");
		}
		log.info("Loaded a total of " + alignments.size() + " alignments.");
		
		translate(inputFile, alignments, outputPrefix+".vcf");
		printStatistics(outputPrefix + ".info");
		log.info("Process finished");
	}
	
	private void logParameters() {
		ByteArrayOutputStream os = new ByteArrayOutputStream();
		PrintStream out = new PrintStream(os);
		out.println("Input VCF file:"+ inputFile);
		out.println("Loaded reference genome from: "+genome.getFilename());
		out.println("Output prefix:"+ outputPrefix);
		if(filenameAlignmentBAM!=null) out.println("Aligned consensus sequenced will be loaded from:"+ filenameAlignmentBAM);
		else if (filenameConsensusFA!=null)  out.println("Fasta fie with consensus sequences:"+ filenameConsensusFA);
		if(fmIndexFile!=null) out.println("FM-index file:"+ fmIndexFile);
		log.info(os.toString());
	}
	
	

	/**
	 *  For each record on the RelativeVCF, this method translates and rewrites. 
	 * @throws Exception
	 */
	public void translate(String inputFile, Map<String, ReadAlignment> alignments, String outputFile) throws IOException {
		List<VCFRecord> translatedRecords = new ArrayList<>();
		List<Sample> samples;
		VCFFileHeader header = VCFFileHeader.makeDefaultEmptyHeader();
		
		try (VCFFileReader vcfOpenFile = new VCFFileReader(inputFile)) {
			samples = vcfOpenFile.getHeader().getSamples();
			for(Sample s:samples) {
				header.addSample(s,header.getSamplesWithHeaderLine().containsKey(s.getId()));
			}
			
			Iterator<VCFRecord> vcfReader = vcfOpenFile.iterator();
			// Iterate over vcfRecords
			while(vcfReader.hasNext()) {
				VCFRecord translatedRecord = null;
				VCFRecord record = vcfReader.next();
				
				
					
				
				String clusterID = record.getSequenceName();
				ReadAlignment alignment =  alignments.get(clusterID);
				if(alignment != null ) {
					// If variant is SNP
					if(record.getVariant().isSNV()) {
						translatedRecord = translateRecord(alignment, record, header);
					} else {
						notSNV++;
					}
					if(translatedRecord != null) {
						translatedRecords.add(translatedRecord);
						numTranslatedRecords++;
					} else {
						untranslated++;
					}
				} else recordWihoutAlign++;
				totalRecords++;
			}
			
		}
		Collections.sort(translatedRecords,new GenomicRegionComparator(genome.getSequencesMetadata()));
		VCFFileWriter writer = new VCFFileWriter ();
		try (PrintStream mappedVCF = new PrintStream(outputFile)) {
			writer.printHeader(header, mappedVCF);
			writer.printVCFRecords(translatedRecords, mappedVCF);
		}
	}
	
	public void printStatistics (String outputFile) throws IOException {
		try (PrintStream info = new PrintStream(outputFile)) {
			info.println("Total number of records in relative VCF: " + totalRecords);
			info.println("Number of translated records: " + numTranslatedRecords);
			info.println("Number of translater biallelic variants: " + biallelic);
			info.println("Total number of consensus sequences: " + consensusTot);
			info.println("Single Consensus: " + singleConsensus);
			info.println("Paired Consensus: " + pairedConsensus);
			info.println("------ Issues with translation ------");
			info.println("Number of records without an alignment: "+recordWihoutAlign);
			info.println("Number of records with an alignment: " + translatedmapped);
			info.println("Number of records not translated even though they had an alignment: " + untranslated);
			info.println("Number of records that are triallelic variants: " + triallelic);
			info.println("Number of records whose Cluster ID was unmapped: " + unmappedCluster);
			info.println("Number of records where matching reference sequence is not DNA: " + notDNA); 
			info.println("Number of records that are not SNV: " + notSNV);
			info.println("Number of records where no reference sequence was found: " + notRefSeq);
			info.println("Number of records where reference sequence does not exist (-1): " + refSeqLess0);
			info.println("Number of records where no calls found (matches number of triallelic): " + trueCallsNull);
			info.println("------ Issues with mapping ------");
			info.println("Number of unmapped consensus: " + unmappedRead);
			info.println("Number of unmapped single consensus: " + unmappedReadSingle);
			info.println("Number of unmapped paired consensus: " + unmappedReadPaired);
			info.println("------ Issues with paired mapping ------");
			info.println("Unmapped paired consensus: " + unmappedReadPaired);
			info.println("Number of consensus where only reverse consensus aligned: " + singlemapfor);
			info.println("Number of consensus where only forward consensus aligned: " + singlemaprev);
			info.println("Number of consensus where neither forward nor reverse aligned: " + noAlignPaired);
			info.println("Number of unpaired consensus: " + unpaired);
			info.println("Number of partial consensus: " + partial);
			info.println("Odd alignment: "+oddAlign);
		}
	}
	
	/**
	 * This method takes in an alignment to the reference and a record 
	 * and translates each record to match the reference.
	 * @param algn
	 * @param record
	 * @return translated record
	 */
	private VCFRecord translateRecord(ReadAlignment algn, VCFRecord record, VCFFileHeader header) {
		translatedmapped++;
		VCFRecord translatedRecord = null;
		List<CalledGenomicVariant> trueCalls = new ArrayList<>();

		// True position of variant with respect to reference forward strand (1-based)
		int truePos;
		
		// Reference nucleotide as interpreted by the alignment to the reference
		char trueRef;
		
		// Position of variant with respect to the de-novo alignment 1-based
		// Relative reference base (the base at position of variant in the consensus seq)
		int zeroBasedRelativePos = record.getFirst()-1;
		
		// Variant as found in de-novo vcf
		GenomicVariant relativeVar = record.getVariant();
		
		// Alternative Allele based on de-novo VCF
		String[] relativeAlleles = relativeVar.getAlleles();
		
		//For each relative allele, position in the translated list of alleles
		Map<String,Integer> relAllelesRelPos = new HashMap<String, Integer>(relativeAlleles.length);
		Map<String,Integer> relAllelesTranslationPos = new HashMap<String, Integer>(relativeAlleles.length);
		
		int[] filedsFormat = record.getFieldsFormat();
		
		// IF the reference is not in the relativeAlleles, the variant is likely triallelic
		boolean refInRelativeAllels = false;
		
		String seqName = algn.getSequenceName();
		
		if(seqName==null) {
			unmappedCluster++;
			unmappedRead++;
			return null;
		}
		
		
		//-1 if the read position does not align with the reference 
		truePos = algn.getReferencePosition(zeroBasedRelativePos);
		
		if(truePos<=0) {
			System.out.println("Variant without ref pos: "+relativeVar.getSequenceName()+":"+relativeVar.getFirst()+" aln seq name: "+seqName+" pos: "+truePos+" aln "+algn);
			refSeqLess0++;
			return null;
		}
		
		CharSequence trueRefSeq = genome.getReference(seqName, truePos, truePos);
		if(trueRefSeq == null) {
			notRefSeq++;
			return null;
		}
		
		if(!DNASequence.isDNA(trueRefSeq)) {
			notDNA++;
			return null;
		} else {
			trueRef = trueRefSeq.charAt(0);
		}
		//TODO: Translate indels
		List<String> refBasedAlleles = new ArrayList<>();
		refBasedAlleles.add(""+trueRef);
		for(int i=0;i<relativeAlleles.length;i++) {
			String allele = relativeAlleles[i];
			if(!DNASequence.isDNA(allele)) continue;
			if(algn.isNegativeStrand()) {
				allele = DNASequence.getReverseComplement(allele).toString();
			}
			relAllelesRelPos.put(allele, i);
			if(allele.charAt(0)==trueRef) {
				refInRelativeAllels = true;
				relAllelesTranslationPos.put(allele, 0);
			} else if(!refBasedAlleles.contains(allele)) {
				relAllelesTranslationPos.put(allele, refBasedAlleles.size());
				refBasedAlleles.add(allele);
			}
		}
		GenomicVariant variant;
		if(refBasedAlleles.size()== 2) {
			variant = new SNV(seqName, truePos, trueRef, refBasedAlleles.get(1).charAt(0));
			biallelic++;
		} else if (refBasedAlleles.size()>= 3) {
			variant = new GenomicVariantImpl(seqName, truePos, refBasedAlleles);
			triallelic++;
		} else {
			nonVariant++;
			return null;
		}
		variant.setVariantQS(relativeVar.getVariantQS());
		if(!refInRelativeAllels) {
			refNotInAlleles++;
		}
		// get the calls for this record.
		List<CalledGenomicVariant> calls = record.getCalls();		
			
		// Translate genotypes 
		//-1 for undecided, 0 for homozygous reference, 1 for heterozygous, 2 for homozygous variant
		for(CalledGenomicVariant relativeCall: calls) {
			String [] calledAlleles = relativeCall.getCalledAlleles();
			int [] acgtCounts = relativeCall.getAllCounts();
			short [] relACN = relativeCall.getAllelesCopyNumber();
			short totalCN = relativeCall.getCopyNumber();
			if (algn.isNegativeStrand()) {
				for(int i=0;i<calledAlleles.length;i++) {
					calledAlleles[i] = DNASequence.getReverseComplement(calledAlleles[i]).toString();
				}
				if (acgtCounts!=null) {
					int tmp = acgtCounts[0];
					acgtCounts [0] = acgtCounts[3];
					acgtCounts[3] = tmp;
					tmp = acgtCounts[1];
					acgtCounts [1] = acgtCounts[2];
					acgtCounts[2] = tmp;
				}
			}
			short [] acn = new short[refBasedAlleles.size()];
			Arrays.fill(acn, (short)0);
			for(int i=0;i<calledAlleles.length;i++) {
				Integer pos = relAllelesTranslationPos.get(calledAlleles[i]);
				Integer relIdx = relAllelesRelPos.get(calledAlleles[i]);
				//if("0".equals(relativeVar.getSequenceName())) System.out.println("Pos: "+relativeVar.getFirst()+" Indv: "+relativeCall.getSampleId()+" Called allele "+calledAlleles[i]+" relAllPos: "+relIdx+" absPos: "+pos);
				if(pos!=null && relIdx!=null && pos<acn.length) acn[pos] = relACN[relIdx];
			}
			if(variant instanceof SNV) {
				byte genotype = CalledGenomicVariant.GENOTYPE_UNDECIDED;
				if(calledAlleles.length==2) {
					genotype = CalledGenomicVariant.GENOTYPE_HETERO;
					if(acn[0]==0 || acn[1]==0 || acn[0]+acn[1]!=totalCN) log.warning("Weird number of called alleles for variant "+relativeCall.getSequenceName()+":"+relativeCall.getFirst()+" heterozygous genotype for: "+relativeCall.getSampleId()+" acnCounts: "+acn[0]+" "+acn[1]+" relative counts: "+relACN[0]+" "+relACN[1]);
				}
				else if (calledAlleles.length==1) {
					if(calledAlleles[0].charAt(0)!=trueRef) {
						genotype = CalledGenomicVariant.GENOTYPE_HOMOALT;
						if(acn[0]!=0 || acn[1]==0 || acn[1]!=totalCN) log.warning("Weird number of called alleles for variant "+relativeCall.getSequenceName()+":"+relativeCall.getFirst()+" homozygous alternative genotype for: "+relativeCall.getSampleId()+" acnCounts: "+acn[0]+" "+acn[1]+" relative counts: "+relACN[0]+" "+relACN[1]);
						acn[1]=totalCN;
						
					} else {
						genotype = CalledGenomicVariant.GENOTYPE_HOMOREF;
						if(acn[0]==0 || acn[1]!=0 || acn[0]!=totalCN) log.warning("Weird number of called alleles for variant "+relativeCall.getSequenceName()+":"+relativeCall.getFirst()+" homozygous reference genotype for: "+relativeCall.getSampleId()+" acnCounts: "+acn[0]+" "+acn[1]+" relative counts: "+relACN[0]+" "+relACN[1]);
						acn[0]=totalCN;
					}
				}
				CalledSNV trueCall = new CalledSNV((SNV)variant, genotype);
				trueCall.setGenotypeQuality(relativeCall.getGenotypeQuality());
				trueCall.setTotalReadDepth(relativeCall.getTotalReadDepth());
				trueCall.setAllelesCopyNumber(acn);
				if(acgtCounts!=null)  trueCall.setAllBaseCounts(acgtCounts);
				trueCalls.add(trueCall);
			}
		}
		if(trueCalls.size()>0) {
			
			translatedRecord = new VCFRecord(variant, filedsFormat, trueCalls, header);
			translatedRecord.addAnnotation(new GenomicVariantAnnotation(variant, "DENOVOCLUSTER", relativeVar.getSequenceName()));
			translatedRecord.addAnnotation(new GenomicVariantAnnotation(variant, "DENOVOCLUSTERPOS", relativeVar.getFirst()));
			translatedRecord.addAnnotation(new GenomicVariantAnnotation(variant, "DENOVOCLUSTERCONSENSUS", relativeVar.getReference()));
			translatedRecord.updateDiversityStatistics();
		} else trueCallsNull++;
		return translatedRecord;
	}
	
	public Map<String, ReadAlignment> alignConsensusSequences(String consensusSequencesFile) throws IOException {
		String debugConsensusName = null;
		Map<String, ReadAlignment> alignments = new HashMap<String, ReadAlignment>();
		
		ReadsAligner aligner = new ReadsAligner();
		aligner.setGenome(genome);
		if(fmIndexFile!=null) {
			aligner.setFmIndex(ReferenceGenomeFMIndex.load(genome, fmIndexFile));
		}
		else aligner.setFmIndex(new ReferenceGenomeFMIndex(genome, log));
		
		String pairedEndAnchor = ReadCluster.MIDDLE_N_SEQUENCE_PAIRED_END;
		int numNCharsPairedEnd = pairedEndAnchor.length();	
		
		
		try (FastaFileReader reader = new FastaFileReader(consensusSequencesFile)) {
			reader.setLog(log);
			Iterator<QualifiedSequence> it = reader.iterator();
			for(int i=0;it.hasNext();i++) {
				consensusTot++;
				QualifiedSequence consensus = it.next();
				String consensusName = consensus.getName();
				if(consensusName.equals(debugConsensusName)) System.err.println(consensusName);
				if(consensusName.startsWith("Cluster_")) {
					consensusName = consensusName.substring(8);
				}
				String seq = consensus.getCharacters().toString();
				if(seq.length() == 0) continue;
				
				int indexN = seq.indexOf(pairedEndAnchor);
				if(indexN <=30 || indexN>=seq.length()-30) {
					singleConsensus++;
					RawRead read = new RawRead(consensusName, consensus.getCharacters(),RawRead.generateFixedQSString('5', consensus.getLength()));
					List<ReadAlignment> alns = aligner.alignRead(read, true, false);
					if((i+1)%10000==0) log.info("Aligning consensus sequence "+(i+1)+" id: "+consensus.getName()+" sequence: "+seq+" alignments: "+alns.size()+". Total unmapped "+unmappedRead);
					if(alns.size()==0) {
						unmappedReadSingle++;
						unmappedRead++;
						continue;
					}
					ReadAlignment first = alns.get(0);
					if(consensusName.equals(debugConsensusName)) System.err.println("First algn single: " + first.getReadCharacters() + "\t" + first.getFlags() + "\t" + first.getFirst() + "\t" + first.getCigarString());
					alignments.put(consensusName,first);
				} else {
					DNAMaskedSequence seq1 = new DNAMaskedSequence(seq.substring(0,indexN));
					DNAMaskedSequence seq2 = new DNAMaskedSequence(seq.substring(indexN+numNCharsPairedEnd));
					
					RawRead read1 = new RawRead(consensusName, seq1, RawRead.generateFixedQSString('5', seq1.length()));
					RawRead read2 = new RawRead(consensusName, seq2.getReverseComplement(), RawRead.generateFixedQSString('5', seq2.length()));
					List<ReadAlignment> alns1 = aligner.alignRead(read1, true, true);
					List<ReadAlignment> alns2 = aligner.alignRead(read2, true, true);
					if((i+1)%10000==0) log.info("Aligning consensus sequence "+(i+1)+" id: "+consensus.getName()+" sequence: "+seq+" alignments 1: "+alns1.size()+" alignments 2: "+alns2.size()+" Total unmapped "+unmappedRead);
					pairedConsensus++;
					if(alns1.size()==0|| alns2.size()==0) {
						if(alns1.size()==0 && alns2.size()!=0) singlemapfor++;
						else if(alns1.size()!=0 && alns2.size()==0) singlemaprev++;
						else noAlignPaired++;
						unmappedReadPaired++;
						unmappedRead++;
						continue;
					}
					List<ReadAlignmentPair> pairs = aligner.findPairs(alns1, alns2, true);
					if(pairs.size()==0) {
						unpaired++;
						unmappedReadPaired++;
						unmappedRead++;
						continue;
					}
					ReadAlignmentPair first = pairs.get(0);
					ReadAlignment aln1 = first.getAln1();
					ReadAlignment aln2 = first.getAln2();
					if(consensusName.equals(debugConsensusName)) System.err.println("First algn aln1: " + aln1.getReadCharacters() + "\t" + aln1.getFlags() + "\t" + aln1.getFirst() + "\t" + aln1.getCigarString());
					if(consensusName.equals(debugConsensusName)) System.err.println("Second algn aln2: " + aln2.getReadCharacters() + "\t" + aln2.getFlags() + "\t" + aln2.getFirst() + "\t" + aln2.getCigarString());
					if(aln1.isPartialAlignment(10) || aln2.isPartialAlignment(10)) {
						partial++;
						unmappedReadPaired++;
						unmappedRead++;
						continue;
					}
					int last1 = aln1.getLast(); 
					int last2 = aln2.getLast();
					if(aln1.getFirst()<last2) {
				 		int internalSoftClip1 = aln1.getSoftClipEnd();
				 		int internalSoftClip2 = aln2.getSoftClipStart();
						String cigar = aln1.getCigarString();
						int last = last2;
						
						int nextRef = last1+internalSoftClip1+numNCharsPairedEnd+1;
						int noSoftClipAln2First = aln2.getFirst()-internalSoftClip2;
						if(nextRef<=noSoftClipAln2First) {
							if(internalSoftClip1>0) {
								//Transform into M
								cigar = cigar.substring(0,cigar.length()-1);
								cigar+="M";
							}
							//Add N matches
							cigar += ""+numNCharsPairedEnd+"M";
							//Add skip ref if needed
							if (nextRef<noSoftClipAln2First) cigar+=(noSoftClipAln2First-nextRef)+"N";
							//Build last part of the CIGAR
							String cigar2 = aln2.getCigarString();
							if(internalSoftClip2==0) cigar+= cigar2;
							else {
								char [] cigar2Chr = cigar2.toCharArray();
								int skipPos = cigar2.indexOf('S');
								cigar2Chr[skipPos] = 'M';
								cigar+=new String(cigar2Chr);
							}
						} else {
							//Skip second read
							last = last1;
							cigar+=""+(seq2.length()+numNCharsPairedEnd)+"S";
						}
						//System.out.println("Combined CIGAR: "+cigar);
						ReadAlignment combined = new ReadAlignment(aln1.getSequenceIndex(), aln1.getFirst(), last, seq.length(), 0);
						combined.setReadName(consensusName);
						combined.setSequenceName(aln1.getSequenceName());
						combined.setReadCharacters(seq);
						combined.setCigarString(cigar);
						combined.setAlignmentQuality(aln1.getAlignmentQuality());
						
						if(consensusName.equals(debugConsensusName)) System.err.println("Combined alignment (aln1 < aln2): " + combined.getReadCharacters() + "\t" +  combined.getFirst() + "\t" + combined.getCigarString());
						if(pairs.size()>1) combined.setAlignmentQuality((byte)10);
						alignments.put(consensusName,combined);
					} else if (aln2.getFirst() < aln1.getLast()) {
				 		int internalSoftClip2 = aln2.getSoftClipEnd();
				 		int internalSoftClip1 = aln1.getSoftClipStart();
						String cigarNonOverlap = aln2.getCigarString();
						if(internalSoftClip2>0) {
							//Transform into M
							cigarNonOverlap = cigarNonOverlap.substring(0,cigarNonOverlap.length()-1);
							cigarNonOverlap+="M";
						}
						//Add N matches
						cigarNonOverlap += ""+numNCharsPairedEnd+"M";
						int nextRef = last2+internalSoftClip2+numNCharsPairedEnd+1;
						int noSoftClipAln1First = aln1.getFirst()-internalSoftClip1;
						
						int firstPosCombined;
						String cigar;
						
						//Build last part of the CIGAR
						String cigar1 = aln1.getCigarString();
						if(nextRef<=noSoftClipAln1First) {
							firstPosCombined = aln2.getFirst();
							cigar = cigarNonOverlap;
							//Add skip ref
							if (nextRef<noSoftClipAln1First) cigar+=(noSoftClipAln1First-nextRef)+"N";
						} else {
							//Skip second read
							firstPosCombined = aln1.getFirst();
							cigar=""+(seq2.length()+numNCharsPairedEnd)+"S";
						}
						if(internalSoftClip1==0 || firstPosCombined==aln1.getFirst()) cigar+= cigar1;
						else {
							char [] cigar1Chr = cigar1.toCharArray();
							int skipPos = cigar1.indexOf('S');
							cigar1Chr[skipPos] = 'M';
							cigar+=new String(cigar1Chr);
						}
						//System.out.println("Combined CIGAR: "+cigar);
						ReadAlignment combined = new ReadAlignment(aln1.getSequenceIndex(), firstPosCombined, aln1.getLast(), seq.length(), 16);
						combined.setReadName(consensusName);
						combined.setSequenceName(aln1.getSequenceName());
						combined.setReadCharacters(DNAMaskedSequence.getReverseComplement(consensus.getCharacters()));
						combined.setCigarString(cigar);
						combined.setAlignmentQuality(aln1.getAlignmentQuality());
						if(consensusName.equals(debugConsensusName)) System.err.println("Combined alignment aln2 < aln1): " + combined.getReadCharacters() + "\t" +  combined.getFirst() + "\t" + combined.getCigarString());
						if(pairs.size()>1) combined.setAlignmentQuality((byte)10);
						alignments.put(consensusName,combined);
					} else {
						oddAlign++;
						unmappedReadPaired++;
						unmappedRead++;
					}
					if(consensusName.equals(debugConsensusName)) System.err.println();
				}
			}
		}
		return alignments;
	}
	
	public void saveAlignments(Map<String, ReadAlignment> alignments, String outputFile) throws IOException {	
		try (PrintStream outAlns = new PrintStream(outputFile);
			 ReadAlignmentFileWriter writer = new ReadAlignmentFileWriter(genome.getSequencesMetadata(), outAlns)) {
			for(ReadAlignment aln:alignments.values()) {
				writer.write(aln);
			}
		}
	}
	
	/**
	 * This method builds a hash with the provided alignments to be consulted at translation time
	 * @throws IOException
	 */
	public Map<String, ReadAlignment> loadAlignmentsFromBAM(String bamFilename) throws IOException {
		Map<String, ReadAlignment> alignments = new HashMap<>();
		try (ReadAlignmentFileReader bamOpenFile = new ReadAlignmentFileReader(bamFilename)) {
			Iterator<ReadAlignment> bamReader = bamOpenFile.iterator();
			while(bamReader.hasNext()) {
				ReadAlignment algn = bamReader.next();
				String algnName = "";
				if(algn.getReadName().contains("_")) {
					algnName = algn.getReadName().split("_")[1];
				} else {
					algnName = algn.getReadName();
				}
				//int clusterId = Integer.parseInt(algnName);
				if(algn.isReadUnmapped()) {
					unmappedRead++;
				} else if(!algn.isSecondary()) {
					alignments.put(algnName, algn);
				}
			}
		}
		return alignments;
	}
	
}
