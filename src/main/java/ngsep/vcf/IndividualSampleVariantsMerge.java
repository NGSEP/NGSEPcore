/*******************************************************************************
 * NGSEP - Next Generation Sequencing Experience Platform
 * Copyright 2016 Jorge Duitama
 *
 * This file is part of NGSEP.
 *
 *     NGSEP is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     NGSEP is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with NGSEP.  If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/
package ngsep.vcf;

import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;
import java.util.logging.Logger;

import ngsep.genome.GenomicRegionComparator;
import ngsep.genome.GenomicRegionPositionComparator;
import ngsep.genome.GenomicRegionSortedCollection;
import ngsep.genome.GenomicRegionSpanComparator;
import ngsep.main.CommandsDescriptor;
import ngsep.main.ProgressNotifier;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;
import ngsep.sequences.io.SimpleSequenceListLoader;
import ngsep.variants.GenomicVariant;
import ngsep.variants.GenomicVariantImpl;
import ngsep.variants.SNV;

public class IndividualSampleVariantsMerge {
	
	// Constants for default values
	
	// Logging and progress
	private Logger log = Logger.getLogger(IndividualSampleVariantsMerge.class.getName());
	private ProgressNotifier progressNotifier=null;
	
	// Parameters
	private String sequenceNamesFile = null;
	private String outputFile = null;
	private List<String> vcfFiles = new ArrayList<String>();
	
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
	
	
	public String getSequenceNamesFile() {
		return sequenceNamesFile;
	}
	public void setSequenceNamesFile(String sequenceNamesFile) {
		this.sequenceNamesFile = sequenceNamesFile;
	}
	
	public String getOutputFile() {
		return outputFile;
	}
	public void setOutputFile(String outputFile) {
		this.outputFile = outputFile;
	}
	
	public static void main(String[] args) throws Exception {
		IndividualSampleVariantsMerge instance = new IndividualSampleVariantsMerge();
		int i = CommandsDescriptor.getInstance().loadOptions(instance, args);
		instance.vcfFiles = new ArrayList<String>();
		for(;i<args.length;i++) {
			instance.vcfFiles.add(args[i]);
		}
		instance.run();
	}
	
	public void run () throws IOException {
		if (sequenceNamesFile==null || sequenceNamesFile.isEmpty()) throw new IOException("The file with sequence names is a required parameter");
		if (outputFile==null || outputFile.isEmpty()) throw new IOException("The output file is a required parameter");
		SimpleSequenceListLoader seqNameHandler = new SimpleSequenceListLoader();
		QualifiedSequenceList sequenceNames = seqNameHandler.loadSequences(sequenceNamesFile);
		log.info("Loaded "+sequenceNames.size()+" sequence names from file "+sequenceNamesFile);
		GenomicRegionSortedCollection<GenomicVariant> variants = mergeVariants(vcfFiles,sequenceNames);
		log.info("Merged variants. Total "+variants.size());
		printVariants(outputFile, variants);
	}

	public GenomicRegionSortedCollection<GenomicVariant> mergeVariants(List<String> vcfFiles, QualifiedSequenceList sequenceNames) throws IOException {
		log.info("Loading and merging variants");
		GenomicRegionComparator comparator = new GenomicRegionComparator(sequenceNames);
		GenomicRegionSortedCollection<GenomicVariant> variants = new GenomicRegionSortedCollection<GenomicVariant>(sequenceNames);
		List<GenomicRegionSortedCollection<GenomicVariant>> consolidatedPerSample = new ArrayList<GenomicRegionSortedCollection<GenomicVariant>>();
		int n=0;
		for(String f:vcfFiles) {
			log.info("Merging variants from file: "+f); 
			List <GenomicVariant> sampleVariants = VCFFileReader.loadVariants(f,true);
			log.info("Loaded "+sampleVariants.size()+" variants");
			consolidatedPerSample.add(mergeVariants (variants,sampleVariants,comparator));
			log.info("Merged variants from file: "+f+". Total variants: "+variants.size());
			n++;
			if (progressNotifier!=null) {
				if (!progressNotifier.keepRunning(n)) return new GenomicRegionSortedCollection<GenomicVariant>();
			}
		}
		log.info("Merging variants in the same location");
		variants = mergeOverlappingVariants(variants,consolidatedPerSample);
		log.info("Variants after last merging: "+variants.size());
		return variants;
	}
	
	private GenomicRegionSortedCollection<GenomicVariant> mergeVariants(GenomicRegionSortedCollection<GenomicVariant> variants,List <GenomicVariant> sampleVariants, GenomicRegionComparator comparator ) throws IOException {
		List<GenomicVariant> answerList = new ArrayList<GenomicVariant>();
		List<GenomicVariant> variantsList = variants.asList();
		int i=0;
		int j=0;
		while(j<sampleVariants.size()) {
			GenomicVariant c2 = sampleVariants.get(j);
			if(i<variantsList.size()) {
				GenomicVariant c1 = variantsList.get(i);
				int cmp = comparator.compare(c1, c2);
				if(cmp<0) {
					i++;
				} else if (cmp==0) {
					boolean add = true;
					for(int k=i;k<variantsList.size();k++) {
						GenomicVariant c3 = variantsList.get(k);
						int cmp2 = comparator.compare(c3, c2);
						if(cmp2==0 && c3.isCompatible(c2)) {
							add = false;
							if(c3.getVariantQS()<c2.getVariantQS()) c3.setVariantQS(c2.getVariantQS());
							answerList.add(c3);
							break;
						} else if (cmp2!=0) {
							break;
						}
					}
					if(add) {
						variants.add(c2);
						answerList.add(c2);
					}
					i++;
					j++;
				} else {
					variants.add(c2);
					answerList.add(c2);
					j++;
				}
			} else {
				variants.add(c2);
				answerList.add(c2);
				j++;
			}
		}
		GenomicRegionSortedCollection<GenomicVariant> answer = new GenomicRegionSortedCollection<GenomicVariant>(variants.getSequenceNames());
		answer.addAll(answerList);
		return answer;
	}

	private GenomicRegionSortedCollection<GenomicVariant> mergeOverlappingVariants(GenomicRegionSortedCollection<GenomicVariant> variants,List<GenomicRegionSortedCollection<GenomicVariant>> consolidatedPerSample) {
		QualifiedSequenceList names = variants.getSequenceNames();
		GenomicRegionSortedCollection<GenomicVariant> answer = new GenomicRegionSortedCollection<GenomicVariant>(names);
		for(QualifiedSequence seq:names) {
			GenomicRegionSortedCollection<GenomicVariant> seqVariants = variants.getSequenceRegions(seq.getName());
			List<GenomicVariant> overlappingVariants = new ArrayList<GenomicVariant>();
			int first = -1;
			int last = -1;
			for(GenomicVariant variant:seqVariants) {
				//System.out.println("Next Variant "+variant.getSequenceName()+":"+variant.getFirst()+" "+variant.getReference()+" "+variant.getAlleles()[1]);
				if(overlappingVariants.size() > 0 && !GenomicRegionSpanComparator.getInstance().span(variant, first, last)) {
					answer.addAll(mergeOverlappingVariants(overlappingVariants,first,last,consolidatedPerSample));
					overlappingVariants.clear();
				}
				overlappingVariants.add(variant);
				if(overlappingVariants.size()==1) {
					first = variant.getFirst();
					last = variant.getLast();
				} else if(last<variant.getLast()) last =variant.getLast();
			}
			if(overlappingVariants.size() > 0) {
				answer.addAll(mergeOverlappingVariants(overlappingVariants,first,last,consolidatedPerSample));
			}
		}
		return answer;
	}
	//PRE: list is not empty and only have overlapping variants in the same sequence; first is the start of the first variant
	private List<GenomicVariant> mergeOverlappingVariants(List<GenomicVariant> overlappingVariants, int first, int last, List<GenomicRegionSortedCollection<GenomicVariant>> consolidatedPerSample ) {
		if(overlappingVariants.size()<=1) return overlappingVariants;
		log.info("Merging "+overlappingVariants.size()+" overlapping variants at "+overlappingVariants.get(0).getSequenceName()+": "+first+"-"+last);
		List<GenomicVariant> answer = new ArrayList<GenomicVariant>();
		answer.addAll(extractSNVs(overlappingVariants));
		if(first < last || containsIndels(overlappingVariants)) answer.add(buildMultiallelicIndel(overlappingVariants,first, last, consolidatedPerSample));
		Collections.sort(answer,GenomicRegionPositionComparator.getInstance());
		return answer;
	}

	

	//PRE: list is not empty and only have overlapping variants in the same sequence;
	private List<GenomicVariant> extractSNVs(List<GenomicVariant> overlappingVariants) {
		List<GenomicVariant> snvs = new ArrayList<GenomicVariant>();
		List<GenomicVariant> overlappingSNVs = new ArrayList<GenomicVariant>();
		int pos = -1;
		for(GenomicVariant var:overlappingVariants) {
			if(!isSNV(var)) continue;
			if(overlappingSNVs.size()>0 && var.getFirst()!=pos) {
				snvs.add(mergeOverlappingSNVs(overlappingSNVs));
				overlappingSNVs.clear();
			}
			overlappingSNVs.add(var);
			pos = var.getFirst();
		}
		if(overlappingSNVs.size()>0) {
			snvs.add(mergeOverlappingSNVs(overlappingSNVs));
		}
		return snvs;
	}

	private boolean isSNV(GenomicVariant var) {
		if(var instanceof SNV) return true;
		String [] alleles = var.getAlleles();
		for(int i=0;i<alleles.length;i++) {
			if(alleles[i].length()>1) return false;
		}
		return true;
	}
	
	//PRE: list is not empty and only have overlapping SNVs in the same position;
	private GenomicVariant mergeOverlappingSNVs(List<GenomicVariant> overlappingSNVs) {
		GenomicVariant firstSNV = overlappingSNVs.get(0);
		if(overlappingSNVs.size()==1) return firstSNV;
		log.info("Creating multiallelic SNP");
		Set<String> sortedAlleles = new TreeSet<String>(); 
		String reference = firstSNV.getReference();
		short maxQS = 0;
		byte type = firstSNV.getType();
		for(GenomicVariant var:overlappingSNVs) {
			sortedAlleles.addAll(Arrays.asList(var.getAlleles()));
			if(maxQS<var.getVariantQS()) maxQS = var.getVariantQS();
			type = GenomicVariantImpl.mergeType(type, var.getType());
		}
		List<String> alleles = new ArrayList<String>();
		alleles.add(reference);
		sortedAlleles.remove(reference);
		alleles.addAll(sortedAlleles);
		log.info("Alleles: "+alleles);
		GenomicVariantImpl multiallelicVar = new GenomicVariantImpl(firstSNV.getSequenceName(), firstSNV.getFirst(), alleles);
		multiallelicVar.setVariantQS(maxQS);
		multiallelicVar.setType(type);
		return multiallelicVar;
	}

	private boolean containsIndels(List<GenomicVariant> overlappingVariants) {
		for(GenomicVariant var:overlappingVariants) {
			if(!isSNV(var)) return true;
		}
		return false;
	}
	
	//PRE: list is not empty and only have overlapping variants in the same sequence; first is the start of the first variant
	private GenomicVariant buildMultiallelicIndel(List<GenomicVariant> overlappingVariants, int first, int last, List<GenomicRegionSortedCollection<GenomicVariant>> consolidatedPerSample) {
		log.info("Creating multiallelic indel");
		String reference = buildReferenceFromVariants (overlappingVariants,first);
		String seqName = overlappingVariants.get(0).getSequenceName();
		Set<String> sortedAlleles = new TreeSet<String>(); 
		short maxQS = 0;
		byte type = GenomicVariant.TYPE_INDEL;
		for(GenomicVariant var:overlappingVariants) {
			if(maxQS<var.getVariantQS()) maxQS = var.getVariantQS();
			type = GenomicVariantImpl.mergeType(type, var.getType());
		}
		//Create one mutated allele per sample assuming that all alternative alleles appear at the same haplotype
		for(GenomicRegionSortedCollection<GenomicVariant> varsSample:consolidatedPerSample) {
			GenomicRegionSortedCollection<GenomicVariant> sampleVarsRegion = varsSample.findSpanningRegions(seqName, first, last);
			int nVarsReg = sampleVarsRegion.size(); 
			//log.info("Creating allele from "+sampleVarsRegion.size()+" sample variants");
			if(nVarsReg>1 ) {
				String allele = makeAlternativeAllele(reference, first, sampleVarsRegion);
				if(allele!=null) sortedAlleles.add(allele);
			} else if (nVarsReg==1) {
				GenomicVariant varSample = sampleVarsRegion.asList().get(0);
				String [] varAlleles = varSample.getAlleles();
				for(int i=1;i<varAlleles.length;i++) {
					String allele = makeAlternativeAllele(reference, first, varAlleles[i], varSample.getFirst(),varSample.getLast());
					if(allele!=null) sortedAlleles.add(allele);
					if(sortedAlleles.size()==GenomicVariant.MAX_NUM_ALLELES-2) break;
				}
				
			}
			if(sortedAlleles.size()==GenomicVariant.MAX_NUM_ALLELES-2) {
				log.severe("WARN:Maximum number of alleles reached for variant at "+seqName+": "+first+"-"+last);
				break;
			}
		}
		//Create variant with the consolidated alleles
		List<String> alleles = new ArrayList<String>();
		alleles.add(reference);
		sortedAlleles.remove(reference);
		alleles.addAll(sortedAlleles);
		log.info("Alleles: "+alleles);
		GenomicVariantImpl multiallelicVar = new GenomicVariantImpl(overlappingVariants.get(0).getSequenceName(), first, alleles);
		multiallelicVar.setVariantQS(maxQS);
		multiallelicVar.setType(type);
		return multiallelicVar;
	}

	private String makeAlternativeAllele(String reference, int first, String alternative, int firstAlt, int lastAlt) {
		StringBuilder answer = new StringBuilder();
		int relEnd = firstAlt - first;
		answer.append(reference.substring(0,relEnd));
		answer.append(alternative);
		int relStart = lastAlt - first + 1;
		answer.append(reference.substring(relStart));
		return answer.toString();
	}

	private String makeAlternativeAllele(String reference, int first, GenomicRegionSortedCollection<GenomicVariant> variants) {
		if(variants==null || variants.size()==0) return null;
		StringBuilder answer = new StringBuilder();
		int next = first;
		for(GenomicVariant variant:variants) {
			if(variant.getLast()<next) continue;
			int cut = next - variant.getFirst();
			if (cut>0) log.severe("WARN: Partial allele will be printed for variant at "+variant.getSequenceName()+": "+variant.getFirst()+" Cut: "+cut);
			else if(cut < 0) {
				int relStart = next - first;
				int relEnd = variant.getFirst() - first;
				answer.append(reference.substring(relStart,relEnd));
			}
			String [] alleles = variant.getAlleles();
			if(alleles.length>1) {
				String allele = alleles[1];
				if(cut < allele.length()) {
					if(cut > 0) allele = allele.substring(cut);
					answer.append(allele);
				}
			}
			next = variant.getLast() + 1;
		}
		int relNext = next-first;
		if(relNext<reference.length()) answer.append(reference.substring(relNext));
		return answer.toString();
	}

	private String buildReferenceFromVariants(List<GenomicVariant> overlappingVariants,int first) {
		StringBuilder refAllele = new StringBuilder(overlappingVariants.get(0).getReference());
		for(GenomicVariant var:overlappingVariants) {
			int relFirst = var.getFirst() - first;
			String refVar = var.getReference();
			if(relFirst+refVar.length()>refAllele.length()) {
				refAllele.append(refVar.substring(refAllele.length()-relFirst));
			}
		}
		
		return refAllele.toString();
	}
	
	public void printVariants(String outFile, GenomicRegionSortedCollection<GenomicVariant> variants) throws IOException {
		VCFFileWriter writer = new VCFFileWriter();
		
		try (PrintStream out = new PrintStream(outFile)){
			VCFFileHeader header = VCFFileHeader.makeDefaultEmptyHeader();
			writer.printHeader(header,out);
			for(GenomicVariant variant:variants) {
				VCFRecord record = new VCFRecord(variant, VCFRecord.DEF_FORMAT_ARRAY_NONE, new ArrayList<>(), header);
				writer.printVCFRecord(record, out);
			}
		}
	}

}