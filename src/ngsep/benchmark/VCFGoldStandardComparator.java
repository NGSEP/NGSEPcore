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
package ngsep.benchmark;
import java.io.IOException;
import java.io.PrintStream;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import ngsep.genome.GenomicRegion;
import ngsep.genome.GenomicRegionComparator;
import ngsep.genome.ReferenceGenome;
import ngsep.genome.io.SimpleGenomicRegionFileHandler;
import ngsep.main.CommandsDescriptor;
import ngsep.math.LogMath;
import ngsep.math.PhredScoreHelper;
import ngsep.sequences.QualifiedSequenceList;
import ngsep.variants.CalledGenomicVariant;
import ngsep.variants.CalledGenomicVariantImpl;
import ngsep.variants.GenomicVariant;
import ngsep.variants.GenomicVariantImpl;
import ngsep.variants.VariantCallReport;
import ngsep.vcf.VCFFileReader;
import ngsep.vcf.VCFFileWriter;
import ngsep.vcf.VCFRecord;

/**
 * Compares a gold standard VCF or gVCF with a VCF file with predicted calls for the same sample  
 * @author Jorge Duitama
 *
 */
public class VCFGoldStandardComparator {
	private ReferenceGenome genome;
	private Map<String, List<GenomicRegion>> complexRegions;
	private Map<String, List<GenomicRegion>> confidenceRegions;
	private String outVCF = null;
	private int mode = 0;
	private short minQuality = 0;
	private int posPrint = -1;
	private boolean genomicVCF = true;
	
	
	private Map<Byte, GoldStandardComparisonCounts> countsPerType = new HashMap<>();
	private long confidenceRegionsLength = 0;

	VCFFileWriter writer = new VCFFileWriter();
	public static void main(String[] args) throws Exception {
		VCFGoldStandardComparator instance = new VCFGoldStandardComparator();
		int i = CommandsDescriptor.getInstance().loadOptions(instance, args);
		String referenceFile = args[i++];
		
		
		String file1 = args[i++];
		String file2 = args[i++];
		
		if(args.length>=i+3) {
			instance.outVCF = args[i++];
			instance.mode = Integer.parseInt(args[i++]);
			instance.minQuality = Short.parseShort(args[i++]);
		}
		
		instance.genome = new ReferenceGenome(referenceFile);
		instance.compareFiles(file1,file2);
		instance.printStatistics (System.out);
	}

	/**
	 * @return the genome
	 */
	public ReferenceGenome getGenome() {
		return genome;
	}
	
	/**
	 * @param genome the genome to set
	 */
	public void setGenome(ReferenceGenome genome) {
		this.genome = genome;
	}
	
	

	/**
	 * @return the genomicVCF
	 */
	public boolean isGenomicVCF() {
		return genomicVCF;
	}

	/**
	 * @param genomicVCF the genomicVCF to set
	 */
	public void setGenomicVCF(boolean genomicVCF) {
		this.genomicVCF = genomicVCF;
	}

	public void setGenomicVCF(Boolean genomicVCF) {
		this.setGenomicVCF(genomicVCF.booleanValue());
	}
	/**
	 * @return the complexRegions of a sequence
	 */
	public List<GenomicRegion> getComplexRegions(String sequenceName) {
		return complexRegions.get(sequenceName);
	}

	/**
	 * @param complexRegionsFile File wit the complexRegions to set
	 */
	public void setComplexRegions(String complexRegionsFile) throws IOException {
		SimpleGenomicRegionFileHandler handler = new SimpleGenomicRegionFileHandler();
		this.complexRegions = handler.loadRegionsAsMap(complexRegionsFile);
	}

	/**
	 * @return the confidenceRegions
	 */
	public List<GenomicRegion> getConfidenceRegions(String sequenceName) {
		return confidenceRegions.get(sequenceName);
	}

	/**
	 * @param confidenceRegions the confidenceRegions to set
	 */
	public void setConfidenceRegions(String regionsFile) throws IOException {
		SimpleGenomicRegionFileHandler handler = new SimpleGenomicRegionFileHandler();
		this.confidenceRegions = handler.loadRegionsAsMap(regionsFile);
		confidenceRegionsLength = 0;
		for(List<GenomicRegion> regions: confidenceRegions.values()) {
			for(GenomicRegion g:regions) {
				confidenceRegionsLength+=g.length();
			}
		}
	}
	
	private void loadConfidenceRegionsFromVCF(String filename, QualifiedSequenceList sequenceNames) throws IOException {
		confidenceRegions = new HashMap<>();
		confidenceRegionsLength = 0;
		try (VCFFileReader in = new VCFFileReader(filename)) {
			in.setLoadMode(VCFFileReader.LOAD_MODE_MINIMAL);
			in.setSequences(sequenceNames);
			Iterator<VCFRecord> it = in.iterator();
			while(it.hasNext()) {
				VCFRecord record = it.next();
				CalledGenomicVariant call = record.getCalls().get(0); 
				if(call.isHomozygousReference()) {
					List<GenomicRegion> regionsSeq = confidenceRegions.get(call.getSequenceName());
					if(regionsSeq==null) {
						regionsSeq = new ArrayList<>();
						confidenceRegions.put(call.getSequenceName(), regionsSeq);
					}
					regionsSeq.add(record.getVariant());
					confidenceRegionsLength+=record.length();
				}
			}
		}
		
	}

	public void compareFiles(String vcfGS, String vcfTest) throws IOException {
		QualifiedSequenceList sequenceNames = genome.getSequencesMetadata();
		if(confidenceRegions==null && genomicVCF) loadConfidenceRegionsFromVCF(vcfGS,sequenceNames);
		int nSeqs = sequenceNames.size();
		countsPerType.put(GenomicVariant.TYPE_BIALLELIC_SNV, new GoldStandardComparisonCounts());
		countsPerType.put(GenomicVariant.TYPE_INDEL, new GoldStandardComparisonCounts());
		countsPerType.put(GenomicVariant.TYPE_STR, new GoldStandardComparisonCounts());
		try (VCFFileReader inGS = new VCFFileReader(vcfGS);
			 VCFFileReader inTest = new VCFFileReader(vcfTest)) {
			inGS.setLoadMode(VCFFileReader.LOAD_MODE_MINIMAL);
			inGS.setSequences(sequenceNames);
			inTest.setSequences(sequenceNames);
			LinkedList<GenomicRegion> confidenceRegionsSeq=null;
			List<GenomicRegion> complexRegionsSeq=null;
			int pIdx = 0;
			List<CalledGenomicVariant> gsCalls = new ArrayList<>();
			List<CalledGenomicVariant> testCalls = new ArrayList<>();
			int sequenceIdx = -1;
			int clusterFirst = 0;
			int clusterLast = 0;
			byte clusterType = GenomicVariant.TYPE_UNDETERMINED;
			
			Iterator<VCFRecord> itGS = inGS.iterator();
			VCFRecord recordGS = loadNextRecord(itGS, false);
			
			
			Iterator<VCFRecord> itTest = inTest.iterator();
			VCFRecord recordTest = recordGS = loadNextRecord(itTest, false);
			
			while(true) {
				int nextSeqIdxGS = nSeqs;
				if(recordGS!=null) nextSeqIdxGS = sequenceNames.indexOf(recordGS.getSequenceName());
				int nextSeqIdxTest = nSeqs;
				if(recordTest!=null) nextSeqIdxTest = sequenceNames.indexOf(recordTest.getSequenceName());
				int nextSeqIdx = Math.min(nextSeqIdxGS, nextSeqIdxTest);
				if(nextSeqIdx == nSeqs) {
					processClusterCalls (gsCalls, testCalls, clusterFirst, sequenceNames.get(sequenceIdx).getLength(), clusterType, confidenceRegionsSeq );
					break;
				}
				
				if(nextSeqIdx>sequenceIdx) {
					if(sequenceIdx>=0) processClusterCalls (gsCalls, testCalls, clusterFirst, sequenceNames.get(sequenceIdx).getLength(), clusterType, confidenceRegionsSeq);
					sequenceIdx = nextSeqIdx;
					String sequenceName = sequenceNames.get(sequenceIdx).getName();
					if(confidenceRegions!=null) confidenceRegionsSeq = new LinkedList<>(confidenceRegions.get(sequenceName));
					if(complexRegions!=null) complexRegionsSeq = complexRegions.get(sequenceName);
					pIdx = 0;
					gsCalls.clear();
					testCalls.clear();
					clusterFirst=clusterLast=0;
					clusterType = GenomicVariant.TYPE_UNDETERMINED;
				}
				if(complexRegionsSeq!=null) {
					for(;pIdx<complexRegionsSeq.size();pIdx++) {
						GenomicRegion region = complexRegionsSeq.get(pIdx);
						if(clusterLast+10<region.getFirst()) break;
						if(clusterFirst<=region.getLast() && clusterLast>=region.getFirst()) {
							clusterType = GenomicVariant.TYPE_STR;
							clusterLast=Math.max(clusterLast, region.getLast());
						}
					}
				}
				int nextClusterFirst = clusterLast+Math.max(10, clusterLast-clusterFirst+1);
				nextClusterFirst = Math.min(nextClusterFirst, clusterLast+100);
				boolean gsClose = nextSeqIdxGS == sequenceIdx && nextClusterFirst>recordGS.getFirst();
				boolean testClose = nextSeqIdxTest == sequenceIdx && nextClusterFirst>recordTest.getFirst();
				if (!gsClose && !testClose) {
					processClusterCalls (gsCalls, testCalls, clusterFirst, clusterLast, clusterType, confidenceRegionsSeq);
					gsCalls.clear();
					testCalls.clear();
					clusterFirst = 0;
					clusterType = GenomicVariant.TYPE_UNDETERMINED;
					gsClose = nextSeqIdxGS == sequenceIdx;
					testClose = nextSeqIdxTest == sequenceIdx;
				}
				if(gsClose) {
					CalledGenomicVariant callGS = recordGS.getCalls().get(0);
					gsCalls.add(callGS);
					if(clusterType==GenomicVariant.TYPE_UNDETERMINED) clusterType = loadType(callGS);
					//Default type for cluster calls
					if(gsCalls.size()>1) clusterType= GenomicVariant.TYPE_STR;
					
					if(clusterFirst == 0 ) clusterFirst = recordGS.getFirst();
					else clusterFirst = Math.min(clusterFirst, recordGS.getFirst());
					clusterLast = Math.max(clusterLast, recordGS.getLast());
					recordGS = loadNextRecord(itGS, false);
				}
				if(testClose) {
					testCalls.add(recordTest.getCalls().get(0));
					if(clusterFirst == 0 ) clusterFirst = recordTest.getFirst();
					else clusterFirst = Math.min(clusterFirst, recordTest.getFirst());
					clusterLast = Math.max(clusterLast, recordTest.getLast());
					recordTest = loadNextRecord(itTest, false);
				}
			}
		}
	}
	
	
	

	private VCFRecord loadNextRecord(Iterator<VCFRecord> it, boolean loadHomozygousReference) {
		VCFRecord record = null;
		while(it.hasNext()) {
			VCFRecord record2 = it.next();
			CalledGenomicVariant call = record2.getCalls().get(0); 
			if (call.isUndecided()) continue;
			if(!loadHomozygousReference && call.isHomozygousReference()) continue;
			record = record2;
			break;
		}
		return record;
	}

	private void processClusterCalls(List<CalledGenomicVariant> gsCalls, List<CalledGenomicVariant> testCalls, int clusterFirst, int clusterLast, byte clusterGSType, LinkedList<GenomicRegion> confidenceRegionsSeq) {
		int lastRowCounts = GoldStandardComparisonCounts.NUM_ROWS_COUNTS-1;
		if(gsCalls.size()==0 && testCalls.size()==0) return;
		//if(gsCalls.size()>0) System.out.println("Calls: "+gsCalls.size()+" first: "+clusterFirst+" first gs call: "+gsCalls.get(0).getFirst());
		while(confidenceRegionsSeq!=null && confidenceRegionsSeq.size()>0) {
			GenomicRegion nextConfident = confidenceRegionsSeq.peek();
			if(nextConfident.getLast()<clusterFirst-10) confidenceRegionsSeq.removeFirst();
			else break;
		}
		if(gsCalls.size()==0) {
			processPossibleFalsePositives(testCalls,confidenceRegionsSeq, lastRowCounts);
		}
		else if(testCalls.size()==0) {
			processFalseNegatives(gsCalls, lastRowCounts);
		} else {
			CalledGenomicVariant firstGS = gsCalls.get(0);
			CalledGenomicVariant firstTest = testCalls.get(0);
			int genotypeFirstGS = getGenotypeNumber(firstGS);
			int genotypeFirstTest = getGenotypeNumber(firstTest);
			short qualFirstTest = loadGenotypeQuality(firstTest);
			GoldStandardComparisonCounts counts = countsPerType.get(clusterGSType);
			if(clusterGSType == GenomicVariant.TYPE_BIALLELIC_SNV) {
				//Process SNVs
				if(firstGS.getFirst()==firstTest.getFirst() && firstGS.getAlleles()[1].equals(firstTest.getAlleles()[1])) {
					// Match between gold standard and test SNV
					int k = 3*genotypeFirstGS+genotypeFirstTest;
					int row = Math.min(qualFirstTest/10, lastRowCounts); 
					counts.update(0,row,k);
					counts.update(row+1,lastRowCounts,9+genotypeFirstGS);
				} else {
					//Isolated SNV calls in different close positions or with different alternative alleles
					processPossibleFalsePositives(testCalls,confidenceRegionsSeq, lastRowCounts);
					processFalseNegatives(gsCalls, lastRowCounts);
				}
			} else {
				GoldStandardHaplotypeReconstruction gsHaps = new GoldStandardHaplotypeReconstruction(genome.getReference(firstGS.getSequenceName(), clusterFirst-5, clusterLast).toString(), gsCalls, clusterFirst-5);
				String [] gsHaplotypes = gsHaps.getPhasedAlleles();
				int genotypeGS = gsHaps.getGenotypeNumber();
				String [] matchingHaplotypes = gsHaps.buildMatchingHaplotypes(testCalls);
				int genotypeTest = CalledGenomicVariant.GENOTYPE_HOMOALT;
				if(!matchingHaplotypes[0].equals(matchingHaplotypes[1])) genotypeTest = CalledGenomicVariant.GENOTYPE_HETERO;
				if(genotypeGS==genotypeTest) {
					//Record genotype errors if haplotype sequences do not match even if the genotypes match 
					if(genotypeGS==CalledGenomicVariant.GENOTYPE_HOMOALT && !gsHaplotypes[0].equals(matchingHaplotypes[0])) {
						genotypeTest = CalledGenomicVariant.GENOTYPE_HETERO;
						//System.out.println("Changing genotype test to heterozygous due to error in sequence. GS hap: "+gsHaplotypes[0]+" test hap: "+matchingHaplotypes[0]);
					} else if (genotypeGS==CalledGenomicVariant.GENOTYPE_HETERO && (!gsHaplotypes[0].equals(matchingHaplotypes[0]) || (!gsHaplotypes[1].equals(matchingHaplotypes[1])))) {
						genotypeTest = CalledGenomicVariant.GENOTYPE_HOMOALT;
					}
				}
				short qualTest = calculateQuality(testCalls);
				int k = 3*genotypeGS+genotypeTest;
				int row = Math.min(qualTest/10, lastRowCounts); 
				counts.update(0,row,k);
				counts.update(row+1,lastRowCounts,9+genotypeGS);
				if(mode == 3 && genotypeGS!=genotypeTest &&  qualTest>=minQuality) {
					System.out.println("Variant "+gsCalls.get(0).getSequenceName()+": "+gsHaps.getFirst());
					System.out.println(gsHaplotypes[0]);
					System.out.println(gsHaplotypes[1]);
					//writer.printVCFRecord(new VCFRecord(gsHaps, VCFRecord.DEF_FORMAT_ARRAY_MINIMAL, gsHaps, null), System.out);
					System.out.println(matchingHaplotypes[0]);
					System.out.println(matchingHaplotypes[1]);
				}
				
			}
		}
		while(confidenceRegionsSeq!=null && confidenceRegionsSeq.size()>0) {
			GenomicRegion nextConfident = confidenceRegionsSeq.peek();
			if(nextConfident.getLast()<clusterLast) confidenceRegionsSeq.removeFirst();
			else break;
		}
		
		
	}

	

	private void processPossibleFalsePositives(List<CalledGenomicVariant> testCalls, LinkedList<GenomicRegion> confidenceRegionsSeq, int lastRowCounts) {
		List<GenomicRegion> confidenceRegionsList = new ArrayList<>();
		if(confidenceRegionsSeq!=null) confidenceRegionsList.addAll(confidenceRegionsSeq);
		for(CalledGenomicVariant call:testCalls) {
			int genotypeTest = getGenotypeNumber(call);
			short qualTest = loadType(call);
			byte typeTest = loadType(call);
			int n = 12;
			for(GenomicRegion r: confidenceRegionsList) {
				if(r.getFirst()<=call.getFirst() && r.getLast()>=call.getLast()) {
					n=0;
					break;
				}
			}
			
			countsPerType.get(typeTest).update(0,Math.min(qualTest/10, lastRowCounts),n+genotypeTest);
		}
	}

	private void processFalseNegatives(List<CalledGenomicVariant> gsCalls, int lastRowCounts) {
		for(CalledGenomicVariant call:gsCalls) {
			int genotypeGS = getGenotypeNumber(call);
			byte type = loadType(call);
			countsPerType.get(type).update(0,lastRowCounts,9+genotypeGS);
		}
	}
	
	private void printStatistics(PrintStream out) {
		double confidentMbp = confidenceRegionsLength/1000000;
		out.println("Confident MBP: "+confidentMbp);
		for(GoldStandardComparisonCounts counts:countsPerType.values()) counts.setConfidentMbp(confidentMbp);
		out.println("SNVs");
		countsPerType.get(GenomicVariant.TYPE_BIALLELIC_SNV).print(System.out);
		out.println("Indels");
		countsPerType.get(GenomicVariant.TYPE_INDEL).print(System.out);
		out.println("STRs/OTHER");
		countsPerType.get(GenomicVariant.TYPE_STR).print(System.out);
		
	}

	public void compareFilesOld(String vcfGS, String vcfTest) throws IOException {
		QualifiedSequenceList sequenceNames = genome.getSequencesMetadata();
		GenomicRegionComparator compRegion = new GenomicRegionComparator(sequenceNames);
		
		countsPerType.put(GenomicVariant.TYPE_BIALLELIC_SNV, new GoldStandardComparisonCounts());
		countsPerType.put(GenomicVariant.TYPE_INDEL, new GoldStandardComparisonCounts());
		countsPerType.put(GenomicVariant.TYPE_STR, new GoldStandardComparisonCounts());
		//countsPerType.put(GenomicVariant.TYPE_UNDETERMINED, new GoldStandardComparisonCounts());
		List<CalledGenomicVariant> callsTest = VCFFileReader.loadCalledVariantsSingleIndividualVCF(vcfTest);
		int lastRowCounts = GoldStandardComparisonCounts.NUM_ROWS_COUNTS-1;
		PrintStream out = null;
		try (VCFFileReader inGS = new VCFFileReader(vcfGS)){ 
			QualifiedSequenceList seqNames = genome.getSequencesMetadata();
			inGS.setSequences(seqNames);
			inGS.setLoadMode(VCFFileReader.LOAD_MODE_MINIMAL);
			if(outVCF!=null) {
				out = new PrintStream(outVCF);
				writer.printHeader(inGS.getHeader(), out);
			}
			Iterator<VCFRecord> itGS = inGS.iterator();
			int idxTest = 0;
			while (itGS.hasNext()) {
				VCFRecord rGS = itGS.next();
				CalledGenomicVariant callGS = rGS.getCalls().get(0);
				if(!callGS.isUndecided())confidenceRegionsLength+=callGS.getReference().length();
				boolean referenceRegion = callGS.isHomozygousReference();
				byte typeGS = loadType(callGS);
				List<CalledGenomicVariant> callsIntersect = new ArrayList<>();
				while(idxTest<callsTest.size()) {
					CalledGenomicVariant callTest = callsTest.get(idxTest);
					if(callTest.isUndecided()) {
						idxTest++;
						continue;
					}
					short qualTest = loadGenotypeQuality(callTest);
					int lastBefore = (idxTest == 0 || !callsTest.get(idxTest-1).getSequenceName().equals(callTest.getSequenceName()))? 0:callsTest.get(idxTest-1).getLast();
					int firstAfter = (idxTest == callsTest.size()-1 || !callsTest.get(idxTest+1).getSequenceName().equals(callTest.getSequenceName()))? seqNames.get(callTest.getSequenceName()).getLength():callsTest.get(idxTest+1).getFirst();
					
					byte type = loadType(callTest);
					int column = getGenotypeNumber(callTest);
					CalledGenomicVariant expandedCallTest = expandReferenceIndels(callTest,lastBefore,firstAfter);
					int cmp = compRegion.compare(callGS, expandedCallTest);
					if(posPrint==callGS.getFirst()) System.out.println("Call GS: "+callGS.getFirst()+"-"+callGS.getLast()+" call test: "+callTest.getFirst()+"-"+callTest.getLast()+" expanded call: "+expandedCallTest.getFirst()+"-"+expandedCallTest.getLast()+"type: "+type+" column: "+column+" comparison: "+cmp);					
					if (cmp<-1) break;
					if (cmp<=1) {
						//Overlap between GS and extended region
						if(!referenceRegion) {
							callsIntersect.add(expandedCallTest);
							idxTest++;
							continue;
						}
						//Variant falling out of reference region will be processed by the next gs region
						int remainder = expandedCallTest.getLast() - callGS.getLast();	
						if(posPrint==callGS.getFirst()) System.out.println("Call GS: "+callGS.getFirst()+"-"+callGS.getLast()+" call test: "+callTest.getFirst()+"-"+callTest.getLast()+" remainder: "+remainder);
						if(!callTest.isHomozygousReference() && remainder>0) break;
						if(callTest.getFirst()<callGS.getFirst() || callGS.getLast()<callTest.getLast()) column+=12;
					} else column+=12;
					countsPerType.get(type).update(0,Math.min(qualTest/10, lastRowCounts),column);
					VCFRecord record = new VCFRecord(callTest, VCFRecord.DEF_FORMAT_ARRAY_QUALITY, callTest, null);
					if(mode == 2 && qualTest>=minQuality && column < 12 && !callTest.isHomozygousReference()) writer.printVCFRecord(record, out);
					idxTest++;
				}
				if(posPrint==callGS.getFirst()) System.out.println("GS call at "+callGS.getSequenceName()+":"+callGS.getFirst()+" type: "+typeGS+" undecided: "+callGS.isUndecided()+" Calls intersect: "+callsIntersect.size()+" next test call: "+callsTest.get(idxTest).getFirst());
				if(callGS.isUndecided() || referenceRegion ) {
					//Ignore variants that intersect with unknown gold standard sites
					continue;
				}
				//Process non-reference call in the gold standard
				int n1 = getGenotypeNumber(callGS);
				if(callsIntersect.size()==0) {
					//Gold standard variant not found
					countsPerType.get(typeGS).update(0,lastRowCounts,9+n1);
					if(mode == 1) writer.printVCFRecord(rGS, out);
					continue;
				}
				boolean covered = false;
				for(CalledGenomicVariant callTest:callsIntersect) {
					boolean consistent = isConsistent(callGS, callTest);
					byte typeTest = loadType(callTest);
					short qualTest = callTest.getGenotypeQuality();
					int n2 = getGenotypeNumber(callTest);
					VCFRecord rTest = new VCFRecord(callTest, VCFRecord.DEF_FORMAT_ARRAY_QUALITY, callTest, null);
					if(consistent) {
						GoldStandardComparisonCounts counts = countsPerType.get(typeGS);
						int k = 3*n1+n2;
						int row = Math.min(qualTest/10, lastRowCounts); 
						counts.update(0,row,k);
						counts.update(row+1,lastRowCounts,9+n1);
						if(mode == 1 &&  qualTest<minQuality) {
							writer.printVCFRecord(rGS, out);
							writer.printVCFRecord(rTest, out);
						}
						if(mode==3 && n1!=n2 &&  qualTest>=minQuality) {
							writer.printVCFRecord(rGS, out);
							writer.printVCFRecord(rTest, out);
						}
						covered= true;
					} else {
						//False positive for inconsistent variant
						countsPerType.get(typeTest).update(0,Math.min(qualTest/10, lastRowCounts),12+n2);
						if(mode==3 && qualTest>=minQuality) {
							writer.printVCFRecord(rGS, out);
							writer.printVCFRecord(rTest, out);
						}
					}
				}
				if(!covered) {
					GoldStandardComparisonCounts counts = countsPerType.get(typeGS);
					//False negative for GS not properly called
					counts.update(0,lastRowCounts,9+n1);	
				}
			}
			while(idxTest<callsTest.size()) {
				//Process test variants at the end
				CalledGenomicVariant callTest = callsTest.get(idxTest);
				if(!callTest.isUndecided()) {
					byte typeTest = loadType(callTest);
					short qualTest = loadGenotypeQuality(callTest);
					int n2 = getGenotypeNumber(callTest);
					countsPerType.get(typeTest).update(0,Math.min(qualTest/10, lastRowCounts),12+n2);
					VCFRecord rTest = new VCFRecord(callTest, VCFRecord.DEF_FORMAT_ARRAY_QUALITY, callTest, null);
					if(mode == 2 && qualTest>=minQuality && !callTest.isHomozygousReference()) writer.printVCFRecord(rTest, out);
				}
				idxTest++;
			}
		} finally {
			if(out!=null) {
				out.flush();
				out.close();
			}
		}
		printStatistics(System.out);
	}

	private CalledGenomicVariant expandReferenceIndels(CalledGenomicVariant callTest, int lastBefore, int firstAfter) {
		String [] alleles = callTest.getAlleles();
		if(alleles.length<2) return callTest;
		Set<Integer> lengths = new HashSet<>();
		for(int i=0;i<alleles.length;i++) lengths.add(alleles[i].length());
		if(lengths.size()==1 && lengths.iterator().next()==1) return callTest;
		//Add 1 reference bp at the end
		int refBefore = Math.max(callTest.getFirst()-5, lastBefore+2);
		
		CharSequence left = null;
		if(refBefore<callTest.getFirst()) left = genome.getReference(callTest.getSequenceName(), refBefore, callTest.getFirst()-1);
		int refAfter = Math.min(callTest.getLast()+5, firstAfter-2);
		CharSequence right = null;
		if(callTest.getLast()<refAfter) right = genome.getReference(callTest.getSequenceName(), callTest.getLast()+1,refAfter);
		if(left == null && right==null) return callTest;
		int first = callTest.getFirst();
		if(left!=null) first-=left.length();
		List<String> extendedAlleles = new ArrayList<>();
		for(int i=0;i<alleles.length;i++) {
			String a = alleles[i];
			if(left!=null) a= left+a;
			if(right!=null) a= a+right;
			extendedAlleles.add(a);
		}
		GenomicVariant newVariant = new GenomicVariantImpl(callTest.getSequenceName(), first, extendedAlleles);
		newVariant.setType(callTest.getType());
		newVariant.setVariantQS(callTest.getVariantQS());
		CalledGenomicVariantImpl answer = new CalledGenomicVariantImpl(newVariant, callTest.getIndexesCalledAlleles());
		answer.setGenotypeQuality(callTest.getGenotypeQuality());
		return answer;
	}

	public boolean isConsistent(CalledGenomicVariant callGS, CalledGenomicVariant callTest) {
		boolean consistent = true;
		List<String> allelesGS = buildExtendedAlleles(callGS, callTest.getFirst(), callTest.getLast());
		if(posPrint==callGS.getFirst()) System.out.println("AllelesGS "+allelesGS);
		
		
		String [] allelesTest = callTest.getAlleles();
		//Expected start of alternative alleles in reference alleles
		int offset = 0;
		if(callGS.getFirst()<callTest.getFirst()) offset = callTest.getFirst()-callGS.getFirst();
		String reference = allelesGS.get(0);
		if(offset>0 && offset<reference.length()) {
			reference = reference.substring(offset);
		}
		consistent = reference.startsWith(allelesTest[0]);
		if(!consistent) System.err.println("WARN: Inconsistent reference for comparison between "+callGS.getSequenceName()+":"+callGS.getFirst()+ " reference: "+reference+" and "+callTest.getSequenceName()+":"+callTest.getFirst()+" reference: "+allelesTest[0]+" offset: "+offset );
		for(int i=1;i<allelesTest.length && consistent;i++) {
			String alleleTest = allelesTest[i];
			if(posPrint==callGS.getFirst()) System.out.println("Next allele test "+alleleTest);
			boolean found = false;
			for(int j=1;j<allelesGS.size() && !found;j++) {
				String alleleGS = allelesGS.get(j);
				if(offset>0 && offset<alleleGS.length()) {
					alleleGS = alleleGS.substring(offset);
				}
				found = alleleGS.startsWith(alleleTest);
			}
			consistent = found;
		}
		return consistent;
	}
	
	private List<String> buildExtendedAlleles(GenomicVariant v, int firstTest, int lastTest) {
		String seqName = v.getSequenceName();
		String left = null;
		if(firstTest<v.getFirst()) left = genome.getReference(seqName, firstTest, v.getFirst()-1).toString();
		String right = null;
		if(v.getLast()<lastTest) right = genome.getReference(seqName, v.getLast()+1, lastTest).toString();
		List<String> answer = new ArrayList<>();
		String [] alleles = v.getAlleles();
		for(int i=0;i<alleles.length;i++) {
			String allele = alleles[i];
			String allO = "";
			if(left!=null) allO += left;
			allO+=allele;
			if(right!=null) allO += right;
			answer.add(allO.toUpperCase());
		}
		return answer;
	}
	
	private byte loadType(CalledGenomicVariant call) {
		byte type = call.getType();
		if(type == GenomicVariant.TYPE_BIALLELIC_SNV) return type;
		if(type == GenomicVariant.TYPE_INDEL) return type;
		if(type == GenomicVariant.TYPE_STR) return type;
		String [] alleles = call.getAlleles();
		if(alleles.length==2 && alleles[0].length()!=alleles[1].length()) return GenomicVariant.TYPE_INDEL;
		else if(alleles.length==2 && alleles[0].length()==1) return GenomicVariant.TYPE_BIALLELIC_SNV;
		return GenomicVariant.TYPE_STR;
	}

	private static int getGenotypeNumber(CalledGenomicVariant call) {
		if(call.isHeterozygous()) return CalledGenomicVariant.GENOTYPE_HETERO;
		else if (!call.isHomozygousReference()) return CalledGenomicVariant.GENOTYPE_HOMOALT;
		return CalledGenomicVariant.GENOTYPE_HOMOREF;
	}
	private short loadGenotypeQuality(CalledGenomicVariant call) {
		short q = call.getGenotypeQuality();
		if(q>0) return q;
		String [] alleles = call.getAlleles();
		String [] calledAlleles = call.getCalledAlleles();
		VariantCallReport report = call.getCallReport();
		if(report==null || calledAlleles.length==0 || !report.logConditionalsPresent()) return 0;
		double logP;
		Double sum=null;
		if(calledAlleles.length==1) logP = report.getLogConditionalProbability(calledAlleles[0], calledAlleles[0]);
		else logP = report.getLogConditionalProbability(calledAlleles[0], calledAlleles[1]);
		for(int i=0;i<alleles.length;i++) {
			for(int j=i;j<alleles.length;j++) {
				sum = LogMath.logSum(sum, report.getLogConditionalProbability(alleles[i], alleles[j]));
			}
		}
		double logPos = LogMath.logProduct(logP, -sum);
		double pos = LogMath.power10(logPos);
		q = PhredScoreHelper.calculatePhredScore(1-pos);
		call.setGenotypeQuality(q);
		//if(call.getFirst()==376) System.out.println("Conditional P: "+logP+" sum: "+sum+" logpos: "+logPos+" posterior: "+pos+" q: "+q);
		return q;
	}
	private short calculateQuality(List<CalledGenomicVariant> testCalls) {
		short quality = 255;
		for(CalledGenomicVariant call:testCalls) {
			short q = loadGenotypeQuality(call);
			if(q<quality) quality = q;
		}
		return quality;
	}
}
class GoldStandardComparisonCounts {
	public static final int NUM_ROWS_COUNTS = 10;
	private int [][] counts;
	private static final DecimalFormat DF = new DecimalFormat("0.0000");
	private boolean countNonGSAsFP = false;
	private double confidentMbp=3000;
	public GoldStandardComparisonCounts () {
		counts = new int [NUM_ROWS_COUNTS][15];
	}
	
	
	/**
	 * @return the countNonGSAsFP
	 */
	public boolean isCountNonGSAsFP() {
		return countNonGSAsFP;
	}


	/**
	 * @param countNonGSAsFP the countNonGSAsFP to set
	 */
	public void setCountNonGSAsFP(boolean countNonGSAsFP) {
		this.countNonGSAsFP = countNonGSAsFP;
	}
	
	/**
	 * @return the confidentMbp
	 */
	public double getConfidentMbp() {
		return confidentMbp;
	}


	/**
	 * @param confidentMbp the confidentMbp to set
	 */
	public void setConfidentMbp(double confidentMbp) {
		this.confidentMbp = confidentMbp;
	}


	public void update(int firstRow, int lastRow, int column) {
		for(int i=firstRow;i<=lastRow;i++) {
			counts[i][column]++;
		}	
	}
	
	public void print(PrintStream out) {
		for (int i=0;i<counts.length; i++) {
			out.print(""+(i*10));
			for(int j=0;j<counts[0].length;j++) {
				out.print("\t"+counts[i][j]);
			}
			printComparisonStats(counts[i], out);
			out.println();
		}
	}

	private void printComparisonStats(int[] row, PrintStream out) {
		int gsTotal0 = row[0]+row[1]+row[2]+row[9];
		int gsTotal1 = row[3]+row[4]+row[5]+row[10];
		int gsTotal2 = row[6]+row[7]+row[8]+row[11];
		out.print("\t"+gsTotal0+"\t"+gsTotal1+"\t"+gsTotal2);
		
		int testTotal0 = row[0]+row[3]+row[6]+row[12];
		int testTotal1 = row[1]+row[4]+row[7]+row[13];
		int testTotal2 = row[2]+row[5]+row[8]+row[14];
		out.print("\t"+testTotal0+"\t"+testTotal1+"\t"+testTotal2);
		
		double recall1 = 0;
		if(gsTotal1>0) {
			recall1 = (double)row[4] / gsTotal1;
		}
		int fd1 = row[1] + row [7];
		double denomFDR1 = testTotal1;
		if(countNonGSAsFP) fd1 += row[13];
		else denomFDR1 -= row[13];
		double fppm1 = fd1/confidentMbp;
		double fdr1 = 0;
		double precision1 = 1;
		if(denomFDR1>1) {
			fdr1 = (double)fd1 / denomFDR1;
			precision1 = (double)row[4] / denomFDR1;
		}
		double f1 = 0;
		if(precision1 + recall1 > 0) {
			f1 = 2.0*precision1*recall1/(precision1+recall1);
		}
		
		out.print("\t"+DF.format(recall1)+"\t"+fd1+"\t"+DF.format(fppm1)+"\t"+DF.format(fdr1)+"\t"+DF.format(precision1)+"\t"+DF.format(f1));
		
		double recall2 = 0;
		if(gsTotal2>0) {
			recall2 = (double)row[8] / gsTotal2;
		}
		int fd2 = row[2] + row [5];
		double denomFDR2 = testTotal2;
		if(countNonGSAsFP) fd2 += row[14];
		else denomFDR2 -= row[14];
		double fppm2 = fd2/confidentMbp;
		double fdr2 = 0;
		double precision2 = 1;
		if(denomFDR2>0) {
			fdr2 = (double)fd2 / denomFDR2;
			precision2 = (double)row[8] / denomFDR2;
		}
		double f2 = 0;
		if(precision2 + recall2 > 0) {
			f2 = 2.0*precision2*recall2/(precision2+recall2);
		}
		out.print("\t"+DF.format(recall2)+"\t"+fd2+"\t"+DF.format(fppm2)+"\t"+DF.format(fdr2)+"\t"+DF.format(precision2)+"\t"+DF.format(f2));
		
	}
	
}
class GoldStandardHaplotypeReconstruction implements CalledGenomicVariant {
	private List<CalledGenomicVariant> calls;
	private String sequenceName;
	private String reference;
	private String haplotype0;
	private String haplotype1;
	private int numAlleles = 1;
	private byte genotypeN = 0;
	private int first;
	private int last;
	private Map<Integer, Integer> startsH0;
	private Map<Integer, Integer> startsH1;
	
	public GoldStandardHaplotypeReconstruction (String reference, List<CalledGenomicVariant> calls, int first) {
		this.calls = calls;
		this.first = first;
		this.last = first+reference.length()-1;
		this.reference = reference;
		buildHaplotypes();
		boolean b0 = haplotype0.equals(reference);
		boolean b1 = haplotype1.equals(reference);
		boolean b3 = haplotype0.equals(haplotype1);
		
		if(!b0) numAlleles++;
		if(!b1) {
			if(!b3) numAlleles++;
		}
		genotypeN = CalledGenomicVariant.GENOTYPE_HOMOREF;
		if(!b3) genotypeN = CalledGenomicVariant.GENOTYPE_HETERO;
		else if (!b0 || !b1) genotypeN = CalledGenomicVariant.GENOTYPE_HOMOALT; 
		
	}
	
	private void buildHaplotypes() {
		StringBuilder hap0 = new StringBuilder();
		StringBuilder hap1 = new StringBuilder();
		startsH0 = new HashMap<>();
		startsH1 = new HashMap<>();
		int pos = first;
		for(CalledGenomicVariant call:calls) {
			//if(!call.isPhased()) System.out.println("Gold standard call: "+call.getSequenceName()+": "+call.getFirst()+" not phased");
			String [] calledAlleles = call.getCalledAlleles();
			String [] phasedAlleles = call.getPhasedAlleles();
			//System.out.println("Pos: "+pos+" first call: "+call.getFirst()+" first this: "+first);
			if(call.getFirst()>pos) {
				
				String refString = reference.substring(pos-first, call.getFirst()-first); 
				hap0.append(refString);
				hap1.append(refString);
			}
			startsH0.put(call.getFirst(), hap0.length());
			startsH1.put(call.getFirst(), hap1.length());
			if(calledAlleles.length==1) {
				hap0.append(calledAlleles[0]);
				hap1.append(calledAlleles[0]);
			} else if(phasedAlleles!=null) {
				hap0.append(phasedAlleles[0]);
				hap1.append(phasedAlleles[1]);
			} else {
				System.err.println("WARN: Unphased heterozygous gold standard call at "+call.getSequenceName()+": "+call.getFirst());
				hap0.append(calledAlleles[0]);
				hap1.append(calledAlleles[1]);
			}
			pos = call.getLast()+1;
		}
		if(pos<=last) {
			String refEnd = reference.substring(pos-first);
			hap0.append(refEnd);
			hap1.append(refEnd);
		}
		haplotype0 = hap0.toString();
		haplotype1 = hap1.toString();	
	}

	public String [] buildMatchingHaplotypes(List<CalledGenomicVariant> testCalls) {
		StringBuilder hap0 = new StringBuilder();
		StringBuilder hap1 = new StringBuilder();
		int pos = first;
		for(CalledGenomicVariant call:testCalls) {
			String [] calledAlleles = call.getCalledAlleles();
			if(call.getFirst()>pos) {
				String refString = reference.substring(pos-first, call.getFirst()-first); 
				hap0.append(refString);
				hap1.append(refString);
			}
			if(calledAlleles.length==1) {
				//Homozygous alternative
				hap0.append(calledAlleles[0]);
				hap1.append(calledAlleles[0]);
			} else {
				
				String search = hap0+calledAlleles[0];
				String search2 = hap0+calledAlleles[1];
				if(call.getFirst()==97007) {
					System.out.println("Matching variant at "+call.getFirst());
					System.out.println(haplotype0);
					System.out.println(haplotype1);
					System.out.println(search);
					System.out.println(search2);
				}
				
				
				if(haplotype0.startsWith(search) || haplotype1.startsWith(search2)) {
					hap0.append(calledAlleles[0]);
					hap1.append(calledAlleles[1]);
				} else if (haplotype1.startsWith(search) || haplotype0.startsWith(search2)) {
					hap0.append(calledAlleles[1]);
					hap1.append(calledAlleles[0]);
				} else {
					// TODO: In principle it will be an error no matter which phase is used
					hap0.append(calledAlleles[0]);
					hap1.append(calledAlleles[1]);
				}
			}
			pos = call.getLast()+1;
		}
		if(pos<=last) {
			String refEnd = reference.substring(pos-first);
			hap0.append(refEnd);
			hap1.append(refEnd);
		}
		String [] answer = {hap0.toString(),hap1.toString()}; 
		return answer;
	}

	

	/**
	 * @return the haplotype0
	 */
	public String getHaplotype0() {
		return haplotype0;
	}

	/**
	 * @return the haplotype1
	 */
	public String getHaplotype1() {
		return haplotype1;
	}
	
	public int getGenotypeNumber () {
		return genotypeN;
	}

	/**
	 * @return the startsH0
	 */
	public Integer getStartH0(int referencePos) {
		return startsH0.get(referencePos);
	}

	/**
	 * @return the startsH1
	 */
	public Integer getStartsH1(int referencePos) {
		return startsH1.get(referencePos);
	}

	@Override
	public String getSequenceName() {
		return sequenceName;
	}

	@Override
	public int getFirst() {
		return first;
	}

	@Override
	public int getLast() {
		return last;
	}

	@Override
	public int length() {
		return reference.length();
	}

	@Override
	public boolean isPositiveStrand() {
		return true;
	}

	@Override
	public boolean isNegativeStrand() {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public String[] getAlleles() {
		Set <String> allelesSet = new HashSet<>();
		allelesSet.add(reference);
		allelesSet.add(haplotype0);
		allelesSet.add(haplotype1);
		return allelesSet.toArray(new String [0]);
	}

	@Override
	public String getReference() {
		return reference;
	}

	@Override
	public String getId() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void setId(String id) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public short getVariantQS() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public void setVariantQS(short qualityScore) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public boolean isCompatible(GenomicVariant variant) {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public boolean isBiallelic() {
		return numAlleles==2;
	}

	@Override
	public boolean isSNV() {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public byte getType() {
		// TODO Auto-generated method stub
		return GenomicVariant.TYPE_STR;
	}

	@Override
	public void setType(byte type) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public String getSampleId() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void setSampleId(String sampleId) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public String[] getCalledAlleles() {
		Set <String> allelesSet = new HashSet<>();
		allelesSet.add(haplotype0);
		allelesSet.add(haplotype1);
		return allelesSet.toArray(new String [0]);
	}

	@Override
	public byte[] getIndexesCalledAlleles() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public byte[] getAllelesCopyNumber() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public byte getCopyNumber() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public int getTotalReadDepth() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public void setTotalReadDepth(int depth) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public short getGenotypeQuality() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public void setGenotypeQuality(short genotypeQuality) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void makeUndecided() {
		// TODO Auto-generated method stub
		
	}

	@Override
	public boolean isUndecided() {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public boolean isHeterozygous() {
		return genotypeN == 1;
	}

	@Override
	public boolean isHomozygous() {
		return genotypeN != 1;
	}

	@Override
	public boolean isHomozygousReference() {
		return genotypeN == 0;
	}

	@Override
	public void updateAllelesCopyNumberFromCounts(byte totalCopyNumber) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void setAllelesCopyNumber(byte[] allelesCN) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public VariantCallReport getCallReport() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public int[] getAllCounts() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public boolean isPhased() {
		return true;
	}

	@Override
	public String[] getPhasedAlleles() {
		String [] haplotypes = new String[2];
		haplotypes[0] = haplotype0;
		haplotypes[1] = haplotype1;
		return haplotypes;
	}

	@Override
	public byte[] getIndexesPhasedAlleles() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public byte getStrandBiasScore() {
		// TODO Auto-generated method stub
		return 0;
	}
	
}
