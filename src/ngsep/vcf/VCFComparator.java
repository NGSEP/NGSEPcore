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
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.logging.Logger;

import ngsep.genome.GenomicRegionComparator;
import ngsep.genome.ReferenceGenome;
import ngsep.main.CommandsDescriptor;
import ngsep.main.ProgressNotifier;
import ngsep.sequences.QualifiedSequenceList;
import ngsep.variants.CalledGenomicVariant;
import ngsep.variants.GenomicVariant;
import ngsep.variants.SNV;

/**
 * COmpares two VCF files looking for differences between the genotypes
 * @author Jorge Duitama
 *
 */
public class VCFComparator {
	private Logger log = Logger.getLogger(VCFComparator.class.getName());
	private ProgressNotifier progressNotifier=null;
	// Sequence names
	private ReferenceGenome genome;
	// Genotyped filter
	private double minPCTGenotyped = 50;
	// Differences filter
	private double maxPCTDiffs = 5;
	
	private List<String> samples1;
	private List<String> samples2;
	// Single file statistics
	private int [] genotypedF1;
	private int [] genotypedF2;
	// Comparison statistics
	private int [] [] genotypedBothFiles;
	private int [] [] homozygousDiffs;
	private int [] [] heterozygousDiffs;
	

	public static void main(String[] args) throws Exception {
		VCFComparator instance = new VCFComparator();
		
		int i = CommandsDescriptor.getInstance().loadOptions(instance, args);
		String referenceGenome = args[i++];
		String vcf1 = args[i++];
		String vcf2 = args[i++];
		instance.genome = new ReferenceGenome(referenceGenome);
		instance.calculateDifferences(vcf1, vcf2);
		instance.printReport(System.out);
	}

	
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


	public double getMinPCTGenotyped() {
		return minPCTGenotyped;
	}

	public void setMinPCTGenotyped(double minPCTGenotyped) {
		this.minPCTGenotyped = minPCTGenotyped;
	}
	
	public void setMinPCTGenotyped(Double minPCTGenotyped) {
		this.setMinPCTGenotyped(minPCTGenotyped.doubleValue());
	}

	public double getMaxPCTDiffs() {
		return maxPCTDiffs;
	}

	public void setMaxPCTDiffs(double maxPCTDiffs) {
		this.maxPCTDiffs = maxPCTDiffs;
	}
	
	public void setMaxPCTDiffs(Double maxPCTDiffs) {
		this.setMaxPCTDiffs(maxPCTDiffs.doubleValue());
	}

	


	public ReferenceGenome getGenome() {
		return genome;
	}


	public void setGenome(ReferenceGenome genome) {
		this.genome = genome;
	}


	public void calculateDifferences(String vcf1, String vcf2) throws IOException {
		// VCF files
		
		boolean debug = false;
		try (VCFFileReader in1 = new VCFFileReader(vcf1);
			 VCFFileReader in2 = new VCFFileReader(vcf2);){ 
			if(log!=null) {
				in1.setLog(log);
				in2.setLog(log);
			}
			in1.setLoadMode(VCFFileReader.LOAD_MODE_MINIMAL);
			in2.setLoadMode(VCFFileReader.LOAD_MODE_MINIMAL);
			QualifiedSequenceList seqNames = genome.getSequencesMetadata();
			in1.setSequences(seqNames);
			in2.setSequences(seqNames);
		
			AlleleCompatibilityGenomicVariantComparator comparator = new AlleleCompatibilityGenomicVariantComparator(seqNames);
			// Samples ids (name sample)
			samples1 = in1.getHeader().getSampleIds();
			samples2 = in2.getHeader().getSampleIds();
			startCounts();
			// Iterator VCF files
			Iterator<VCFRecord> it1 = in1.iterator();
			Iterator<VCFRecord> it2 = in2.iterator();
			VCFRecord r1 = null;
			VCFRecord r2 = null;
			if(it1.hasNext()) r1 = it1.next();
			if(it2.hasNext()) r2 = it2.next();
			int n=0;
			while (r1!=null && r2!=null) {
				GenomicVariant g1 = r1.getVariant();
				GenomicVariant g2 = r2.getVariant();
				int cmp = comparator.compare(g1, g2);
				if(cmp < 0) {
					processSingleRecord(r1,1);
					n++;
					if(it1.hasNext()) r1 = it1.next();
					else {
						processSingleRecord(r2,2);
						break;
					}
				} else if (cmp>0) {
					processSingleRecord(r2,2);
					if(debug) throw new RuntimeException("Variant at "+r2.getSequenceName()+":"+r2.getFirst()+"-"+r2.getLast()+" not found in first file. Next var f1 at "+r1.getSequenceName()+":"+r1.getFirst()+"-"+r1.getLast());
					if(it2.hasNext()) r2 = it2.next();
					else {
						processSingleRecord(r1,1);
						break;
					}
				} else {
					n++;
					processSingleRecord(r1,1);
					processSingleRecord(r2,2);
					compareRecords(r1,r2);
					if(it1.hasNext()) r1 = it1.next();
					else break;
					if(it2.hasNext()) r2 = it2.next();
					else break;
				}
				if (progressNotifier!=null && n%1000==0 && !progressNotifier.keepRunning(n/1000) ) return;
			}
			while(it1.hasNext()) {
				r1 = it1.next();
				processSingleRecord(r1,1);
				n++;
				if (progressNotifier!=null && n%1000==0 && !progressNotifier.keepRunning(n/1000) ) return;
			}
			while(it2.hasNext()) {
				r2 = it2.next();
				processSingleRecord(r2,2);
				n++;
				if (progressNotifier!=null && n%1000==0 && !progressNotifier.keepRunning(n/1000) ) return;
			}
		}
		
	}

	private void startCounts() {
		genotypedF1 = new int [samples1.size()];
		Arrays.fill(genotypedF1, 0);
		genotypedF2 = new int [samples2.size()];
		Arrays.fill(genotypedF2, 0);
		genotypedBothFiles = new int [samples1.size()][samples2.size()];
		heterozygousDiffs = new int [samples1.size()][samples2.size()];
		homozygousDiffs = new int [samples1.size()][samples2.size()];
		for(int i=0;i<samples1.size();i++) {
			for(int j=0;j<samples2.size();j++) {
				genotypedBothFiles[i][j] = 0;
				heterozygousDiffs[i][j] = 0;
				homozygousDiffs[i][j] = 0;
			}
		}
		
	}

	private void processSingleRecord(VCFRecord r, int fileId) {
		List<CalledGenomicVariant> calls = r.getCalls();
		for(int i=0;i<calls.size();i++) {
			CalledGenomicVariant call = calls.get(i);
			if(!call.isUndecided()) {
				if(fileId == 1) genotypedF1[i]++;
				else genotypedF2[i]++;
			}
		}
		
	}

	private void compareRecords(VCFRecord r1, VCFRecord r2) {
		List<CalledGenomicVariant> calls1 = r1.getCalls();
		List<CalledGenomicVariant> calls2 = r2.getCalls();
		for(int i=0;i<calls1.size();i++) {
			CalledGenomicVariant call1 = calls1.get(i);
			if(call1.isUndecided()) continue;
			for(int j=0;j<calls2.size();j++) {
				CalledGenomicVariant call2 = calls2.get(j);
				if(call2.isUndecided()) continue;
				genotypedBothFiles[i][j]++;
				if(call1.isHeterozygous() || call2.isHeterozygous()) {
					if(call1.isHeterozygous()!=call2.isHeterozygous()) {
						heterozygousDiffs[i][j]++;
					}
				} else if (call1.isHomozygousReference()!=call2.isHomozygousReference()) {
					homozygousDiffs[i][j]++;
				}
			}
		}	
	}
	
	public void printReport(PrintStream out) {
		
		out.println("Sample1" + "\t" + "Sample2" + "\t" + "SNPsSample1" + "\t" + "SNPsSample2" + "\t"+ "SNPsBoth" + "\t" + "No.HeteroDifferences" + "\t" + "%HeteroDifferences"  +"\t" + "No.HomoDifferences" + "\t" + "%HomoDifferences"+ "\t" + "No.TotalDifferences" + "\t" + "%TotalDifferences");
		for(int i=0;i<samples1.size();i++) {
			Map<Integer, List<String>> linesSample = new TreeMap<Integer, List<String>>();
			int gen1 = genotypedF1[i];
			for(int j=0;j<samples2.size();j++) {
				int gen2 = genotypedF2[j];
				int genBoth = genotypedBothFiles[i][j];
				if(0.01*minPCTGenotyped*gen1>genBoth) continue;
				if(0.01*minPCTGenotyped*gen2>genBoth) continue;
				int homoDiff = homozygousDiffs[i][j];
				int heteroDiff = heterozygousDiffs[i][j];
				double heteroDiffPercent = 100.0*((double)heteroDiff)/genBoth;
				double homoDiffPercent = 100.0*((double)homoDiff)/genBoth;
				int totalDiff = homoDiff+heteroDiff;
				double totalDiffPercent = 100.0*((double)(homoDiff + heteroDiff))/genBoth;
				if(0.01*maxPCTDiffs*genBoth<totalDiff) continue;
				String line = samples1.get(i)+"\t"+samples2.get(j)+"\t"+gen1+"\t"+gen2+"\t"+genBoth+"\t"+heteroDiff+"\t"+ heteroDiffPercent+"\t"+homoDiff+"\t"+homoDiffPercent+"\t"+totalDiff+"\t"+ totalDiffPercent;

				List<String> linesDiffs = linesSample.get(totalDiff);
				if(linesDiffs==null) {
					linesDiffs = new ArrayList<String>();
					linesSample.put(totalDiff, linesDiffs);
				}
				linesDiffs.add(line);
			}
			for(Integer d:linesSample.keySet()) {
				for(String line:linesSample.get(d)) out.println(line);
			}
		}
	}
	private class AlleleCompatibilityGenomicVariantComparator implements Comparator<GenomicVariant> {

		private QualifiedSequenceList seqNames;
		private GenomicRegionComparator internalComparator;
		
		public AlleleCompatibilityGenomicVariantComparator(QualifiedSequenceList sequenceNames) {
			seqNames = sequenceNames;
			internalComparator = new GenomicRegionComparator(sequenceNames);
		}
		@Override
		public int compare(GenomicVariant v1, GenomicVariant v2) {
			int cmp1 = internalComparator.compare(v1, v2);
			if(cmp1>2 || cmp1<-2) return cmp1;
			if(v1 instanceof SNV && v2 instanceof SNV) return cmp1;
			//if(v2.getFirst()==753845) System.out.println("Variant 1: "+v1.getSequenceName()+":"+v1.getFirst()+"-"+v1.getLast()+"Comparator result: "+cmp1);
			if(v2.getFirst() - v1.getLast() > 2 || v1.getFirst() - v2.getLast() > 2) return cmp1;
			int firstRegion = Math.min(v1.getFirst(), v2.getFirst());
			int lastRegion = Math.max(v1.getLast(), v2.getLast());
			//if(v2.getFirst()==753845) System.out.println("Getting alleles between: "+firstRegion+" and "+lastRegion);
			List<String> alleleStrs1 = new ArrayList<>(buildAlleleStrings(v1,firstRegion,lastRegion));
			List<String> alleleStrs2 = new ArrayList<>(buildAlleleStrings(v2,firstRegion,lastRegion));
			//if(v2.getFirst()==753845) System.out.println("Alleles 1: "+alleleStrs1.toString()+" Alleles 2: "+alleleStrs2.toString());
			if(alleleStrs1.size()!=alleleStrs2.size()) return cmp1!=0?cmp1: alleleStrs1.size()-alleleStrs2.size();
			for(int i=0;i<alleleStrs1.size();i++) {
				String a1 = alleleStrs1.get(i);
				String a2 = alleleStrs2.get(i);
				if(!a1.equals(a2))  return cmp1!=0?cmp1:a1.compareTo(a2);
			}
			return 0;
		}
		private Set<String> buildAlleleStrings(GenomicVariant v, int firstRegion, int lastRegion) {
			String seqName = v.getSequenceName();
			String left = genome.getReference(seqName, Math.max(1, firstRegion-3), Math.max(1, v.getFirst()-1)).toString();
			int length = seqNames.get(seqName).getLength();
			String right = genome.getReference(seqName, Math.min(length, v.getLast()+1), Math.min(length, lastRegion+3)).toString();
			Set<String> allelesSet = new TreeSet<>();
			String [] alleles = v.getAlleles();
			for(int i=0;i<alleles.length;i++) {
				allelesSet.add((left+alleles[i]+right).toUpperCase());
			}
			return allelesSet;
		}
		
	}
}

