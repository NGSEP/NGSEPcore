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

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;
import java.util.logging.Logger;

import ngsep.genome.GenomicRegion;
import ngsep.genome.GenomicRegionSortedCollection;
import ngsep.genome.ReferenceGenome;
import ngsep.genome.io.SimpleGenomicRegionFileHandler;
import ngsep.main.CommandsDescriptor;
import ngsep.main.ProgressNotifier;
import ngsep.variants.CalledGenomicVariant;
import ngsep.variants.DiversityStatistics;
import ngsep.variants.GenomicVariant;
import ngsep.variants.GenomicVariantAnnotation;
import ngsep.variants.SNV;
import ngsep.variants.Sample;


public class VCFFilter {
	
	private Logger log = Logger.getLogger(VCFFilter.class.getName());
    // Genotype filters
    private int minGenotypeQuality = 0;
    private int minCoverage = 0;
    private int minDistance = 0;
    private int minIndividualsGenotyped = 0;
    private boolean filterInvariant = false;
    private boolean filterInvariantReference = false;
    private boolean filterInvariantAlternative = false;
    private boolean keepOnlySNVs = false;
    private double minMAF = 0;
    private double maxMAF = 0.5;
    private double minOH = 0;
    private double maxOH =1;
    private double minGCContent = 40.0;
    private double maxGCContent = 65.0;
    private int maxCNVs = -1;
    private String geneId = null;
    private Set <String> annotations = null;
    
    private Set<String> sampleIds = null;
    private boolean filterSamples = false;
    
    private GenomicRegionSortedCollection<GenomicRegion> regionsToFilter = null;
    private GenomicRegionSortedCollection<GenomicRegion> regionsToSelect = null;
    private ReferenceGenome genome = null;
    
    private ProgressNotifier progressNotifier=null;

    public static void main(String[] args) throws Exception {
		VCFFilter filter = new VCFFilter();
		int i=CommandsDescriptor.getInstance().loadOptions(filter, args);
		boolean systemInput = "-".equals(args[i]);

		if(systemInput) {
			filter.processVariantsFile(System.in, System.out);
		} else {
			String vcfFile = args[i];
			filter.processVariantsFile(vcfFile,System.out);
		}
    }
    
    public ProgressNotifier getProgressNotifier() {
		return progressNotifier;
	}

	public void setProgressNotifier(ProgressNotifier progressNotifier) {
		this.progressNotifier = progressNotifier;
	}
    
    public Logger getLog() {
		return log;
	}

	public void setLog(Logger log) {
		this.log = log;
	}

	public int getMinGenotypeQuality() {
		return minGenotypeQuality;
	}

	public void setMinGenotypeQuality(int minGenotypeQuality) {
		this.minGenotypeQuality = minGenotypeQuality;
	}

	public void setMinGenotypeQuality(Integer minGenotypeQuality) {
		this.setMinGenotypeQuality(minGenotypeQuality.intValue());
	}
	
	public int getMinCoverage() {
		return minCoverage;
	}

	public void setMinCoverage(int minCoverage) {
		this.minCoverage = minCoverage;
	}

	public void setMinCoverage(Integer minCoverage) {
		this.setMinCoverage(minCoverage.intValue());
	}
	
	public int getMinDistance() {
		return minDistance;
	}
	
	public void setMinDistance(int minDistance) {
		this.minDistance = minDistance;
	}
	
	public void setMinDistance(Integer minDistance) {
		this.setMinDistance(minDistance.intValue());
	}

	public int getMinIndividualsGenotyped() {
		return minIndividualsGenotyped;
	}

	public void setMinIndividualsGenotyped(int minIndividualsGenotyped) {
		this.minIndividualsGenotyped = minIndividualsGenotyped;
	}
	
	public void setMinIndividualsGenotyped(Integer minIndividualsGenotyped) {
		this.setMinIndividualsGenotyped(minIndividualsGenotyped.intValue());
	}

	public boolean isKeepOnlySNVs() {
		return keepOnlySNVs;
	}

	public void setKeepOnlySNVs(boolean keepOnlySNVs) {
		this.keepOnlySNVs = keepOnlySNVs;
	}
	
	public void setKeepOnlySNVs(Boolean keepOnlySNVs) {
		this.setKeepOnlySNVs(keepOnlySNVs.booleanValue());
	}

	public boolean isFilterInvariant() {
		return filterInvariant;
	}

	public void setFilterInvariant(boolean filterInvariant) {
		this.filterInvariant = filterInvariant;
	}
	
	public void setFilterInvariant(Boolean filterInvariant) {
		this.setFilterInvariant(filterInvariant.booleanValue());
	}

	public boolean isFilterInvariantReference() {
		return filterInvariantReference;
	}

	public void setFilterInvariantReference(boolean filterInvariantReference) {
		this.filterInvariantReference = filterInvariantReference;
	}
	
	public void setFilterInvariantReference(Boolean filterInvariantReference) {
		this.setFilterInvariantReference(filterInvariantReference.booleanValue());
	}

	public boolean isFilterInvariantAlternative() {
		return filterInvariantAlternative;
	}

	public void setFilterInvariantAlternative(boolean filterInvariantAlternative) {
		this.filterInvariantAlternative = filterInvariantAlternative;
	}
	
	public void setFilterInvariantAlternative(Boolean filterInvariantAlternative) {
		this.filterInvariantAlternative = filterInvariantAlternative;
	}
	
	public boolean isFilterSamples() {
		return filterSamples;
	}

	public void setFilterSamples(boolean filterSamples) {
		this.filterSamples = filterSamples;
	}
	
	public void setFilterSamples(Boolean filterSamples) {
		this.setFilterSamples(filterSamples.booleanValue());
	}

	public double getMinMAF() {
		return minMAF;
	}

	public void setMinMAF(double minMAF) {
		this.minMAF = minMAF;
	}
	
	public void setMinMAF(Double minMAF) {
		this.setMinMAF(minMAF.doubleValue());
	}

	public double getMaxMAF() {
		return maxMAF;
	}

	public void setMaxMAF(double maxMAF) {
		this.maxMAF = maxMAF;
	}
	
	public void setMaxMAF(Double maxMAF) {
		this.setMaxMAF(maxMAF.doubleValue());
	}

	public double getMinOH() {
		return minOH;
	}

	public void setMinOH(double minOH) {
		this.minOH = minOH;
	}
	
	public void setMinOH(Double minOH) {
		this.setMinOH(minOH.doubleValue());
	}

	public double getMaxOH() {
		return maxOH;
	}

	public void setMaxOH(double maxOH) {
		this.maxOH = maxOH;
	}
	
	public void setMaxOH(Double maxOH) {
		this.setMaxOH(maxOH.doubleValue());
	}

	public double getMinGCContent() {
		return minGCContent;
	}

	public void setMinGCContent(double minGCContent) {
		this.minGCContent = minGCContent;
	}
	
	public void setMinGCContent(Double minGCContent) {
		this.setMinGCContent(minGCContent.doubleValue());
	}

	public double getMaxGCContent() {
		return maxGCContent;
	}

	public void setMaxGCContent(double maxGCContent) {
		this.maxGCContent = maxGCContent;
	}
	
	public void setMaxGCContent(Double maxGCContent) {
		this.setMaxGCContent(maxGCContent.doubleValue());
	}
	
	public int getMaxCNVs() {
		return maxCNVs;
	}

	public void setMaxCNVs(int maxCNVs) {
		this.maxCNVs = maxCNVs;
	}
	
	public void setMaxCNVs(Integer maxCNVs) {
		this.setMaxCNVs(maxCNVs.intValue());
	}

	public String getGeneId() {
		return geneId;
	}

	public void setGeneId(String geneId) {
		this.geneId = geneId;
	}
	
	public Set<String> getAnnotations() {
		return annotations;
	}

	public void setAnnotations(Set<String> annotations) {
		this.annotations = annotations;
	}
	
	public void setAnnotations(String csAnns) {
		annotations = new TreeSet<String>();
		annotations.addAll(Arrays.asList(csAnns.split(",")));
	}

	public ReferenceGenome getGenome() {
		return genome;
	}

	public void setGenome(ReferenceGenome genome) {
		this.genome = genome;
	}
	
	public void setGenome(String genomeFile) throws IOException {
		if(genomeFile==null || genomeFile.length()==0) this.genome = null;
		else this.genome = new ReferenceGenome(genomeFile);
	}
	
	public List<GenomicRegion> getRegionsToFilter() {
		return regionsToFilter.asList();
	}

	public void setRegionsToFilter(List<GenomicRegion> regions) {
		this.regionsToFilter = new GenomicRegionSortedCollection<GenomicRegion>(regions);
	}
	
	public void setRegionsToFilter(String regionsFile) throws IOException {
		if(regionsFile==null || regionsFile.length()==0) {
			this.regionsToFilter = null;
			return;
		}
		SimpleGenomicRegionFileHandler regionFileHandler = new SimpleGenomicRegionFileHandler();
		List<GenomicRegion> regions = regionFileHandler.loadRegions(regionsFile);
		this.regionsToFilter = new GenomicRegionSortedCollection<GenomicRegion>(regions);
	}

	public List<GenomicRegion> getRegionsToSelect() {
		return regionsToSelect.asList();
	}

	public void setRegionsToSelect(List<GenomicRegion> regions) {
		this.regionsToSelect = new GenomicRegionSortedCollection<GenomicRegion>(regions);
	}
	
	public void setRegionsToSelect(String regionsFile) throws IOException {
		if(regionsFile==null || regionsFile.length()==0) {
			this.regionsToSelect = null;
			return;
		}
		SimpleGenomicRegionFileHandler regionFileHandler = new SimpleGenomicRegionFileHandler();
		List<GenomicRegion> regions = regionFileHandler.loadRegions(regionsFile);
		this.regionsToSelect = new GenomicRegionSortedCollection<GenomicRegion>(regions);
	}

	

	public Set<String> getSampleIds() {
		return sampleIds;
	}

	public void setSampleIds(Set<String> sampleIds) {
		this.sampleIds = sampleIds;
	}
	
	public void setSampleIds(String sampleIdsFile) throws IOException {
		if(sampleIdsFile==null || sampleIdsFile.length()==0) {
			sampleIds = null;
			return;
		}
		sampleIds = new TreeSet<String>();
		try ( FileReader fr = new FileReader(sampleIdsFile);
				BufferedReader in = new BufferedReader(fr);
		) {
			String line = in.readLine();
			while (line != null) {
				String[] items = line.split("\t| ");
				sampleIds.add(items[0]);
				line = in.readLine();
			}
		} catch (IOException e) {
			sampleIds = null;
			throw e;
		}
	}

	public void processVariantsFile(String vcfFile, PrintStream out) throws IOException {
		VCFFileReader reader = null;
		try {
			reader = new VCFFileReader(vcfFile);
			processVariantsFile(reader, out);
		} finally {
			if(reader!=null) reader.close();
		}
	}
	public void processVariantsFile(InputStream in, PrintStream out) throws IOException {
		VCFFileReader reader = null;
		try {
			reader = new VCFFileReader(in);
			processVariantsFile(reader, out);
		} finally {
			if(reader!=null) reader.close();
		}
	}
	public void processVariantsFile(VCFFileReader reader, PrintStream out) throws IOException {
		VCFFileWriter writer = new VCFFileWriter();
		
		if(log!=null) reader.setLog(log);
		VCFFileHeader header = reader.getHeader();
		VCFFileHeader outHeader = header.cloneEmpty();
		List<Sample> inSamples = header.getSamples();
		Set<String> inSampleIdsWithHeader = header.getSamplesWithHeaderLine().keySet();
		boolean [] selectedSamples = calculateSelectedSamples (inSamples); 
		for(int i=0;i<selectedSamples.length;i++) {
			Sample s = inSamples.get(i);
			if(selectedSamples[i]) outHeader.addSample(s,inSampleIdsWithHeader.contains(s.getId()));
		}
		writer.printHeader(outHeader, out);
		Iterator<VCFRecord> it = reader.iterator();
		int lastPos = -minDistance;
		String lastSeqName = null;
		VCFRecord lastRecord = null;
		int n=0;
		while (it.hasNext()) {
			VCFRecord record = it.next();
			VCFRecord vr = filterSamplesAndGenotypes(record, selectedSamples, outHeader);
			GenomicVariant gv = vr.getVariant();
			if (!gv.getSequenceName().equals(lastSeqName)) {
				if (passFilters(lastRecord)) writer.printVCFRecord(lastRecord, out);
				lastRecord = vr;
				lastSeqName = gv.getSequenceName();
				lastPos = gv.getLast();
				continue;
			}
			// Distance filter
			int first = gv.getFirst();
			if (minDistance <= 0 || first - lastPos > minDistance) {
				if (passFilters(lastRecord)) writer.printVCFRecord(lastRecord, out);
				lastRecord = vr;
			} else {
				lastRecord = null;
			}
			lastSeqName = gv.getSequenceName();
			lastPos = gv.getLast();
			n++;
			if (progressNotifier!=null && n%1000==0) {
				int progress = n/1000;
				if (!progressNotifier.keepRunning(progress)) {
					out.flush();
					return;
				}
			}
		}
		if (passFilters(lastRecord)) writer.printVCFRecord(lastRecord, out);
		out.flush();
    }

	private boolean[] calculateSelectedSamples(List<Sample> inSamples) {
		boolean[]  answer = new boolean[inSamples.size()];
		if(sampleIds==null || sampleIds.size()==0) {
			//If no samples were selected, then keep everybody
			Arrays.fill(answer, true);
			return answer;
		}
		for(int i=0;i<inSamples.size();i++) {
			answer[i] = (sampleIds.contains(inSamples.get(i).getId())!=filterSamples);
		}
		return answer;
	}

	private boolean passFilters(VCFRecord record) {
    	if (record == null) return false;
    	GenomicVariant var = record.getVariant();
    	if(keepOnlySNVs && !(var instanceof SNV)) return false;
    	//System.out.println("Passing filters for record: "+record.getVariant().getSequenceName()+": "+record.getVariant().getFirst());
		if(maxCNVs>=0 && calculateNumCNVs(record)>maxCNVs) return false;
		//System.out.println("Passed CNVs");
    	if(geneId!=null && !isInGene(record)) return false;
    	//System.out.println("Passed Gene");
    	if(annotations!=null && !hasAnnotation(record)) return false;
    	//System.out.println("Passed Annotation");
    	//TODO: Do not recalculate by default
    	DiversityStatistics divStats = DiversityStatistics.calculateDiversityStatistics(record.getCalls(), false);
    	int numCalledAlleles = divStats.getNumCalledAlleles();
    	int [] counts = divStats.getAlleleCounts();
    	double maf = divStats.getMaf();
    	double oh = divStats.getObservedHeterozygosity();
    	//System.out.println("Count 0: "+counts[0]+" Count 1: "+counts[1]+" Alleles: "+numCalledAlleles+". MAF: "+maf);
    	if (filterInvariant && numCalledAlleles < 2) return false;
    	//System.out.println("Passed invariant");
    	//Only reference or undecided.
    	if (filterInvariantReference && numCalledAlleles == 1 && counts[0]>0) return false;
    	//System.out.println("Passed invariant reference");
    	//Only alternative alleles
    	if (filterInvariantAlternative && numCalledAlleles == 1 && counts[0]==0) return false;
    	//System.out.println("Passed invariant alternative");
    	if (divStats.getNumSamplesGenotyped() < minIndividualsGenotyped) return false;
    	//System.out.println("Passed minInd. MAF: "+maf);
    	if (maf < minMAF || maf> maxMAF) return false;
    	//System.out.println("Passed MAF");
    	if (oh < minOH || oh> maxOH) return false;
    	//System.out.println("Passed OH");
    	if (regionsToFilter!=null && regionsToFilter.findSpanningRegions(var).size()>0) return false;
    	//System.out.println("Passed regions to filter");
    	if (regionsToSelect!=null && regionsToSelect.findSpanningRegions(var).size()==0) return false;
    	//System.out.println("Passed regions to select");
    	if (filterGCContent(var)) return false;
    	//System.out.println("Passed GCContent");
    	//Update annotations if passes filters
    	record.addAnnotation(new GenomicVariantAnnotation(var, GenomicVariantAnnotation.ATTRIBUTE_SAMPLES_GENOTYPED, divStats.getNumSamplesGenotyped()));
		record.addAnnotation(new GenomicVariantAnnotation(var, GenomicVariantAnnotation.ATTRIBUTE_NUMBER_ALLELES, divStats.getNumCalledAlleles()));
		record.addAnnotation(new GenomicVariantAnnotation(var, GenomicVariantAnnotation.ATTRIBUTE_ALLELE_FREQUENCY_SPECTRUM, format(divStats.getAlleleCounts())));
		if(divStats.getNumCalledAlleles()==2) record.addAnnotation(new GenomicVariantAnnotation(var, GenomicVariantAnnotation.ATTRIBUTE_MAF, divStats.getMaf()));
    	return true;
    }

	private String format(int[] alleleCounts) {
		StringBuilder answer = new StringBuilder(""+alleleCounts[0]);
		for(int i=1;i<alleleCounts.length;i++) answer.append(","+alleleCounts[i]);
		return answer.toString();
	}

	private boolean hasAnnotation(VCFRecord record) {
		GenomicVariantAnnotation ann = record.getInfoField(GenomicVariantAnnotation.ATTRIBUTE_TRANSCRIPT_ANNOTATION);
		if(ann!=null && annotations.contains(ann.getValue())) return true;
		return false;
	}

	private boolean isInGene(VCFRecord record) {
		GenomicVariantAnnotation ann1 = record.getInfoField(GenomicVariantAnnotation.ATTRIBUTE_TRANSCRIPT_ID);
		if(ann1!=null && geneId.equals(ann1.getValue())) return true;
		GenomicVariantAnnotation ann2 = record.getInfoField(GenomicVariantAnnotation.ATTRIBUTE_GENE_NAME);
		if(ann2!=null && geneId.equals(ann2.getValue())) return true;
		return false;
	}

	public int calculateNumCNVs(VCFRecord record) {
		GenomicVariantAnnotation ann = record.getInfoField(GenomicVariantAnnotation.ATTRIBUTE_IN_CNV);
    	int numCNVs = 0;
		//Support for old boolean format kept for compatibility with old CNV files
		if(ann!=null && ann.getValue()!=null) {
			if((ann.getValue() instanceof Boolean) && ((Boolean)ann.getValue()).booleanValue()) numCNVs = 1;
			else if ((ann.getValue() instanceof Integer)) numCNVs= (Integer)ann.getValue();
			else if ((ann.getValue() instanceof String)) numCNVs+= Integer.parseInt((String)ann.getValue());
		}
		return numCNVs;
	}

    /**
     * 
     * @param v
     * @return
     */
    private boolean filterGCContent(GenomicVariant v) {
    	if (genome != null) {
    		CharSequence seq = genome.getReference(v.getSequenceName(), v.getFirst() - 100, v.getLast() + 100);
    		if (seq == null) return true;
    		String segment = seq.toString();
    		double gcContent = getGCContent(segment);
    		if (gcContent < minGCContent || gcContent > maxGCContent) return true;
	    }
    	return false;
    }

    public double getGCContent(String segment) {
		double countGC = 0;
		double countACGT = 0;
		segment = segment.toUpperCase();
		for (int i = 0; i < segment.length(); i++) {
		    char c = segment.charAt(i);
		    if (c == 'C' || c == 'G')
			countGC++;
		    if (c == 'C' || c == 'G' || c == 'A' || c == 'T')
			countACGT++;
		}
		if (countACGT == 0)
		    return 0;
		return countGC * 100.0 / countACGT;
    }

    public VCFRecord filterSamplesAndGenotypes(VCFRecord record, boolean [] selectedSamples, VCFFileHeader outHeader ) {
    	List<CalledGenomicVariant> genotyped = record.getCalls();
    	List<CalledGenomicVariant> newList = new ArrayList<CalledGenomicVariant>();
    	boolean recalCNVs = false;
    	for (int i=0;i<genotyped.size();i++) {
    		CalledGenomicVariant cv = genotyped.get(i);
    		if(!selectedSamples[i]) {
    			recalCNVs = true;
    			continue; 
    		}
    		short q = cv.getGenotypeQuality();
    		int depth = cv.getTotalReadDepth();
    		if (q < minGenotypeQuality || depth < minCoverage) {
    			cv.makeUndecided();
    		}
    		newList.add(cv);
    	}
    	VCFRecord newRecord = new VCFRecord(record.getVariant(),record.getFilters(), record.getInfoFields(), record.getFieldsFormat(), newList, outHeader);
    	if(recalCNVs) updateSamplesWithCNVs(newRecord, outHeader.getSamples());
    	return newRecord;
    }

	private void updateSamplesWithCNVs(VCFRecord record, List<Sample> outSamples) {
		GenomicVariant var = record.getVariant();
		List<CalledGenomicVariant> calls = record.getCalls();
		int numCNVs = 0;
		for(int i=0;i<calls.size();i++) {
			CalledGenomicVariant cv = calls.get(i);
			byte normalPloidy = outSamples.get(i).getNormalPloidy(); 
			byte copyNumber = cv.getCopyNumber();
			if(copyNumber!=normalPloidy) {
				numCNVs++;
			}
		}
		if(numCNVs==0) {
			record.removeAnnotation(GenomicVariantAnnotation.ATTRIBUTE_IN_CNV);
		} else {
			record.addAnnotation(new GenomicVariantAnnotation(var, GenomicVariantAnnotation.ATTRIBUTE_IN_CNV, numCNVs));
		}
	}
}
