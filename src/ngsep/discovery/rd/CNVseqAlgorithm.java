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
package ngsep.discovery.rd;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.logging.Logger;

import JSci.maths.statistics.NormalDistribution;
import ngsep.genome.ReferenceGenome;
import ngsep.main.CommandsDescriptor;
import ngsep.main.ProgressNotifier;
import ngsep.main.io.ParseUtils;
import ngsep.math.PhredScoreHelper;
import ngsep.variants.CalledCNV;
import ngsep.variants.GenomicVariant;
import ngsep.variants.GenomicVariantImpl;


/**
 * Implementation of the CNV-seq algorithm, as published by Xie and Tammi, 2009:
 * "CNV-seq, a new method to detect copy number variation using high-throughput sequencing". BMC bioinformatics, 10:80.
 * @author Juan Fernando de la Hoz
 */
public class CNVseqAlgorithm {
	//------------------------------------------------------------------------------------------------------------------------------------------------
	//								ATTRIBUTES
	//------------------------------------------------------------------------------------------------------------------------------------------------
	// file management
	// input
	private String reference;
	private String bamXfile;
	ReadDepthDistribution sampleX;
	private String bamYfile;
	ReadDepthDistribution sampleY;
	// output
	private static Logger log = Logger.getLogger(CNVseqAlgorithm.class.getName());
	private String outFile;
	public static final String EMPTY="*";
	public static final String SEP="\t";
    private ProgressNotifier progressNotifier=null;
    private int progress = 0;

	// parameters needed for CNV ratio calculation
	private long genomeSize;
	private double readNumX;
	private double readNumY;
	private double lambdaX;
	private double lambdaY;
	
	// lists to manage the read count for each bin 
	private List<ReadDepthBin> seqBinsX;
	private List<ReadDepthBin> seqBinsY;
	private List<Double> rdListX;
	private List<Double> rdListY;
	private List<Double> ratioRDList;
	private List<Double> ratioCNVList;
		
	// optional parameters
	private boolean gcCorrection = false;
	private boolean wholeGenomePrnt = false;
	private boolean bonferroni = false;
	private double pValue = 0.001;
	private int binSize = 100;
	
	//------------------------------------------------------------------------------------------------------------------------------------------------
	//								MAIN METHODS
	//------------------------------------------------------------------------------------------------------------------------------------------------
	/**
	 * Receives the parameters from the command line interface and distributes the duties
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		if (args.length == 0 || args[0].equals("-h") || args[0].equals("--help")){
			CommandsDescriptor.getInstance().printHelp(CNVseqAlgorithm.class);
			return;
		}
		CNVseqAlgorithm cnvSeq = new CNVseqAlgorithm();

		// parse arguments from input
		int i = 0;
		CNVseqAlgorithm.log.info("Reading from input.");
		while(i < args.length && args[i].charAt(0) == '-'){
			if(args[i].equals("-binSize")){
				i ++;
				cnvSeq.binSize = Integer.parseInt(args[i]);
			} else if(args[i].equals("-p")){
				i ++;
				cnvSeq.pValue = Double.parseDouble(args[i]);
			} else if(args[i].equals("-w")){
				cnvSeq.wholeGenomePrnt = true;
			} else if(args[i].equals("-g")){
				cnvSeq.gcCorrection = true;
			} else if(args[i].equals("-b")){
				cnvSeq.bonferroni = true;
			} else {
				System.err.println("Unrecognized option "+args[i]);
				CommandsDescriptor.getInstance().printHelp(CNVseqAlgorithm.class);
				return;
			}
			i ++;
		}
		// obtain BAMs and reference file paths
		cnvSeq.bamXfile = args[i++];
		cnvSeq.bamYfile = args[i++];
		cnvSeq.reference = args[i++];
		cnvSeq.loadFiles(log);
		
		// create the output file with the extension .cnvSeq
		cnvSeq.outFile = args[i++]+".cnvSeq";
		cnvSeq.runCNVseq(log);
		
		CNVseqAlgorithm.log.info("Done.");
	}
	
	/**
	 * Constructor
	 */
	public CNVseqAlgorithm(){
		log.info("Starting CNV-seq");
	}

	//------------------------------------------------------------------------------------------------------------------------------------------------
	//								EXECUTION METHODS
	//------------------------------------------------------------------------------------------------------------------------------------------------
	/**
	 * Creates the objects necessary for the implementation:
	 * a genome assembly using the reference genome,
	 * both ReadDepthDistributions using the bam for each sample
	 * @throws IOException
	 */
	public void loadFiles(Logger log) throws IOException{
		// create the object GenomeAssembly using the reference genome file and set the genomeSize attribute
		log.info("Loading reference genome.");
		advanceNotifier();
		ReferenceGenome refGenome = new ReferenceGenome(reference);
		log.info("Loaded reference genome. Total number of sequences: "+refGenome.getNumSequences());
		genomeSize = refGenome.getTotalLength();
		
		// create both instances of ReadDepthDistribution, one for each BAM file
		log.info("Loading bam file for sample X. This can take a couple of minutes, please wait...");
		advanceNotifier();
		sampleX = new ReadDepthDistribution(refGenome, binSize);
		sampleX.processAlignments(bamXfile);
		log.info("Loading bam file for sample Y. This can take a couple of minutes, please wait...");
		advanceNotifier();
		sampleY = new ReadDepthDistribution(refGenome, binSize);
		sampleY.processAlignments(bamYfile);
		readNumX = sampleX.getTotalReads();
		readNumY = sampleY.getTotalReads();
		log.info("Both bam files loaded."+"\n"+ "bamX with "+ readNumX +" reads, bamY with "+ readNumY +" reads.");
		advanceNotifier();
	}
	
	/**
	 * This is the principal method, which executes the whole algorithm
	 * and call the other methods. 
	 * @param out
	 * @throws FileNotFoundException 
	 */
	public void runCNVseq(Logger log) throws FileNotFoundException{		
		// perform GC correction
		if(gcCorrection){
			log.info("Correcting by GC content.");
			advanceNotifier();
			sampleX.correctDepthByGCContent();
			sampleY.correctDepthByGCContent();
		}
		
		// get lists of bins
		seqBinsX = new ArrayList<ReadDepthBin>();
		seqBinsX = sampleX.getAllBins();
		seqBinsY = new ArrayList<ReadDepthBin>();
		seqBinsY = sampleY.getAllBins();
		advanceNotifier();

		// get the list of the CNV ratio between sampleX:sampleY
		log.info("Obtaining read depth for each bin.");
		advanceNotifier();
		rdListX = getRDList(seqBinsX);
		rdListY = getRDList(seqBinsY);
		double totalCountRatio = (readNumY/readNumX);
		log.info("Calculating read depth ratio for each bin.");
		advanceNotifier();
		ratioRDList = getRDratios(rdListX, rdListY);
		log.info("Estimating CNV ratio for each bin.");
		advanceNotifier();
		ratioCNVList = calculateCNVratios(ratioRDList, totalCountRatio);	

		// perform statistical significance test
		log.info("Calculating p-value for each CNV ratio.");
		advanceNotifier();
		lambdaX = readNumX*binSize/genomeSize; 
		lambdaY = readNumY*binSize/genomeSize; 
		List<Double> pValueList = calculatePvalue(ratioCNVList, ratioRDList, lambdaX, lambdaY);
		
		// print the output
		if(bonferroni) pValue = pValue / pValueList.size();
		if(wholeGenomePrnt) pValue = 0.5;
		log.info("Writing CNV list.");
		advanceNotifier();
		PrintStream out = new PrintStream(outFile);
		printCNVList(seqBinsX, rdListX, rdListY, ratioCNVList, pValueList, out);
		out.flush();
		out.close();
		advanceNotifier();
	}

	//------------------------------------------------------------------------------------------------------------------------------------------------
	//								AUXILIARY METHODS
	//------------------------------------------------------------------------------------------------------------------------------------------------

	/**
	 * This method creates a list with read counts for each bin in a list of bins
	 * takes into account GC correction
	 * @param List containing all bins in the genome
	 * @param boolean whether to perform CG correction or not
	 * @return List of the read depth per bin 
	 */
	public List<Double> getRDList(List<ReadDepthBin> bins){
		List<Double> binsRD = new ArrayList<Double>();
		
		for(ReadDepthBin bin:bins){
			if(gcCorrection){
				binsRD.add(bin.getCorrectedReadDepth());
			} else{
				binsRD.add(bin.getRawReadDepth());
			}
		}
		return binsRD;
	}

	/**
	 * This method calculates the ratio between read counts for each bin given in the lists
	 * @param List of the read depth for each bin in sample X
	 * @param List of the read depth for each bin in sample Y
	 * @return List containing the read depth ratio for each bin X:Y
	 */
	public List<Double> getRDratios(List<Double> countSampleX, List<Double> countSampleY){
		List<Double> readCountRatios = new ArrayList<Double>();
		for(int i=0;i<countSampleX.size();i++){
			double countRatio = countSampleX.get(i)/countSampleY.get(i);
			readCountRatios.add(countRatio);
		}
		return readCountRatios;
	}
	
	/**
	 * This method calculates the CNV ratio for each bin,
	 * taking into account the number of reads for each sample
	 * The CNV ratio is given without any logarithmic transformation, mean should be 1.
	 * @param List containing the read depth ratio for each bin X:Y
	 * @param double constant to normalize by the total amount of reads for each sample
	 * @return List containing predicted CNV ratios for each bin X:Y
	 */
	public List<Double> calculateCNVratios(List<Double> binCountRatio, double totalCountRatio){
		List<Double> cnvRatios = new ArrayList<Double>();
		for(int i=0; i < binCountRatio.size();i++){
			double predictedCNVratio = binCountRatio.get(i) * totalCountRatio;
			cnvRatios.add(predictedCNVratio);
		}
		return cnvRatios;
	}

	/**
	 * This method calculates the p-value for each read counts ratio,
	 * taking into account the distribution of both read counts and the Geary-Hinkley transformation
	 * @param List containing predicted CNV ratios for each bin X:Y
	 * @param List containing the read depth ratio for each bin X:Y
	 * @param double Average number of reads per window in random sequencing of sample X
	 * @param double Average number of reads per window in random sequencing of sample Y
	 * @return List of p-values are given for each CNV ratio in the list 
	 */
	public List<Double> calculatePvalue(List<Double> cnvRatios, List<Double> rdRatios, double lambX, double lambY){
		List<Double> pValues = new ArrayList<Double>();
		NormalDistribution normDist = new NormalDistribution();
		// to get p-value, the cumulative normal distribution is used
		for(int i = 0; i < rdRatios.size(); i++){
			double t = z2tTransform(rdRatios.get(i),lambX,lambY);
			if(cnvRatios.get(i) >= 1){
				pValues.add(1-normDist.cumulative(t));
			}
			else {
				pValues.add(normDist.cumulative(t));
			}
		}
		return pValues;
	}
	
	/**
	 * this method performs the Geary-Hinkley trasnformation of z to t
	 * @param z read depth ratio for one bin 
	 * @param lambdaX mean and variance in sample X read depth distribution
	 * @param lambdaY mean and variance in sample Y read depth distribution
	 * @return z transformed to t
	 */
	private double z2tTransform(double z, double lambdaX, double lambdaY){
		return ((lambdaY*z)-lambdaX)/(StrictMath.sqrt((lambdaY*(z*z))+lambdaX));
	}
	
	/**
	 * This method calculates the p-value for each read counts ratio,
	 * @param List posList, ordered bins to extract their position in genome
	 * @param List readDepthX, read depth for every bin of the sample X
	 * @param List readDepthX, read depth for every bin of the sample Y
	 * @param List CNVratios, estimated CNV ratio for each bin
	 * @param List pvalueList, calculated p-value for each CNV ratio
	 * @param String outFile
	 * @throws FileNotFoundException 
	 * @post a table is created listing all the information from CNV-seq algorithm. each line is a bin in the genome.
	 */
	public void printCNVList(List<ReadDepthBin> posList, List<Double> readDepthX, List<Double> readDepthY, List<Double> cnvRatioList, List<Double> pValueList, PrintStream out) {
		log.info("The maximum p-value reported is: "+ pValue);
		DecimalFormat df = ParseUtils.ENGLISHFMT;
		for ( int i = 0 ; i < posList.size() ; i++) {
			if(pValueList.get(i) <= pValue){
				out.print(posList.get(i).getSequenceName()+SEP);
				out.print(posList.get(i).getFirst()+SEP);
				out.print(posList.get(i).getLast()+SEP);
				if(readDepthX.get(i) == null) out.print(EMPTY+SEP);
				else out.print(df.format(readDepthX.get(i))+SEP);
				if(readDepthY.get(i) == null) out.print(EMPTY+SEP);
				else out.print(df.format(readDepthY.get(i))+SEP);
				if(cnvRatioList.get(i) == null) out.print(EMPTY+SEP);
				else out.print(cnvRatioList.get(i)+SEP);
				if(pValueList.get(i) == null) out.print(EMPTY);
				else out.print(pValueList.get(i));
				out.println();
			}
		}
	}
	
	/**
	 * TODO
	 * AT THE MOMENT NO ONE CALLS THIS METHOD, BUT IT WILL SOON BE IMPLEMENTED TO MODIFY THE OUTPUT FORMAT
	 * merges large CNVs if several, continuous bins, have low p-values
	 * @param List posList, ordered bins to extract their position in genome
	 * @param List CNVratios, estimated CNV ratio for each bin
	 * @param List pvalueList, calculated p-value for each CNV ratio
	 * @return a list of all the detected CNVs with their average p-value and their average CNV
	 */
	public List<CalledCNV> mergeCNV(List<ReadDepthBin> posList, List<Double> cnvRatioList, List<Double> pValueList){
		List<CalledCNV> mergedCNVs = new ArrayList<CalledCNV>();
		for(int i = 0 ; i < posList.size() ; i++){
			if(pValueList.get(i) < pValue){
				int startBin = i;
				double avrgPval = 0;
				double avrgCNVratio = 0;
				while(pValueList.get(i) < pValue && posList.get(i).getSequenceName().equals(posList.get(startBin).getSequenceName())){
					avrgPval += pValueList.get(i);
					avrgCNVratio += cnvRatioList.get(i);
					i++;
				}
				int endBin = i-1;
				GenomicVariantImpl cnv = new GenomicVariantImpl(posList.get(startBin).getSequenceName(),posList.get(startBin).getFirst(),posList.get(endBin).getLast(),GenomicVariant.TYPE_CNV);
				CalledCNV largeCNV = new CalledCNV(cnv,(float)avrgCNVratio);
				avrgPval /= (endBin-startBin);
				avrgCNVratio /= (endBin-startBin);
				largeCNV.setGenotypeQuality(PhredScoreHelper.calculatePhredScore(avrgPval));
				mergedCNVs.add(largeCNV);
			}
		}
		return mergedCNVs;
	}
	
	public void advanceNotifier(){
		progress ++;
		if (progressNotifier!=null && !progressNotifier.keepRunning(progress)) {
			return;
		}
	}
	
	
	//------------------------------------------------------------------------------------------------------------------------------------------------
	//								SETTERS AND GETTERS
	//------------------------------------------------------------------------------------------------------------------------------------------------

	public String getReference() {
		return reference;
	}	
	
	public void setReference(String reference) {
		this.reference = reference;
	}	
	
	public String getBamXfile() {
		return bamXfile;
	}	
	
	public void setBamXfile(String bamXfile) {
		this.bamXfile = bamXfile;
	}	
	
	public String getBamYfile() {
		return bamYfile;
	}	
	
	public void setBamYfile(String bamYfile) {
		this.bamYfile = bamYfile;
	}	
	
	public String getOutFile() {
		return outFile;
	}	
	
	public void setOutFile(String outFile) {
		this.outFile = outFile;
	}	
	
	public int getBinSize() {
		return binSize;
	}	
	
	public void setBinSize(int binSize) {
		this.binSize = binSize;
	}	
	
	public boolean isGcCorrection() {
		return gcCorrection;
	}	
	
	public void setGcCorrection(boolean gcCorrection) {
		this.gcCorrection = gcCorrection;
	}	
	
	public boolean isBonferroni() {
		return bonferroni;
	}	
	
	public void setBonferroni(boolean bonferroni) {
		this.bonferroni = bonferroni;
	}	
	
	public boolean isWholeGenomePrnt() {
		return wholeGenomePrnt;
	}	
	
	public void setWholeGenomePrnt(boolean wholeGenomePrnt) {
		this.wholeGenomePrnt = wholeGenomePrnt;
	}	
	
	public double getpValue() {
		return pValue;
	}	
	
	public void setpValue(double pValue) {
		this.pValue = pValue;
	}	
	
	public ProgressNotifier getProgressNotifier() {
		return progressNotifier;
	}	
	
	public void setProgressNotifier(ProgressNotifier progressNotifier) {
		this.progressNotifier = progressNotifier;
	}
	
	/*
	 * this methods are thought to be implemented when the optimal window size calculation is done

	public int getBinSize() {
		return binSize;
	}
	public void setBinSize(int binSize) {
		this.binSize = binSize;
	}
	*/
}
