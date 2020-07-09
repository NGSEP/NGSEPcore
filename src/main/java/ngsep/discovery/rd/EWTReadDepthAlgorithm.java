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

import java.io.PrintStream;	
import java.lang.Math;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.logging.Logger;

import JSci.maths.statistics.NormalDistribution;
import ngsep.genome.GenomicRegionComparator;
import ngsep.genome.GenomicRegionImpl;
import ngsep.genome.ReferenceGenome;
import ngsep.math.PhredScoreHelper;
import ngsep.sequences.QualifiedSequenceList;
import ngsep.variants.CalledCNV;
import ngsep.variants.GenomicVariant;
import ngsep.variants.GenomicVariantImpl;
import ngsep.variants.io.GFFVariantsFileHandler;

/**
 * @author Juan Fernando De la Hoz
 * This class implements the Event-Wise-Testing algorithm proposed by Yoon et al. 2009:
 * "Sensitive and accurate detection of copy number variants using read depth of coverage". Genome Research, 19:1586-1592.
 * and integrates it with NGSEP
 */
public class EWTReadDepthAlgorithm implements SingleSampleReadDepthAlgorithm {

	public static final String SOURCE_EWT = "EWT";
	//----------------------------------------------------------------
	//								ATTRIBUTES
	//----------------------------------------------------------------

	private Logger log = Logger.getLogger(EWTReadDepthAlgorithm.class.getName());
	private ReadDepthDistribution readDepthDistribution;
	private	GenomicRegionComparator comparator;
	private QualifiedSequenceList sequences;
	
	private Map <String, List<RDbinProbabilities>> probabilities;
	private double readDepthMean;
	private double readDepthSDeviation;
	private double falsePositiveRate = 0.05;
	private int binSize = ReadDepthDistribution.DEFAULT_BIN_SIZE;
	private byte normalPloidy = 2;
	private long genomeSize = 0;
	private boolean merge = true;
	private boolean filter = true;
	
	//------------------------------------------------------------------
	//								MAIN METHODS
	//------------------------------------------------------------------
	
	/**
	 * For calling CNVs using the EWT algorithm through the Command Line,
	 * receives as input the reference FASTA file, the alignments in BAM format,
	 * and the prefix for the output. Outputs a GFF file with CNVs information. 
	 * @param args "--fpr <INT falsePositiveRate> -b <INT binSize> -p <INT ploidy> -noMerge -noFilter [Reference] [Alignments] [OutputPrefix]"
	 * @throws Exception
	 */
	public static void main(String[] args) throws Exception {
		EWTReadDepthAlgorithm ewt = new EWTReadDepthAlgorithm();
		ewt.log.info("EWT algorithm");
		ewt.log.info("loading parameters");
		int i=0;
		while(i<args.length && args[i].charAt(0)=='-') {
			if("--fpr".equals(args[i])) {
				i++;
				ewt.setFalsePositiveRate(Double.parseDouble(args[i]));
				ewt.log.info("False Positive Rate set to " + ewt.getFalsePositiveRate());
			} else if ("-b".equals(args[i])) {
				i++;
				ewt.setBinSize(Integer.parseInt(args[i]));
				ewt.log.info("Bin Size set to " + ewt.getBinSize());
			} else if ("-p".equals(args[i])) {
				i++;
				ewt.setNormalPloidy(Byte.parseByte(args[i]));
				ewt.log.info("Normal Ploidy set to " + ewt.getNormalPloidy());
			} else if ("-noMerge".equals(args[i])) {
				ewt.setMerge(false);
				ewt.log.info("Events will not be merged");
			} else if ("-noFilter".equals(args[i])) {
				ewt.setFilter(false);
				ewt.log.info("Events will not be filtered for low difference with mean");
			}
			i++;
		}
		String reference = args[i++];
		String alignment = args[i++];
		String outFile = args[i++];
		
		// load input
		ewt.log.info( "Loading reference genome" );
		ReferenceGenome genome = new ReferenceGenome( reference );
		ewt.log.info( "Loaded genome reference. Sequences: " + genome.getNumSequences() );
		ewt.log.info( "Dividing genome into " + ewt.binSize + "bp windows" );
		ewt.readDepthDistribution = new ReadDepthDistribution( genome, ewt.binSize );
		ewt.log.info( "Loading alignment file" );
		ewt.readDepthDistribution.processAlignments(alignment);
		ewt.log.info( "Processed alignment file: " + alignment );
		ewt.readDepthDistribution.correctDepthByGCContent();
		ewt.log.info( "GC-content corrected" );
		ewt.readDepthDistribution.calculateReadDepthDistParameters();
		
		// call CNVs
		List<CalledCNV> cnvs = ewt.callCNVs();
		
		// print CNVs
		GFFVariantsFileHandler cnvFH = new GFFVariantsFileHandler();
		PrintStream out = new PrintStream( outFile + ".gff" );
		cnvFH.saveVariants( cnvs, out );
		out.flush();
		out.close();
		ewt.log.info( "Saved cnvs" );
	}

	//------------------------------------------------------------------------
	//								EXECUTION METHODS
	//------------------------------------------------------------------------
	
	@Override
	public List<CalledCNV> callCNVs() {
		log.info( "Calling CNVs using EWT algorithm" );
		// set variables
		sequences = readDepthDistribution.getSequences();
		comparator = new GenomicRegionComparator(sequences);
		readDepthMean = readDepthDistribution.getMeanReadDepth();
		readDepthSDeviation = readDepthDistribution.getSigmaReadDepth();
		List<Interval> detectedCNVs = new ArrayList<Interval>();
		
		// obtain z-score, upperLimit probability and lowerLimit probability for each bin, from the read depth
		log.info( "Transforming read depth to z-score" );
		rdBinsToRDProb();
		
		// iterate over each chromosome
		for ( String seqName : sequences.getNamesStringList() ) {
			log.info( "Calling CNVs for sequence " + seqName );
			
			// get all bins in chromosome
			List<RDbinProbabilities> seqProbs = probabilities.get(seqName);
			int numProbs = seqProbs.size();
			
			// calculate N, and make intervals of all suggested lengths (2 <= l <= N)
			double significance;																											
			for ( int l = 2 ; (significance = Math.pow( (falsePositiveRate / (numProbs / l)) , (1.0 / l) )) < 0.50  ; l++ ) {				
				//log.info( "Searching CNVs in " + l + "00bp-sized intervals" );
				List<Interval> intervals = getIntervals(seqName, seqProbs, l);
				for ( int i = 0 ; i < intervals.size() ; i++ ) {
					
					// test for duplications and deletions
					Interval event = intervals.get(i);
					if ( event.getMaxUpperProb() < significance ) detectedCNVs.add(event);
					else if ( event.getMaxLowerProb() < significance ) detectedCNVs.add(event);
				}
			}
		}
		
		// sort, filter, merge and output
		log.info( detectedCNVs.size() + " total events detected by EWT algorithm" );
		Collections.sort(detectedCNVs, comparator);
		List<Interval> filteredEvents = filter ? filterEvents(detectedCNVs) : detectedCNVs;
		List<Interval> reportCNVs = merge ? mergeEvents(filteredEvents) : filteredEvents;
		log.info( "Called " + reportCNVs.size() + " CNVs" );
		
		return intervalToCalledCNV(reportCNVs);
	}

	//------------------------------------------------------------------------------
	//								AUXILIARY METHODS
	//------------------------------------------------------------------------------
	
	/**
	 * Fills in the probabilities TreeMap, saving one list of bins with their associated
	 * probabilities for each chromosome (or sequence in the FASTA file).
	 * Takes advantage of the bins Map present in the readDepthDistribution object.
	 */
	public void rdBinsToRDProb () {
		probabilities = new TreeMap<String, List<RDbinProbabilities>>();

		// iterate over all the bins in the genome
		for ( String seqName : sequences.getNamesStringList() ) {
			List<ReadDepthBin> seqBins = readDepthDistribution.getBins(seqName);
			List<RDbinProbabilities> seqProbs = new ArrayList<RDbinProbabilities>();
			
			log.info( "normalizing read depth for bins in " + seqName );
			for ( int i = 0 ; i < seqBins.size() ; i++ ) {
				ReadDepthBin bin = seqBins.get(i);

				// change the nature of each bin to its probability
				RDbinProbabilities binP = new RDbinProbabilities(bin.getSequenceName(), bin.getFirst(), bin.getLast(), bin.getGcContent(), bin.getCorrectedReadDepth());

				// calculate Z-score and add to the new list
				binP.setzScore( (bin.getCorrectedReadDepth() - readDepthMean) / readDepthSDeviation );
				seqProbs.add(binP);
			}
			
			probabilities.put(seqName, seqProbs);
		}
	}
	
	/**
	 * Generate fixed-sized intervals from a list of bins in a single chromosome.
	 * @param seqName String, the name of the chromosome or scaffold
	 * @param windows List, the list of bins in the sequence
	 * @param l	INT, the number of bins in the interval
	 * @return List a list of intervals of the given length
	 */
	public List<Interval> getIntervals( String seqName, List<RDbinProbabilities> windows, int l ) {
		List<Interval> intervals = new ArrayList<Interval>();
		
		for ( int i = 0 ; i < windows.size()-l ; i += l ) {
			List<RDbinProbabilities> bins = new ArrayList<RDbinProbabilities>();
			for ( int j = 0 ; j < l ; j++ ) {
				bins.add( windows.get(i+j) );
			}
			Interval interval = new Interval(seqName, bins.get(0).getFirst(), bins.get(l-1).getLast(), bins);
			intervals.add(interval);
		}
		
		return intervals;
	}
	
	/**
	 * Filter events with low absolute difference from the average Read Depth, and the ones in repetitive regions,
	 * this means, the events containing reads with multiple alignments
	 * @param intervals List
	 * @return filteredEvents List  
	 */
	public List<Interval> filterEvents( List<Interval> intervals ) {
		List<Interval> filteredEvents = new ArrayList<Interval>();
		
		for ( int i = 0 ; i < intervals.size() ; i++ ) {
			if ( intervals.get(i).getMedianRD() > (1.25 * readDepthMean) || intervals.get(i).getMedianRD() < (0.75 * readDepthMean) ) {
				filteredEvents.add( intervals.get(i) );
			}
		}
		
		log.info( (intervals.size() - filteredEvents.size()) + " events were filtered out because of low absolute difference from mean read depth" );
		return filteredEvents;
	}
	
	/**
	 * Merge events with overlapping or adjacent bins in a single, large event.
	 * @param intervals List
	 * @return mergedIntervals List
	 */
	public List<Interval> mergeEvents ( List<Interval> intervals ) {
		List<Interval> mergedIntervals = new ArrayList<Interval>();
		
		for ( int i = 0 ; i < intervals.size() ; i++ ) {
			Interval next = intervals.get(i);
			Interval last = (mergedIntervals.size()>0) ? mergedIntervals.get(mergedIntervals.size()-1) : null;
			
			// check if it's the first CNV in the list and look for overlapping regions (-1), assuming sorted list.
			int cmp = (last == null) ? -3 : comparator.compare(last, next);	
			// mark for merging if CNVs are adjacent
			if ( cmp == -2 && (last.getLast() + 1 == next.getFirst()) )	cmp = -1;
			// keep mark for merging only if both CNVs go in the same direction
			if ( cmp == -1 && !((next.getMedianRD()>readDepthMean && last.getMedianRD()>readDepthMean) || (next.getMedianRD()<readDepthMean && last.getMedianRD()<readDepthMean)) ) cmp = -2;
			
			// Merge if marked
			if ( cmp == -1 ) {
				List<RDbinProbabilities> mergedBins = last.getBins();
				for ( RDbinProbabilities bin : next.getBins() ) {
					if ( !mergedBins.contains(bin) )
						mergedBins.add(bin);
				}
				last.setLast( mergedBins.get(mergedBins.size()-1).getLast() );
				last.calculateProbs();
			} else {
				// keep if not marked
				mergedIntervals.add(next);
			}
		}
		
		log.info( (intervals.size() - mergedIntervals.size()) + " events were merged because of overlapping or adjacent regions" );
		return mergedIntervals;
	}

	/**
	 * Transforms a list of intervals, to a list of detected CNVs,
	 * takes information from the bins list of the interval and 
	 * adds it to the new CNV, including: length, type, source, 
	 * pValue, number of copies, and number of fragments aligned.
	 * @param detectedCNVs List
	 * @return cnvsList List
	 */
	public List<CalledCNV> intervalToCalledCNV( List<Interval> detectedCNVs ) {
		List<CalledCNV> cnvsList = new ArrayList<CalledCNV>();
		
		for ( int j = 0 ; j < detectedCNVs.size() ; j++ ) {
			Interval interval = detectedCNVs.get(j);
			GenomicVariantImpl cnvObj = new GenomicVariantImpl( interval.getSequenceName(), interval.getFirst(), interval.getLast(), GenomicVariant.TYPE_CNV );
			float copies = (float) (normalPloidy * interval.getAverageRD() / readDepthMean);
			CalledCNV cnv = new CalledCNV( cnvObj, copies );
			cnv.setSource( SOURCE_EWT );
			double pvalue = ( interval.getMedianRD() > readDepthMean ) ? interval.getMaxUpperProb() : interval.getMaxLowerProb() ;
			cnv.setGenotypeQuality(PhredScoreHelper.calculatePhredScore(pvalue));
			cnv.setTotalReadDepth( (int)Math.round(interval.getSum()) );
			cnvsList.add( cnv );
		}
		
		return cnvsList;
	}
	
	//--------------------------------------------------------------------------------------
	//								SETTERS AND GETTERS
	//--------------------------------------------------------------------------------------

	public ReadDepthDistribution getReadDepthDistribution() {
		return readDepthDistribution;
	}
	public void setReadDepthDistribution(ReadDepthDistribution readDepthDistribution) {
		this.readDepthDistribution = readDepthDistribution;
	}
	public Logger getLog() {
		return log;
	}
	public void setLog(Logger log) {
		if (log == null) throw new NullPointerException("Log can not be null");
		this.log = log;
	}
	public double getFalsePositiveRate() {
		return falsePositiveRate;
	}
	public void setFalsePositiveRate(double falsePositiveRate) {
		this.falsePositiveRate = falsePositiveRate;
	}
	public int getBinSize() {
		return binSize;
	}
	public void setBinSize(int binSize) {
		this.binSize = binSize;
	}
	public byte getNormalPloidy() {
		return normalPloidy;
	}
	public void setNormalPloidy(byte normalPloidy) {
		this.normalPloidy = normalPloidy;
	}
	public long getGenomeSize() {
		return genomeSize;
	}
	public void setGenomeSize(long genomeSize) {
		this.genomeSize = genomeSize;
	}
	public boolean isMerge() {
		return merge;
	}
	public void setMerge(boolean merge) {
		this.merge = merge;
	}
	public boolean isFilter() {
		return filter;
	}
	public void setFilter(boolean filter) {
		this.filter = filter;
	}
	
}
	
	//------------------------------------------------------------------------------------------
	//								AUXILIARY CLASSES
	//------------------------------------------------------------------------------------------

/**
 * Genomic region containing a list of bins with normalized z-scores. 
 * When declared, it calculates several statistics about the probability of being a CNV. 
 * If modified, the calculateProbs() method should be called.
 */
class Interval extends GenomicRegionImpl {
	private double maxUpperProb;
	private double maxLowerProb;
	private List<RDbinProbabilities> bins;
	private double medianRD;
	private double averageRD;
	private double sum;
	
	/**
	 * Generates the interval and calculates the related probabilities
	 * @param sequenceName
	 * @param first
	 * @param last
	 * @param bins List of bins with normalized probabilities
	 */
	public Interval(String sequenceName, int first, int last, List<RDbinProbabilities> bins) {
		super(sequenceName, first, last);
		this.bins = bins;
		calculateProbs();
	}
	/**
	 * Calculate the total number of reads mapped to this interval, 
	 * and the following statistics from the contained bins:
	 * mean read depth, median read depth, maximum upper-tail probability,
	 * and maximum lower-tail probability.
	 */
	public void calculateProbs() {
		maxUpperProb = 0.0;
		maxLowerProb = 0.0;
		sum = 0;

		// get the maxUpperTail and maxLowerTail probabilities
		for ( int i = 0 ; i < bins.size() ; i++ ) {
			if ( maxUpperProb < bins.get(i).getUpperTailP() ) maxUpperProb = bins.get(i).getUpperTailP();
			if ( maxLowerProb < bins.get(i).getLowerTailP() ) maxLowerProb = bins.get(i).getLowerTailP();
		}
		
		// calculate mean and median Read Depths
		ArrayList<Double> readDepths = new ArrayList<Double>();
		for ( int i = 0 ; i < bins.size() ; i++ ) {
			readDepths.add( bins.get(i).getCorrectedReadDepth() );
			sum += bins.get(i).getCorrectedReadDepth();
		}
		averageRD = sum/bins.size();
		Collections.sort(readDepths);
	    if ( readDepths.size() % 2 == 1 ) {
	    	medianRD = readDepths.get( ( (readDepths.size() + 1) / 2) -1 );
	    } else {
			double lower = readDepths.get( (readDepths.size() / 2) -1 );
			double upper = readDepths.get( readDepths.size() / 2 );
			medianRD = (lower + upper) / 2;
	    }
	}
	public double getMedianRD() {
		return medianRD;
	}
	public double getAverageRD() {
		return averageRD;
	}
	public double getSum() {
		return sum;
	}
	public List<RDbinProbabilities> getBins() {
		return bins;
	}
	public void setBins(List<RDbinProbabilities> bins) {
		this.bins = bins;
	}
	public double getMaxUpperProb() {
		return maxUpperProb;
	}
	public double getMaxLowerProb() {
		return maxLowerProb;
	}
}

/**
 * Extention of the ReadDepthBin Class, includes more information from each bin:
 * a z-score, the upper-tail probability and the lower-tail probability. 
 */
class RDbinProbabilities extends ReadDepthBin {
	NormalDistribution normDist = new NormalDistribution();
	private double zScore;
	private double upperTailP;
	private double lowerTailP;

	/**
	 * Bins should be initialized with the information of gc content and the corrected read depth 
	 * @param sequenceName
	 * @param first
	 * @param last
	 * @param gcContent
	 * @param readDepth
	 */
	public RDbinProbabilities(String sequenceName, int first, int last, double gcContent, double readDepth) {
		super(sequenceName, first, last, gcContent);
		this.setCorrectedReadDepth(readDepth);
	}
	/**
	 * When the z-score is set, the upper-tail and lower-tail probabilities are calculated 
	 * @param zScore
	 */
	public void setzScore(double zScore) {
		this.zScore = zScore;
		// calculate upperLimit and lowerLimit using the cumulative normal distribution
		lowerTailP = normDist.cumulative(zScore);
		upperTailP = 1.0 - normDist.cumulative(zScore);
	}
	public double getzScore() {
		return zScore;
	}
	public double getUpperTailP() {
		return upperTailP;
	}
	public double getLowerTailP() {
		return lowerTailP;
	}
}
