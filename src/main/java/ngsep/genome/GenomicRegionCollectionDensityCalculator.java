package ngsep.genome;

import java.io.IOException;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;

import ngsep.genome.io.SimpleGenomicRegionFileHandler;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;

public class GenomicRegionCollectionDensityCalculator {

	public static final int DEF_WINDOW_LENGTH = 100000;
	private ReferenceGenome genome;
	private int windowLength = DEF_WINDOW_LENGTH;
	private boolean printTotalBpPerWindow = false;
	public static void main(String[] args) throws Exception {
		GenomicRegionCollectionDensityCalculator instance = new GenomicRegionCollectionDensityCalculator();
		instance.genome = new ReferenceGenome(args[0]);
		instance.calculateDensity (args[1], System.out);
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
	 * @return the windowLength
	 */
	public int getWindowLength() {
		return windowLength;
	}


	/**
	 * @param windowLength the windowLength to set
	 */
	public void setWindowLength(int windowLength) {
		this.windowLength = windowLength;
	}


	/**
	 * @return the printTotalBpPerWindow
	 */
	public boolean isPrintTotalBpPerWindow() {
		return printTotalBpPerWindow;
	}


	/**
	 * @param printTotalBpPerWindow the printTotalBpPerWindow to set
	 */
	public void setPrintTotalBpPerWindow(boolean printTotalBpPerWindow) {
		this.printTotalBpPerWindow = printTotalBpPerWindow;
	}


	public void calculateDensity(String filename, PrintStream out) throws IOException {
		SimpleGenomicRegionFileHandler handler = new SimpleGenomicRegionFileHandler();
		List<GenomicRegion> regions = handler.loadRegions(filename);
		calculateDensity (regions, out);
		
	}
	public void calculateDensity(Collection<? extends GenomicRegion> regions, PrintStream out) {
		QualifiedSequenceList seqNames = genome.getSequencesMetadata();
		GenomicRegionSortedCollection<GenomicRegion> sortedRegions = new GenomicRegionSortedCollection<>(seqNames);
		sortedRegions.addAll(regions);
		for(QualifiedSequence seq:seqNames) {
			int seqLen = seq.getLength();
			//Counts in genomic coordinates
			int [] densityMap = new int [seqLen+1];
			int nw = seqLen/windowLength;
			if(seqLen%windowLength>0)nw++;
			int [] windowCounts = new int [nw];
			Arrays.fill(densityMap, 0);
			Arrays.fill(windowCounts, 0);
			List<GenomicRegion> seqRegions = sortedRegions.getSequenceRegions(seq.getName()).asList();
			for(GenomicRegion region:seqRegions) {
				for(int i=region.getFirst();i<=region.getLast();i++) {
					densityMap[i]++;
				}
				int midPoint = (region.getFirst()+region.getLast())/2;
				int w = (midPoint-1)/windowLength;
				windowCounts[w]++;
			}
			//Print counts per window
			int totalW = 0;
			int countW = 0;
			int w = 0;
			for(int i=1;i<=seqLen;i++) {
				totalW+=densityMap[i];
				if(densityMap[i]>0) countW++;
				if(i%windowLength==0 || i==seqLen) {
					out.print(seq.getName()+"\t"+(w*windowLength+1)+"\t"+i+"\t"+windowCounts[w]+"\t"+countW);
					if(printTotalBpPerWindow)out.print("\t"+totalW);
					out.println();
					totalW = countW = 0;
					w++;
				}
			}
		}
		
	}

}
