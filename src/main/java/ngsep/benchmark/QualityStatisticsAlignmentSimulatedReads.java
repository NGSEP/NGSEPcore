package ngsep.benchmark;

import java.io.IOException;
import java.io.PrintStream;
import java.util.Iterator;

import ngsep.alignments.ReadAlignment;
import ngsep.alignments.io.ReadAlignmentFileReader;

public class QualityStatisticsAlignmentSimulatedReads {

	private String alignmentsFile;
	private int minAlignmentQuality;
	
	//Statistics
	private int alignedReads = 0;
	private int unalignedReads = 0;
	private int alignedReadsPassMQ = 0;
	private int alignedReadsBelowMQ = 0;
	private double squaredError = 0;
	private int alignmentsProperPairs = 0;
	
	public String getAlignmentsFile() {
		return alignmentsFile;
	}
	public void setAlignmentsFile(String alignmentsFile) {
		this.alignmentsFile = alignmentsFile;
	}

	public int getMinAlignmentQuality() {
		return minAlignmentQuality;
	}
	public void setMinAlignmentQuality(int minAlignmentQuality) {
		this.minAlignmentQuality = minAlignmentQuality;
	}
	
	public static void main(String[] args) throws Exception {
		QualityStatisticsAlignmentSimulatedReads instance = new QualityStatisticsAlignmentSimulatedReads();
		instance.alignmentsFile = args[0];
		instance.minAlignmentQuality = Integer.parseInt(args[1]);
		instance.run();

	}
	public void run() throws IOException {
		processFile(alignmentsFile);
		printStatistics(System.out);
	}
	public void processFile(String alignmentsFile) throws IOException {
		try (ReadAlignmentFileReader reader = new ReadAlignmentFileReader(alignmentsFile)){
			reader.setLoadMode(ReadAlignmentFileReader.LOAD_MODE_FULL);
			Iterator<ReadAlignment> iter = reader.iterator();
			while (iter.hasNext())
			{
				ReadAlignment aln = iter.next();
				if(aln.isSecondary()) continue;
				String name = aln.getReadName();
				String [] items = name.split("_");
				int readLength = aln.getReadLength();
				short mq = aln.getAlignmentQuality();
 
				int firstIndex = getFirstIndex (items);
				int expectedStart;
				if (!aln.isPaired() || aln.isPositiveStrand()){
					expectedStart = Integer.parseInt(items[firstIndex]);
				} else {
					expectedStart = Integer.parseInt(items[firstIndex+1]) - readLength + 1;
				}
				if (aln.isReadUnmapped()) {
					unalignedReads++;
				} else {
					alignedReads++;
					if(mq>=minAlignmentQuality) {
						int alignedStart = aln.getFirst(); 
						alignedReadsPassMQ++;
						int difference = expectedStart - alignedStart;
						
						//Distance larger than 2x read length 
						if ((Math.abs(difference)) <= 2*readLength) { 
							squaredError +=  Math.pow(difference, 2);
						} else {
							squaredError += 4*Math.pow(difference, 2);
		                }
						if(aln.isPaired() && aln.isProperPair()) {
					    	alignmentsProperPairs ++;
					    }
					} else {
						alignedReadsBelowMQ++;
					}
				}
			}
		}
	}
	private int getFirstIndex(String[] items) {
		for(int i=1;i<items.length;i++) {
			try {
				Integer.parseInt(items[i]);
				return i;
			} catch (NumberFormatException e) {
				
			}
		}
		return -1;
	}
	public void printStatistics(PrintStream out) {
		double rmse = 0;
		if(alignedReadsPassMQ>0) rmse= squaredError/alignedReadsPassMQ;
		rmse = Math.sqrt(rmse);
		int totalNegatives = alignedReadsBelowMQ+unalignedReads;
		out.println("STATS\t"+alignedReadsPassMQ+"\t"+totalNegatives+"\t"+rmse+"\t"+alignedReads+"\t"+unalignedReads);
		if(alignmentsProperPairs>0) out.println("Reads aligned as proper pairs: "+alignmentsProperPairs);
	}
}
