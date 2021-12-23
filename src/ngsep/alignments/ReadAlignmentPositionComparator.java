package ngsep.alignments;

import java.util.Comparator;

import ngsep.genome.GenomicRegionPositionComparator;

public class ReadAlignmentPositionComparator implements Comparator<ReadAlignment>{

	private static ReadAlignmentPositionComparator instance = new ReadAlignmentPositionComparator();
	private ReadAlignmentPositionComparator () {
		
	}
	@Override
	public int compare(ReadAlignment aln0, ReadAlignment aln1) {
		if(aln0.isReadUnmapped() && aln1.isReadUnmapped()) return aln0.getReadNumber()-aln1.getReadNumber();
		if(aln0.isReadUnmapped()) return 1;
		if(aln1.isReadUnmapped()) return -1;
		return GenomicRegionPositionComparator.getInstance().compare(aln0, aln1);
	}
	public static ReadAlignmentPositionComparator getInstance() {
		return instance;
	}
}
