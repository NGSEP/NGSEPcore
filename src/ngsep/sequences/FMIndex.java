package ngsep.sequences;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

import ngsep.genome.GenomicRegion;

public class FMIndex implements Serializable {
	private List<FMIndexSingleSequence> singleSequenceIndexes = new ArrayList<>();

	/**
	 * 
	 * @param filename
	 */
	public static FMIndex load(String filename) {
		return null;
	}
	
	public void save (String filename) {
		
	}
	
	public List<GenomicRegion> search (String searchSequence) {
		List<GenomicRegion> alignments = new ArrayList<>();
		return alignments;
	}
	public List<GenomicRegion> search (String sequenceName, String searchSequence) {
		List<GenomicRegion> alignments = new ArrayList<>();
		return alignments;
	}
	
}
