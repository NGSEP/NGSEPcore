package ngsep.sequences;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import ngsep.genome.GenomicRegion;

public class FMIndexSingleSequence implements Serializable {
	private String sequenceName;
	private Map<Integer,Integer> partialSuffixArray = new HashMap<>();
	private List<Integer []> tallyIndexes = new ArrayList<>();
	private int tallyDistance;
	private String bwt;
	private int [] characterCounts;
	private String alphabet;
	public FMIndexSingleSequence(String seqName, CharSequence sequence) {
		
	}
	public FMIndexSingleSequence(QualifiedSequence sequence) {
		this (sequence.getName(),sequence.getCharacters());
	}
	public List<GenomicRegion> search (String searchSequence) {
		List<GenomicRegion> alignments = new ArrayList<>();
		return alignments;
	}
}
