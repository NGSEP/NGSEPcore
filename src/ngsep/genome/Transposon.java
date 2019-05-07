package ngsep.genome;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;

import ngsep.alignments.ReadAlignment;
import ngsep.genome.ReferenceGenome;
import ngsep.math.Distribution;
import ngsep.sequences.DefaultKmersMapImpl;
import ngsep.sequences.FMIndex;
import ngsep.sequences.KmersCounter;
import ngsep.sequences.QualifiedSequence;

public class Transposon extends GenomicRegionImpl{
	
	private String type;
	
	private int score;
	
	public Transposon(String sequenceName, int first, int last, String ptype, int pscore) {
		super(sequenceName, first, last);
		type = ptype;
		score = pscore;
	}

	public String getType() {
		return type;
	}

	public void setType(String type) {
		this.type = type;
	}
	
	public int getScore() {
		return score;
	}
	
	public void setScore(int score) {
		this.score = score;
	}
	
}
