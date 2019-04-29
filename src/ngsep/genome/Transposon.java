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

public class Transposon {
	
	private GenomicRegionImpl startRegion;
	
	private GenomicRegionImpl endRegion;
	
	private String type;
	
	public Transposon(GenomicRegionImpl start, GenomicRegionImpl end, String ptype) {
		startRegion = start;
		endRegion = end;
		type = ptype;
	}

	public GenomicRegionImpl getStartRegion() {
		return startRegion;
	}

	public void setStartRegion(GenomicRegionImpl startRegion) {
		this.startRegion = startRegion;
	}

	public GenomicRegionImpl getEndRegion() {
		return endRegion;
	}

	public void setEndRegion(GenomicRegionImpl endRegion) {
		this.endRegion = endRegion;
	}

	public String getType() {
		return type;
	}

	public void setType(String type) {
		this.type = type;
	}
	
}
