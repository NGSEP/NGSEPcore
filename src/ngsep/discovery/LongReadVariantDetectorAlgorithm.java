package ngsep.discovery;

import java.util.List;
import java.util.Map;

import ngsep.genome.GenomicRegionSortedCollection;
import ngsep.genome.ReferenceGenome;
import ngsep.variants.GenomicVariant;

public interface LongReadVariantDetectorAlgorithm {
		
	//public void setSignatures(Map<String, List<GenomicVariant>> signatures);
	
	public void setSignatures(GenomicRegionSortedCollection<GenomicVariant> signatures);
	
	public GenomicRegionSortedCollection<GenomicVariant> callVariants();
}
