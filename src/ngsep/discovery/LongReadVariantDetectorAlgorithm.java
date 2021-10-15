package ngsep.discovery;

import java.util.List;
import java.util.Map;

import ngsep.genome.GenomicRegion;
import ngsep.genome.GenomicRegionSortedCollection;
import ngsep.genome.ReferenceGenome;
import ngsep.variants.GenomicVariant;

public interface LongReadVariantDetectorAlgorithm {
			
	public void setRefGenome(ReferenceGenome genome);
	
	public void setSignatures(GenomicRegionSortedCollection<GenomicRegion> signatures);
	
	public GenomicRegionSortedCollection<GenomicVariant> callVariants();
}
