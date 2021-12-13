package ngsep.discovery;

import java.util.List;
import java.util.Map;

import ngsep.genome.GenomicRegionSortedCollection;
import ngsep.variants.GenomicVariant;
import ngsep.variants.GenomicVariantImpl;

public interface LongReadVariantDetectorClusteringAlgorithm {
	
	public Map<String, List<List<Integer>>> callVariantClusters();
	
}
