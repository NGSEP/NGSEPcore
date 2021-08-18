package ngsep.discovery;

import java.util.List;
import java.util.Map;

import ngsep.genome.ReferenceGenome;
import ngsep.variants.GenomicVariant;

public interface LongReadVariantCallerAlgorithm {
	
	public void setReferenceGenome(ReferenceGenome ref);
	
	public void setPartitionSize(int partitions);
	
	public void setSignatures(Map<String, List<Signature>> signatures);
	
	public List<GenomicVariant> callVariants();
}
