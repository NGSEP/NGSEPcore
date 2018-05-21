package ngsep.alignments;
import java.util.Comparator;

import ngsep.genome.GenomicRegionPositionComparator;
/**
 * Compatarator to sort KmerAlignment by position in reference
 * @author German Andrade
 *
 */
public class KmerAlignmentComparator implements Comparator<KmerAlignment>
{
	private GenomicRegionPositionComparator internalCMP = GenomicRegionPositionComparator.getInstance();
	private static final KmerAlignmentComparator instance = new KmerAlignmentComparator();
	private KmerAlignmentComparator () {
		
	}
	public static KmerAlignmentComparator getInstance() 
	{
		return instance;
	}

	@Override
	public int compare(KmerAlignment o1, KmerAlignment o2) 
	{
		if(o1.isNegativeStrand()!=o2.isNegativeStrand()) return o1.isPositiveStrand()?-1:1;
		return internalCMP.compare(o1, o2);
	}

}
