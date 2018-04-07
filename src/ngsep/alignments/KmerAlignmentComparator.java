package ngsep.alignments;
import java.util.Comparator;
/**
 * Compatarator to sort KmerAlignment by position in reference
 * @author German Andrade
 *
 */
public class KmerAlignmentComparator implements Comparator<KmerAlignment>
{

	public static KmerAlignmentComparator getInstance() 
	{
		return new KmerAlignmentComparator();
	}

	@Override
	public int compare(KmerAlignment o1, KmerAlignment o2) 
	{
		int cmp=0;
		if(o1.getReadAlignment().getFirst()!=o2.getReadAlignment().getFirst())
		{
			cmp = o1.getReadAlignment().getFirst()>o2.getReadAlignment().getFirst()?1:-1;
		}
		return cmp;
	}

}
