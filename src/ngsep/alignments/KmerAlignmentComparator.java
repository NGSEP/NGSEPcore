package ngsep.alignments;

import java.util.Comparator;

public class KmerAlignmentComparator implements Comparator<KmerAlignment>{

	public static KmerAlignmentComparator getInstance() {
		// TODO Auto-generated method stub
		return new KmerAlignmentComparator();
	}

	@Override
	public int compare(KmerAlignment o1, KmerAlignment o2) {
		// TODO Auto-generated method stub
		int cmp=0;
		if(o1.getReadAlignment().getFirst()!=o2.getReadAlignment().getFirst())
		{
			cmp = o1.getReadAlignment().getFirst()>o2.getReadAlignment().getFirst()?1:-1;
		}
		return cmp;
	}

}
