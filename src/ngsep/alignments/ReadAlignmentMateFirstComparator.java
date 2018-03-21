package ngsep.alignments;

import java.util.Comparator;

public class ReadAlignmentMateFirstComparator implements Comparator<ReadAlignment>{

	public static ReadAlignmentMateFirstComparator getInstance() {
		// TODO Auto-generated method stub
		return new ReadAlignmentMateFirstComparator();
	}

	public int compare(ReadAlignment o1, ReadAlignment o2) 
	{
		int cmp=0;
		if(o1.getMateFirst()!=o2.getMateFirst())
		{
			cmp = o1.getMateFirst()>o2.getMateFirst()?1:-1;
		}
		return cmp;
	}

}
