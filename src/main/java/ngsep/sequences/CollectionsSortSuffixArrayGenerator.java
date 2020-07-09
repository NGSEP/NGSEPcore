package ngsep.sequences;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

public class CollectionsSortSuffixArrayGenerator implements SuffixArrayGenerator, Comparator<Integer> {

	private CharSequence sequence;
	private List<Integer> saList;
	public CollectionsSortSuffixArrayGenerator(CharSequence sequence) {
		this.sequence = sequence;
		saList = new ArrayList<Integer>(sequence.length()+1);
		for(int i=0;i<=sequence.length();i++) saList.add(i);
		Collections.sort(saList,this);
	}
	@Override
	public int [] getSuffixArray() {
		int [] sa = new int [saList.size()];
		for(int i=0;i<sa.length;i++) sa[i] = saList.get(i);
		return sa;
	}
	@Override
	public int compare(Integer i1, Integer i2) {
		if(i1==i2) return 0;
		int n = sequence.length();
		while (i1<n && i2<n) {
			char c1 = sequence.charAt(i1);
			char c2 = sequence.charAt(i2);
			if(c1!=c2) return c1-c2;
			i1++;
			i2++;
		}
		if(i1==n) return -1;
		return 1;
	}

}
