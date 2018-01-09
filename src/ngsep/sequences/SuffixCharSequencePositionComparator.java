package ngsep.sequences;
import java.util.Comparator;

public class SuffixCharSequencePositionComparator implements Comparator<Integer> 
{
	private CharSequence word;
	private int length;
	public SuffixCharSequencePositionComparator(CharSequence word) 
	{
		super();
		this.word = word;
		this.length = word.length();
	}
	
	@Override
	public int compare(Integer o1, Integer o2) {
		while (o1<length && o2<length) {
			char c1 = word.charAt(o1);
			char c2 = word.charAt(o2); 
			if(c1!=c2) return c1-c2;
			o1++;
			o2++;
		}
		if( o2<length) return -1;
		if( o1<length) return 1;
		return 0;
	}

}
