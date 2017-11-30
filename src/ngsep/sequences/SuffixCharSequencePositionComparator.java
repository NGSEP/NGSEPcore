package ngsep.sequences;
import java.util.Comparator;
import java.lang.Math;

public class SuffixCharSequencePositionComparator implements Comparator<Integer> 
{
	private CharSequence word;
	public SuffixCharSequencePositionComparator(CharSequence word) 
	{
		super();
		this.word = word;
	}
	
	@Override
	public int compare(Integer o1, Integer o2) 
	{
		//o1>o2?
		int r = 0;
		int max = Math.max(o1, o2);
		for (int i = 0; max < word.length(); i++) 
		{
			if(word.charAt(o1)>word.charAt(o2))
			{
				return 1;
			}
			else if(word.charAt(o1)<word.charAt(o2))
			{
				return -1;
			}
			o1++;
			o2++;
			max++;
			if(o1==word.length() & max< word.length())
			{
				return -1;
			}
			if(o2==word.length() & max< word.length())
			{
				return 1;
			}
			//TODO SI SE ACABA 1 ES MENOR....
			
		}
		return r;
	}

}
