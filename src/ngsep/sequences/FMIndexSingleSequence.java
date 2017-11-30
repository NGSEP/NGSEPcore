package ngsep.sequences;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import ngsep.genome.GenomicRegion;
import ngsep.genome.GenomicRegionImpl;

public class FMIndexSingleSequence implements Serializable 
{
	//Name of the sequence
	private String sequenceName;

	//Original sequence
	private String sequence;

	//startSeq of some indexes  representing a partial suffix array
	private Map<Integer,Integer> partialSuffixArray = new HashMap<>();

	//Ranks in the bwt for each character in the alphabet
	private List<Integer []> tallyIndexes = new ArrayList<>();

	//1 of each tallyDistance is saved
	private int tallyDistance=50;

	// 1/suffixFraction indexes are saved
	private int suffixFraction=10;

	//Burrows Wheeler transform
	private String bwt;

	//times each character appears
	private int [] characterCounts;

	//Inferred alphabet of the sequence ordered lexicographical 
	private String alphabet="";

	public FMIndexSingleSequence(String seqName, CharSequence sequence) 
	{
		sequenceName=seqName;
		this.sequence=sequence.toString();
		this.sequence+="$";

		calculate();
	}

	private void calculate() 
	{
		ArrayList<Integer> sufixes = buildSuffixArray();
		String[] lAndF = getLastAndFirst(sufixes);
		//		System.out.println("l "+lAndF[0]);
		//		System.out.println("f "+lAndF[1]);
		//		printSuffixes();
		bwt=lAndF[0];

		ArrayList<Character> seen = new ArrayList<>();
		ArrayList<Integer> counts= new ArrayList<>();

		//iterate last column to know alphabet and counts...
		int j =-1;
		for (int i = 1; i < lAndF[1].length(); i++) 
		{
			Character c = lAndF[1].charAt(i);
			if(!seen.contains(c))
			{
				j++;
				counts.add(1);
				seen.add(c);
			}
			else
				counts.set(j, counts.get(j)+1);
		}

		buildAlphabet(seen);
		buildCharacterCounts(counts);
		buildTally(sufixes);
		createPartialSufixArray(sufixes);

	}

	private void buildCharacterCounts(ArrayList<Integer> counts) 
	{
		characterCounts=new int[alphabet.length()];
		for (int i = 0; i < counts.size(); i++) 
		{
			characterCounts[i]=counts.get(i);
		}
	}

	private void buildAlphabet(ArrayList<Character> seen) 
	{
		StringBuilder alp=new StringBuilder();
		for (int i = 0; i < seen.size(); i++) 
		{
			alp.append(seen.get(i));
		}
		alphabet=alp.toString();
	}
	private ArrayList<Integer> buildSuffixArray() {
		ArrayList<Integer> sufixes = new ArrayList<Integer>();
		for (int i = 0; i < sequence.length(); i++) 
		{
			sufixes.add(i);
		}
		Collections.sort(sufixes, new SuffixCharSequencePositionComparator(sequence));
		return sufixes;
	}
	private void createPartialSufixArray(ArrayList<Integer> sufixes) 
	{
		partialSuffixArray = new HashMap<Integer,Integer>();
		int n = sufixes.size();
		for(int i=0;i<n;i++) 
		{
			int startSeq = sufixes.get(i);
			if(startSeq%suffixFraction==0) 
			{
				partialSuffixArray.put(i, startSeq);
			}
		}
	}
	private void buildTally(ArrayList<Integer> sufixes) 
	{
		Integer[] arr= new Integer[alphabet.length()];
		for (int i = 0; i < arr.length; i++) 
		{
			arr[i]=0;
		}

		for (int i = 0; i < sufixes.size(); i++) 
		{

			int n = sufixes.get(i)-1;
			if(n>=0)
			{
				Character actual = sequence.charAt(n);
				//				System.out.println(actual);
				arr[alphabet.indexOf(actual)]++;
			}

			if(i%tallyDistance==0)
			{
				Integer[] copy= new Integer[alphabet.length()];
				for (int j = 0; j < copy.length; j++) 
				{
					copy[j]=arr[j];
				}
				tallyIndexes.add(copy);
			}
		}
	}
	private String[] getLastAndFirst(ArrayList<Integer> sufixes) 
	{
		String[] lAndF = new String[2];

		StringBuilder bwt = new StringBuilder();
		StringBuilder f = new StringBuilder();

		for (int i = 0; i < sufixes.size(); i++) 
		{
			int n = sufixes.get(i)-1;
			if(n>=0)
			{
				//				System.out.println(word.charAt(n));
				bwt.append(sequence.charAt(n));
				Character c = sequence.charAt(n+1);
				f.append(c);
			}
			else
			{
				bwt.append("$");
				f.append(sequence.charAt(0));
				//				System.out.println("$");	
			}
		}
		lAndF[0]=bwt.toString();
		lAndF[1]=f.toString();

		return lAndF;
	}
	public FMIndexSingleSequence(QualifiedSequence sequence) 
	{
		this (sequence.getName(),sequence.getCharacters());
	}
	public List<GenomicRegion> search (String searchSequence) 
	{
		List<GenomicRegion> alignments = new ArrayList<>();
		int[] range = getRange(searchSequence);
		if (range !=null)
		{

			for (int i = range[0]; i <= range[1]; i++) 
			{
				int begin=0;
				int actual = i;
				if(partialSuffixArray.containsKey(actual))
				{
					begin = partialSuffixArray.get(actual);
				}
				else
				{
					boolean encontrado = false;
					int posible = 0;
					int steps =0;
					while(!encontrado)
					{
						//					System.out.println("ac " +actual);
						posible = lfMap(actual);
						encontrado=partialSuffixArray.containsKey(posible);
						actual =posible;
						steps++;
					}
					begin = partialSuffixArray.get(posible)+steps;
				}
				alignments.add(new GenomicRegionImpl(sequenceName,begin , begin+searchSequence.length()));

			}
		}

		return alignments;
	}

	public int[] getRange(String query)
	{
		char c = query.charAt(query.length()-1);
		
		int idxS=1;
		int idxF=-1;
		for(int i=0;i<alphabet.length();i++) 
		{
			char ct = alphabet.charAt(i);
			if(ct!=c) idxS += characterCounts[alphabet.indexOf(ct)];
			else 
			{
				idxF = idxS + characterCounts[alphabet.indexOf(ct)] - 1 ;
				break;
			}
		}
		if(idxF==-1) {
			return null;
		}
		for(int j=query.length()-2;j>=0;j--) {
			c = query.charAt(j);
			if(!alphabet.contains(""+c)) {
				return null;
			}
			boolean add1 = (bwt.charAt(idxS)!=c); 
			idxS = lfMapping(c, idxS);
			if(add1) idxS++; 
			idxF = lfMapping(c, idxF);
			if(idxS>idxF) {
				return null;
			}
		}

		return new int[] {idxS, idxF };
	}
	
	
	
	private Character getCharacterOfFirstAt(int j) 
	{
		if(j==0)
			return '$';
		int cuenta=0;
		for (int i = 0; i < characterCounts.length; i++) 
		{
			cuenta+=characterCounts[i];
			if(j<=cuenta)
				return alphabet.charAt(i);
		}
		return null;
	}

	private int getIndexInFirstOf(char c) 
	{
		int indexOfChar = alphabet.indexOf(c);
		int sum = 0;
		for (int i = 0; i < indexOfChar; i++) 
		{
			sum+=characterCounts[i];
		}
		return sum+1;
	}
	public int getTallyOf(Character c, int i)
	{
		int r = 0;

		int a = i/tallyDistance;
		int b = a+1;

		if( i-a*tallyDistance < Math.abs(i-b*tallyDistance) || tallyIndexes.size()<=b) 
		{

			int d = tallyIndexes.get(a)[alphabet.indexOf(c)];

			for (int j = a*tallyDistance+1; j <= i; j++) 
			{
				//				Character cA = l.charAt(j);
				Character cA = getCharacterOfBWTAt(j);
				if(cA==c)
					d++;
			}
			r=d;
		}
		else
		{
			int d = tallyIndexes.get(b)[alphabet.indexOf(c)];

			for (int j = b*tallyDistance; j > i; j--) 
			{
				Character cA = getCharacterOfBWTAt(j);
				if(cA==c)
					d--;
			}
			r=d;
		}


		return r;
	}
	private Character getCharacterOfBWTAt(int j) 
	{
		//		System.out.println(j);
		return j<bwt.length()?bwt.charAt(j):'$';
	}
	
	private int lfMapping(char c, int pos) 
	{
		int rank = getTallyOf(c, pos);
		int idx = 1;
		for(int i=0;i<alphabet.length();i++) 
		{
			char ct = alphabet.charAt(i); 
			if(ct!=c) idx += characterCounts[alphabet.indexOf(ct)];
			else {
				idx += (rank-1);
				break;
			}
		}
		return idx;
	}
	
	int lfMap (int indexOfBwt)
	{
		char c = getCharacterOfBWTAt(indexOfBwt);
		//		System.out.println(""+c);
		int count =getTallyOf(c, indexOfBwt)-1;
		return getIndexInFirstOf(c)+count;
	}
	public String getSequenceName() 
	{
		return sequenceName;
	}
	public int getTallyDistance() {
		return tallyDistance;
	}
	public void setTallyDistance(int tallyDistance) {
		this.tallyDistance = tallyDistance;
	}
	public String getAlphabet() {
		return alphabet;
	}
	public void setAlphabet(String alphabet) {
		this.alphabet = alphabet;
	}
	
	

}