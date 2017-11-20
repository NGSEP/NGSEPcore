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
	private String sequenceName;
	private Map<Integer,Integer> partialSuffixArray = new HashMap<>();
	private List<Integer []> tallyIndexes = new ArrayList<>();
	private int tallyDistance=50;
	private int sufixDistance=30;
	private String bwt;
	private int [] characterCounts;
	private String alphabet="";
	private String sequence;
	
	public FMIndexSingleSequence(String seqName, CharSequence sequence) 
	{
		sequenceName=seqName;
		this.sequence=sequence.toString();
		
		ArrayList<Character> alphabetArrayList = new ArrayList<>();
		int n = this.sequence.length();
		for (int j = 0; j < n; j++) 
		{
			Character a = sequence.charAt(j);
			if(!alphabetArrayList.contains(a))
			{
				alphabetArrayList.add(a);
			}
		}
		
		Collections.sort(alphabetArrayList, new Comparator<Character>() 
		{

			@Override
			public int compare(Character arg0, Character arg1) 
			{
				if(arg0==arg1)
					return 0;
				if(arg0>arg1)
				return 1;
				return -1;
			}
		});
		
		StringBuilder alphabetBuilder = new StringBuilder();
		
		for (int i = 0; i < alphabetArrayList.size(); i++) 
		{
			alphabetBuilder.append(alphabetArrayList.get(i));
		}
		alphabet=alphabetBuilder.toString();
		this.sequence+="$";
		
		calculate();
	}
	private void calculate() 
	{
		ArrayList<Integer> sufixes = new ArrayList<Integer>();
		for (int i = 0; i < sequence.length(); i++) 
		{
			sufixes.add(i);
		}
		Collections.sort(sufixes, new ComparatorString(sequence));
		String[] lAndF = getLastAndFirst(sufixes);
//		System.out.println("l "+lAndF[0]);
//		System.out.println("f "+lAndF[1]);
//		printSuffixes();
		buildTally(sufixes);
		bwt=lAndF[0];
		createPartialSufixArray(sufixes);
		
		ArrayList<Character> vistos = new ArrayList<>();
		characterCounts=new int[alphabet.length()];
		int j =-1;
		for (int i = 1; i < lAndF[1].length(); i++) 
		{
			Character c = lAndF[1].charAt(i);
			if(!vistos.contains(c))
			{
				j++;
				if(j<alphabet.length())
				characterCounts[j]++;
				vistos.add(c);
			}
			else if(j<alphabet.length())
				characterCounts[j]++;
		}
		
//		System.out.println("f comprex:");
		for (int i = 0; i < alphabet.length(); i++) 
		{
//			System.out.println(alphabet.charAt(i)+" "+characterCounts[i]);
		}
	}
	private void createPartialSufixArray(ArrayList<Integer> sufixes) 
	{
		partialSuffixArray = new HashMap<Integer,Integer>();
		
		int n = sequence.length()-1;
		
		while(n>=0)
		{
			partialSuffixArray.put(sufixes.indexOf(n), n);
			n-=sufixDistance;
		}
		partialSuffixArray.put(sufixes.indexOf(0), 0);
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
		
//		System.out.println("buscando: "+searchSequence);
		ArrayList<Integer> range = getRange(searchSequence);

		for (int i = 0; i < range.size(); i++) 
		{
			int begin=0;
			int actual = range.get(i);
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
		
		return alignments;
	}
	private ArrayList<Integer> getRange(String searchSequence) 
	{
		int begin =0;
		int fin = bwt.length()-1;
		boolean parar = false;
//		System.out.println(wantedWord);
		for (int i = searchSequence.length()-1 ; i >=0 && !parar; i--) 
		{
			Character actual = searchSequence.charAt(i);
			if(!alphabet.contains(""+actual))
			{
				return new ArrayList<Integer>();
			}
//			System.out.println(actual);
			//			fin = fin>=f.length()?f.length()-1:fin;

			int finCopy = fin;
			for (int j = begin; j <= finCopy; j++) 
			{
				//				char actual2 = f.charAt(j);
				char actual2 = getCharacterOfFirstAt(j);
				if(actual2==actual)
				{
					if(begin==0)
						begin=j;

					fin=j;
				}
			}
//			System.out.println("[ "+inicio+", "+fin+"]");

			if(i-1>=0)
			{
				char siguiente = searchSequence.charAt(i-1);
				int tallyI=getTallyOf(siguiente, begin-1);
				int tallyF=getTallyOf(siguiente, fin);

				int diferencia=tallyF-tallyI;
				if(diferencia==0)
				{
//					System.out.println("no");
					return new ArrayList<Integer>();
//					throw new Exception("no existe");
				}

				//				if(f.indexOf(siguiente)+tallyI+tallyF-1>=f.length())
				//Acá está el error
				//en el calculo de los tally hay error
				if(getIndexInFirstOf(siguiente)+tallyI-1>=bwt.length())
				{
					parar=true;
				}
				else
				{

					//					inicio=f.indexOf(siguiente)+tallyI;
					begin=getIndexInFirstOf(siguiente)+tallyI;
					fin=begin+diferencia-1;
				}
			}
		}
		ArrayList<Integer> r = new ArrayList<Integer>();
		int finCopy = fin;
		for (int j = begin; j <= finCopy; j++) 
		{
			r.add(j);
		}
		return r;
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

		if( i-a*tallyDistance < Math.abs(i-b*tallyDistance)) 
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
