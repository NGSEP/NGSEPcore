package ngsep.kmerscount;

public class ReverseComplement {

	String sequence="";
	
	public ReverseComplement(String sequence){
		this.sequence = sequence;
	}
	
	public String makeReverseComplement(){
		
		String reverseComplementSequence = "";
		int seqLength = sequence.length();
		
		for(int i = 0; i < seqLength; i++)
		{
			char nucleotide = sequence.charAt(i);
			
			switch(nucleotide){
			case 'A': 
				nucleotide = 'T';
				break;
			case 'T': 
				nucleotide = 'A';
				break;
			case 'G': 
				nucleotide = 'C';
				break;
			case 'C': 
				nucleotide = 'G';
				break;
			case 'a': 
				nucleotide = 't';
				break;
			case 't': 
				nucleotide = 'a';
				break;
			case 'g': 
				nucleotide = 'c';
				break;
			case 'c': 
				nucleotide = 'g';
				break;
			}
			reverseComplementSequence = nucleotide + reverseComplementSequence;
		}
		
		return reverseComplementSequence;
	}
	
}
