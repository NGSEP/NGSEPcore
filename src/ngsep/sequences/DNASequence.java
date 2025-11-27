/*******************************************************************************
 * NGSEP - Next Generation Sequencing Experience Platform
 * Copyright 2016 Jorge Duitama
 *
 * This file is part of NGSEP.
 *
 *     NGSEP is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     NGSEP is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with NGSEP.  If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/
package ngsep.sequences;

import java.util.Random;

/**
 * Implementation of AbstractSequence for DNA Sequences without degenerate bases
 * @author Jorge Duitama
 */
public class DNASequence extends AbstractLimitedSequence{
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	public static final String [] BASES_ARRAY = {"A","C","G","T"};
	public static final String BASES_STRING = "ACGT";
	private static final int [] ARRAY_BASES_INDEXING = {0,-1,1,-1,-1,-1,2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,3,-1,-1,-1,-1,-1,-1};
	public static final DNASequence EMPTY_DNA_SEQUENCE = new DNASequence();
	
	/**
	 * Creates a new DNA sequence with the given bases
	 * @param sequence CharSequence with the bases for the new sequence
	 */
	public DNASequence(CharSequence sequence) {
		this.setSequence(sequence);
	}
	/**
	 * Creates a new empty DNA Sequence
	 */
	public DNASequence() {
	}
	@Override
	public String getAlphabet() {
		return BASES_STRING;
	}
	@Override
	protected int getBitsPerCharacter() {
		return 2;
	}
	
	@Override
	protected byte getMaxHashSize() {
		return (byte)16;
	}
	private static int getDNAIndex (char base) {
		int i = base - 65;
		try {
			return ARRAY_BASES_INDEXING[i];
		} catch (ArrayIndexOutOfBoundsException e) {
			return -1;
		}
	}
	@Override
	public boolean isInAlphabet(char base) {
		return getDNAIndex(base)>=0;
	}
	
	public static boolean isDNA (CharSequence sequence) {
		for(int i=0;i<sequence.length();i++) {
			if(getDNAIndex(sequence.charAt(i))<0) return false;
		}
		return true;
	}
	
	
	/**
	 * @return DNASequence Reverse complement of this sequence
	 */
	public DNASequence getReverseComplement() {
		StringBuffer sb = new StringBuffer();
		DNASequence answer = new DNASequence();
		int l = length();
		for(int i=l-1;i>=0;i--) {
			sb.append(getComplement(charAt(i)));
			if(sb.length()%100000==0) {
				answer.append(sb);
				sb =new StringBuffer();
			}
		}
		answer.append(sb);
		return answer;
	}
	/**
	 * Gets the complement of the given base
	 * @param base to complement
	 * @return char Complement of the given base
	 */
	public static char getComplement(char base) {
		int index = getDNAIndex(base);
		if(index>=0) {
			return BASES_STRING.charAt(BASES_STRING.length()-index-1);
		}
		return base;
	}
	/**
	 * Returns the reverse complement of a DNA sequence encoded in the given String
	 * @param sequence Sequence to be translated
	 * @return CharSequence reverse complement of the given sequence
	 */
	public static CharSequence getReverseComplement(CharSequence sequence) {
		if (sequence instanceof DNASequence) {
			return ((DNASequence)sequence).getReverseComplement();
		}
		return new DNASequence(sequence).getReverseComplement();
	}
	
	/**
	 * Returns the number corresponding with a suitable size substring of the given sequence 
	 * @param seq Sequence to calculate hash
	 * @param start Zero based start position
	 * @param end Zero based end position
	 * @param targetSeq AbstractLimitedSequence having the target alphabet
	 * @return long Positive number representing the substring of seq between start (included) and end (not included)
	 */
	@Override
	public long getLongCode(CharSequence seq, int start, int end) {
		long number =0;
		for(int i=start;i<end;i++) {
			number <<= 2;
			int index = getDNAIndex(seq.charAt(i));
			if(index <0) {
				throw new IllegalArgumentException("Character "+seq.charAt(i)+" not supported by "+DNASequence.class.getName());
			}
			number+=index;
			if(number<0) throw new RuntimeException("Encoding reached a long negative number for sequence: "+seq+" between "+start+" and "+end);
		}
		return number;
	}
	
	/**
	 * Gets the sequence corresponding with the given code and the given length
	 * @param code Number to decode
	 * @param length of the sequence to return
	 * @return char[] Decoded String as a char array
	 */
	@Override
	public char [] getSequenceFromCode(long number, int length) {
		char [] answer = new char[length];
		for(int i=0;i<length;i++) {
			int nextDigit = (int)(number&3);
			int index = length-i-1;
			answer[index] = DNASequence.BASES_STRING.charAt(nextDigit);
			number = number>>2;
		}
		return answer;
	}
	/**
	 * Calculates the code corresponding to the sequence removing the first character of the
	 * sequence represented by the given code and adding the given character at the end
	 * @param code of the original short sequence
	 * @param length of the original and of the new sequence
	 * @param nextCharacter to append at the end of the new sequence
	 * @return long encoding of the new sequence
	 */
	@Override
	public long getNextCode(long code, int length, char nextCharacter) {
		int index = getDNAIndex(nextCharacter);
		if(index <0) {
			throw new IllegalArgumentException("Character "+nextCharacter+" not supported by sequence of type "+DNASequence.class.getName());
		}
		long bitModule = 1L << (2*(length-1));
		bitModule--;
		long answer = code & bitModule;
		answer <<= 2;
		answer+=index;
		return answer;
	}
	/**
	 * Test main for methods of AbstractLimitedSequence
	 * @param args
	 */
	public static void main(String[] args) {
		StringBuilder randomSequence = new StringBuilder();
		Random r = new Random();
		int length = r.nextInt(10000)+100000;
		for(int i=0;i<length;i++) {
			int bpI = r.nextInt(4);
			randomSequence.append(DNASequence.BASES_STRING.charAt(bpI));
		}
		DNASequence dna = new DNASequence(randomSequence);
		if(!randomSequence.toString().equals(dna.toString())) throw new RuntimeException("Sequences not equal");
		//System.out.println("Sequence: "+ randomSequence);
		for(int i=0;i<randomSequence.length();i++) {
			if(randomSequence.charAt(i)!=dna.charAt(i)) {
				throw new RuntimeException("Error with charAt around position: "+i+" . Sequence: "+randomSequence.substring(Math.max(0, i-5),Math.min(randomSequence.length(), i+5))+" expected "+randomSequence.charAt(i)+" given "+dna.charAt(i));
			}
			char mutSeq = DNASequence.BASES_STRING.charAt(r.nextInt(4));
			randomSequence.setCharAt(i, mutSeq);
			dna.setCharAt(i, mutSeq);
		}
		if(!randomSequence.toString().equals(dna.toString())) throw new RuntimeException("Mutated sequences not equal");
		System.out.println("charAt and setCharAt test finished");
		for(int i=0;i<100;i++) {
			int first = r.nextInt(length);
			int end = Math.min(length, r.nextInt(length)+first);
			String subseq1 = randomSequence.substring(first, end);
			DNASequence subseq2DNA = (DNASequence) dna.subSequence(first, end);
			String subseq2 = subseq2DNA.toString();
			if(!subseq1.equals(subseq2)) throw new RuntimeException("Substrings from "+first+" to "+end+" not equal. expected: "+subseq1+ " given "+subseq2);
			randomSequence.append(subseq1);
			dna.append(subseq1);
			
			if(!randomSequence.toString().equals(dna.toString())) throw new RuntimeException("Appended strings are not equal expected: "+randomSequence.toString()+ " given "+dna);
		}
		System.out.println("subSequence and append test finished");
		//Efficiency
		String randomSeq = dna.getReverseComplement().toString();
		long time1 = System.currentTimeMillis();
		randomSequence = new StringBuilder(randomSeq);
		long time2 = System.currentTimeMillis();
		System.out.println("Time StringBuilder construction: "+ (time2 - time1));
		time1 = time2;
		dna = new DNASequence(randomSeq);
		time2 = System.currentTimeMillis();
		System.out.println("Time DNASequence construction: "+ (time2 - time1));
		time1 = time2;
		for(int i=0;i<randomSequence.length();i++) randomSequence.charAt(i);
		time2 = System.currentTimeMillis();
		System.out.println("Time StringBuilder charAt: "+ (time2 - time1));
		time1 = time2;
		for(int i=0;i<randomSequence.length();i++) dna.charAt(i);
		time2 = System.currentTimeMillis();
		System.out.println("Time DNASequence charAt: "+ (time2 - time1));
		time1 = time2;
		for(int i=0;i<100;i++) {
			int first = r.nextInt(length);
			int end = Math.min(length, r.nextInt(length)+first);
			randomSequence.subSequence(first, end);
		}
		time2 = System.currentTimeMillis();
		System.out.println("Time StringBuilder substring: "+ (time2 - time1));
		time1 = time2;
		for(int i=0;i<100;i++) {
			int first = r.nextInt(length);
			int end = Math.min(length, r.nextInt(length)+first);
			dna.subSequence(first, end);
		}
		time2 = System.currentTimeMillis();
		System.out.println("Time DNASequence substring: "+ (time2 - time1));
		int first = r.nextInt(length);
		int end = Math.min(length, r.nextInt(length)+first);
		String seq = randomSequence.substring(first, end);
		time1 = time2;
		for(int i=0;i<100;i++) {
			randomSequence.append(seq);
		}
		time2 = System.currentTimeMillis();
		System.out.println("Time StringBuilder append: "+ (time2 - time1));
		time1 = time2;
		for(int i=0;i<100;i++) {
			dna.append(seq);
		}
		time2 = System.currentTimeMillis();
		System.out.println("Time DNASequence append: "+ (time2 - time1));
		time1 = time2;
		long code = DNASequence.EMPTY_DNA_SEQUENCE.getLongCode(randomSequence, 0, 15);
		for(int i=15;i<randomSequence.length();i++) {
			code = DNASequence.EMPTY_DNA_SEQUENCE.getNextCode(code, 15, randomSequence.charAt(i));
			if(code!= (code & 0x3FFFFFFF)) System.out.println("Code "+code+" larger than 2E30");
			char [] seq2 = DNASequence.EMPTY_DNA_SEQUENCE.getSequenceFromCode(code, 15);
			String k1 = new String (seq2);
			String k2 = randomSequence.substring(i-14, i+1);
			if(!k1.equals(k2)) System.out.println("Error encoding / decoding next character. Expected: "+k2+" given: "+k1); 
		}
		time2 = System.currentTimeMillis();
		System.out.println("Time kmers encoding and decoding: "+ (time2 - time1));
	}
}
