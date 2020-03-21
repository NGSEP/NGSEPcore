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

import java.io.Serializable;
import java.lang.reflect.Constructor;
import java.util.Arrays;

/**
 * Implements the Sequence interface minimizing memory consumption no matter 
 * which is the underlying alphabet
 * @author Jorge Duitama
 *
 */
public abstract class AbstractLimitedSequence implements LimitedSequence, Serializable {
	private static final long serialVersionUID = 1L;
	
	int [] sequence= new int [0];
	private byte lastHashSize=0;
	private int length=0;
	
	/**
	 * Makes this sequence equals to the given String
	 * @param sequence Sequence with data to copy
	 */
	public void setSequence(CharSequence sequence) {
		int l = sequence.length();
		int nHashNumbers = calculateHashNumbers(l);
		this.sequence = new int [nHashNumbers];
		if (l==0) lastHashSize = 0;
		else this.lastHashSize = encodeAndAppendSequence(sequence, 0, this.sequence, 0);
		this.length = l;
	}
	
	/**
	 * Append the given String at the end of this
	 * @param sequence String to append
	 */
	public void append(CharSequence sequence) {
		if(sequence.length()==0) return;
		int newLength = length()+sequence.length();
		//if(length()==0 && sequence.length()<5) System.out.println("Appending "+sequence+" to "+this.toString()+". Seq len: "+sequence.length()+" new length: "+newLength);
		int nHashNumbers = calculateHashNumbers(newLength);
		int maxHashSize = getMaxHashSize();
		//Generate junction number
		StringBuffer junction = new StringBuffer();
		if(lastHashSize>0 && lastHashSize < maxHashSize) {
			junction.append(getSequence(this.sequence[this.sequence.length-1], lastHashSize));
		}
		
		int firstPosAppend = 0;
		if(nHashNumbers == this.sequence.length) {
			junction.append(sequence);
			this.sequence[this.sequence.length-1] = getHash(junction.toString(), 0, junction.length());
			this.lastHashSize = (byte)junction.length();
			this.length = newLength;
			return;
		} else if (junction.length()>0) {
			firstPosAppend = maxHashSize-junction.length();
			//Avoid the use of subSequence or substring to avoid infinite recursion
			for(int i=0;i<firstPosAppend;i++) junction.append(sequence.charAt(i));
		}
		int [] newSequence = Arrays.copyOf(this.sequence, nHashNumbers);
		if(junction.length()>0) {
			newSequence[this.sequence.length-1] = getHash(junction.toString(), 0, junction.length());
		}
		this.lastHashSize = encodeAndAppendSequence(sequence, firstPosAppend, newSequence,this.sequence.length);
		
		this.sequence = newSequence;
		this.length = newLength;
	}
	private int calculateHashNumbers(int length) {
		byte maxHashSize = getMaxHashSize();
		int nHashNumbers = length/maxHashSize;
		if(length%maxHashSize>0) nHashNumbers++;
		return nHashNumbers;
	}
	private byte encodeAndAppendSequence(CharSequence sequence, int firstPosAppend, int[] newSequence, int firstIndex) {
		int nHashNumbers = newSequence.length;
		int maxHashSize = getMaxHashSize();
		int j = firstPosAppend;
		for(int i=firstIndex;i<nHashNumbers-1 && j<sequence.length();i++) {
			newSequence[i]= getHash(sequence,j,j+maxHashSize);
			j+=maxHashSize;
		}
		int lastStart = j;
		byte lastHashS = (byte)(sequence.length() - lastStart);
		if(lastHashS>0) newSequence[nHashNumbers-1]= getHash(sequence,lastStart,sequence.length());
		return lastHashS;
	}
	
	/**
	 * Returns the segment of this between the given indexes.
	 * @param beginIndex Zero based beginning index. Inclusive. 0<= beginIndex < this.length()
	 * @param endIndex Zero based ending index. Exclusive. 0<= endIndex <= this.length()
	 * @return Sequence Subsequence between the given indexes. The sequence type should
	 * be the same as this sequence
	 */
	public CharSequence subSequence (int start, int end) {
		if(start < 0 || start >= length) {
			throw new StringIndexOutOfBoundsException(start);
		}
		if(end<0 || end > length) {
			throw new StringIndexOutOfBoundsException(end);
		}
		if(end<start) {
			throw new StringIndexOutOfBoundsException("End index "+end+" cannot be less than start index: "+start);
		}
		AbstractLimitedSequence answer;
		try {
			Constructor<?> emptyConstructor = this.getClass().getConstructor();
			answer = (AbstractLimitedSequence) emptyConstructor.newInstance();
		} catch (Exception e) {
			throw new RuntimeException(e);
		}
		if(start== end) return answer; 
		byte maxHashSize = getMaxHashSize();
		int relStart = start/maxHashSize;
		int relEnd = end/maxHashSize;
		if (end % maxHashSize > 0) {
			relEnd++;
		}
		answer.length = end-start;
		int hashNumbers = calculateHashNumbers(answer.length);
		answer.sequence = new int [hashNumbers];
		int firstIndex = 0;
		StringBuffer buffer = new StringBuffer();
		for(int i=relStart;i<relEnd;i++) {
			int hashSize = getHashSize(i);
			String s = String.valueOf(getSequence(sequence[i],hashSize));
			int firstPos = maxHashSize*i;
			int endPos = maxHashSize*i + hashSize ;
			int firstSPos = 0;
			int endSPos = s.length();
			if(firstPos < start) {
				firstSPos = (start-firstPos);
			}
			if(endPos > end) {
				endSPos = end-firstPos; 
			}
			buffer.append(s.substring(firstSPos, endSPos));
			int bl = buffer.length();
			int blM = bl%maxHashSize;
			int idxR = bl - blM;
			if(buffer.length()>1000000) {
				answer.encodeAndAppendSequence(buffer.substring(0, idxR), 0, answer.sequence, firstIndex);
				firstIndex+=bl/maxHashSize;
				buffer = new StringBuffer(buffer.substring(idxR));
			}
		}
		answer.lastHashSize = answer.encodeAndAppendSequence(buffer, 0, answer.sequence, firstIndex);
		
		return answer;
	}
	public CharSequence subSequence (int start) {
		return subSequence(start,length);
	}
	/**
	 * Returns the character at the given position
	 * @param position Zero based position to look for
	 * @return char Character at the given position
	 */
	public char charAt(int position) {
		if(position < 0 || position>=this.length()) throw new StringIndexOutOfBoundsException(position);
		byte maxHashSize = getMaxHashSize();
		int relPos = position/maxHashSize;
		int subPos = position%maxHashSize;
		int size = getHashSize(relPos);
		char subSeq [] = getSequence(sequence[relPos], size);
		return subSeq[subPos];
	}
	/**
	 * Sets the given character at the given index
	 * @param index Zero based position to modify
	 * @param ch Character to set at index
	 */
	public void setCharAt(int position, char c) {
		if(getAlphabetIndex(c)<0) return;
		if(position < 0 || position>=this.length()) throw new StringIndexOutOfBoundsException(position);
		byte maxHashSize = getMaxHashSize();
		int relPos = position/maxHashSize;
		int subPos = position%maxHashSize;
		int size = getHashSize(relPos);
		char subSeq [] = getSequence(sequence[relPos], size);
		if(subSeq[subPos] != c) {
			subSeq[subPos] = c;
			int hash = getHash(String.valueOf(subSeq), 0, subSeq.length);
			sequence[relPos]=hash;
		}
	}
	@Override
	public int compareTo (LimitedSequence sequence2) {
		if(!(sequence2 instanceof AbstractLimitedSequence)) throw new RuntimeException("Default method only can compare AbstractLimitedSequence objects");
		AbstractLimitedSequence s2 = (AbstractLimitedSequence)sequence2;
		if(this.getAlphabet()!=s2.getAlphabet()) throw new RuntimeException("Attempt to compare two limited sequences with different alphabets");
		int [] seq1 = this.sequence;
		int [] seq2 = s2.sequence;
		int i;
		for(i=0;i<seq1.length -1 && i < seq2.length -1;i++) {
			if(seq1[i]!=seq2[i]) return seq1[i]-seq2[i];
		}
		int hs1 = this.getHashSize(i);
		int hs2 = s2.getHashSize(i);
		char [] c1 = getSequence(seq1[i], hs1);
		char [] c2 = getSequence(seq2[i], hs2);
		for(int j=0;j<c1.length && j<c2.length;j++) {
			if(c1[j]!=c2[j]) return getAlphabetIndex(c1[j])-getAlphabetIndex(c2[j]);
		}
		if(c1.length!=c2.length) return c1.length - c2.length;
		if(seq1.length!=seq2.length) return seq1.length - seq2.length;
		return 0;
	}
	@Override
	public String toString() {
		if(length()==0) return "";
		StringBuffer answer = new StringBuffer();
		for(int i=0;i<sequence.length;i++) {
			int size = getHashSize(i);
			answer.append(getSequence(sequence[i], size));
		}
		return answer.toString();
	}
	
	/* (non-Javadoc)
	 * @see java.lang.Object#equals(java.lang.Object)
	 */
	@Override
	public boolean equals(Object o) {
		if(!(o instanceof AbstractLimitedSequence)) return false;
		LimitedSequence seq2 = (LimitedSequence) o;
 		return compareTo(seq2)==0;
	}

	/* (non-Javadoc)
	 * @see java.lang.Object#hashCode()
	 */
	@Override
	public int hashCode() {
		if(sequence.length==0) return 0;
		long absoluteNumber = (long)sequence[0]-(long)Integer.MIN_VALUE;
		return ((int)(absoluteNumber%100000000));
	}

	private int getHashSize(int pos) {
		byte maxHashSize = getMaxHashSize();
		if(pos < sequence.length-1) {
			return maxHashSize;
		} else if (pos == sequence.length-1) {
			return this.lastHashSize;
		}
		return 0;
	}
	
	@Override
	public int length() {
		return length;
	}
	private byte maxHashSize = 0;
	/**
	 * Return the maximum number of characters encoded in a single number
	 * @return byte
	 */
	private byte getMaxHashSize() {
		if(maxHashSize==0) maxHashSize = (byte)(32/getBitsPerCharacter());
		return maxHashSize;
	}
	/**
	 * Returns the number corresponding with a suitable size substring of the given sequence 
	 * @param seq Sequence to calculate hash
	 * @param start Zero based first position
	 * @param end Zero based last position
	 * @return int Number representing the substring of seq between start (included) and end (not included)
	 */
	private int getHash(CharSequence seq, int start, int end) {
		long number = getHash(seq, start, end, this);
		return (int)(number+Integer.MIN_VALUE);
	}
	/**
	 * Gets the sequence corresponding with the given hash and the given size
	 * @param number Hash number to decode
	 * @param size Size of the sequence to return
	 * @return char[] Decoded String as a char array
	 */
	private char [] getSequence(int number, int size) {
		long absoluteNumber = (long)number-(long)Integer.MIN_VALUE;
		return getSequence(absoluteNumber, size,this);
	}
	
	/**
	 * Returns the number corresponding with a suitable size substring of the given sequence 
	 * @param seq Sequence to calculate hash
	 * @param start Zero based start position
	 * @param end Zero based end position
	 * @param targetSeq AbstractLimitedSequence having the target alphabet
	 * @return long Positive number representing the substring of seq between start (included) and end (not included)
	 */
	public static long getHash(CharSequence seq, int start, int end, AbstractLimitedSequence targetSeq) {
		long number =0;
		int alpSize = targetSeq.getAlphabetSize();
		for(int i=start;i<end;i++) {
			number*=alpSize;
			int index = targetSeq.getAlphabetIndex(seq.charAt(i));
			if(index <0) {
				index = targetSeq.getDefaultIndex();
			}
			if(index <0) {
				throw new IllegalArgumentException("Character "+seq.charAt(i)+" not supported by sequence of type "+targetSeq.getClass().getName());
			}
			number+=index;
			if(number<0) throw new RuntimeException("Encoding reached a long negative number for sequence: "+seq+" between "+start+" and "+end);
		}
		return number;
	}
	
	/**
	 * Gets the sequence corresponding with the given hash and the given size
	 * @param number Hash number to decode
	 * @param size Size of the sequence to return
	 * @param targetSeq AbstractLimitedSequence having the target alphabet
	 * @return char[] Decoded String as a char array
	 */
	public static char [] getSequence(long number, int size, AbstractLimitedSequence targetSeq) {
		char [] answer = new char[size];
		int alpSize = targetSeq.getAlphabetSize();
		for(int i=0;i<size;i++) {
			int nextDigit = (int)(number%alpSize);
			int index = size-i-1;
			answer[index] = targetSeq.getAlphabetCharacter(nextDigit);
			number = number/alpSize;
		}
		return answer;
	}
	
	/**
	 * Returns the length of the maximum overlap between a suffix of sequence 1 and a prefix of sequence 2
	 * @param sequence1 Sequence to evaluate suffixes
	 * @param sequence2 Sequence to evaluate prefixes
	 * @return int Maximum overlap between a prefix of sequence2 and a suffix of sequence 1
	 */
	public static int getOverlapLength(CharSequence sequence1, CharSequence sequence2) {
		for(int i=0;i<sequence1.length();i++) {
			if(isSuffixAPrefix(sequence1,i,sequence2)) return sequence1.length()-i;
		}
		return 0;
	}

	private static boolean isSuffixAPrefix(CharSequence sequence1, int start1, CharSequence sequence2) {
		int i = start1;
		for(int j=0;i<sequence1.length() && j<sequence2.length();j++) {
			if(sequence1.charAt(i)!=sequence2.charAt(j)) return false;
			i++;
		}
		return i==sequence1.length();
	}
	/**
	 * Returns the character of the alphabet with the given index in the alphabet String.
	 * @param index Zero based index to look for in the alphabet
	 * @return char Character at the given index
	 */
	public char getAlphabetCharacter(int index) {
		String alp = getAlphabet();
		if(index>=0 && index<alp.length()) {
			return alp.charAt(index);
		}
		return 0;
	}
	/**
	 * Gets the index for the given character in the alphabet.
	 * Default implementation performs a linear search on the alphabet
	 * @param ch character to look for
	 * @return in Zero based index or -1 if the character is not in the alphabet
	 */
	public int getAlphabetIndex(char ch) {
		String alp = getAlphabet();
		return alp.indexOf(ch);
	}
	/**
	 * Tells if the character is in the alphabet
	 * @param ch character to look for
	 * @return boolean true if ch is in the alphabet, false otherwise
	 */
	public boolean isInAlphabet(char ch) {
		return getAlphabetIndex(ch)>=0;
	}
	/**
	 * Returns the size of the alphabet
	 * @return int Size of the alphabet
	 */
	public int getAlphabetSize() {
		return getAlphabet().length();
	}
	/**
	 * Gets the minimum number of bits needed to represent uniquely every character of
	 * the alphabet 
	 * @return int minimum number of bits needed to represent uniquely every character of
	 * the alphabet
	 */
	protected int getBitsPerCharacter() {
		return (int)(Math.log10(getAlphabetSize())/Math.log10(2))+1;
	}
	
	protected int getDefaultIndex () {
		return -1;
	}
	
}
