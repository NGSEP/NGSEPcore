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
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Random;
import java.util.TreeSet;
/**
 * 
 * @author Jorge Duitama
 *
 */
public class DNAShortKmer implements CharSequence, Comparable<DNAShortKmer>, Serializable {
	private static final long serialVersionUID = 1L;
	private long index;
	private byte length;
	private static final DNASequence EMPTYDNASEQ = new DNASequence();
	
	public DNAShortKmer(CharSequence kmerSeq) {
		length = (byte) kmerSeq.length();
		if(length>31) throw new IllegalArgumentException("The maximum k-mer size for this class is 31. Input k-mer: "+kmerSeq);
		index = AbstractLimitedSequence.getHash(kmerSeq, 0, length,EMPTYDNASEQ);
	}
	
	@Override
	public char charAt(int i) {
		char [] characters = AbstractLimitedSequence.getSequence(index, length, EMPTYDNASEQ);
		return characters[i];
	}
	@Override
	public int length() {
		return length;
	}
	@Override
	public CharSequence subSequence(int start, int end) {
		String s = toString();
		return s.subSequence(start, end);
	}
	
	
	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		return new String (AbstractLimitedSequence.getSequence(index, length, EMPTYDNASEQ));
	}

	/* (non-Javadoc)
	 * @see java.lang.Object#equals(java.lang.Object)
	 */
	@Override
	public boolean equals(Object o) {
		if(!(o instanceof DNAShortKmer)) return false;
		DNAShortKmer kmer2 = (DNAShortKmer) o;
 		return length== kmer2.length && index==kmer2.index; 
	}

	/* (non-Javadoc)
	 * @see java.lang.Object#hashCode()
	 */
	@Override
	public int hashCode() {
		return (int) (index%100000000);
	}

	@Override
	public int compareTo(DNAShortKmer kmer2) {
		if(length!=kmer2.length) return length - kmer2.length;
		if(index == kmer2.index) return 0;
		return (index<kmer2.index)?-1:1;
	}
	/**
	 * Main method to test the number of short k-mers that can be added to a list
	 * @param args
	 */
	public static void main(String[] args) throws Exception {
		int numKmers = Integer.parseInt(args[0]);
		int length = Integer.parseInt(args[1]);
		int collectiontype = Integer.parseInt(args[2]);
		Collection<DNAShortKmer> kmers = new ArrayList<>();
		boolean isList = true;
		if (collectiontype == 1) {
			kmers = new LinkedList<>();
		} else if (collectiontype == 2) {
			kmers = new TreeSet<>();
			isList = false;
		} else if (collectiontype == 3) {
			kmers = new HashSet<>();
			isList = false;
		}
		
		char [] kmerC = new char [length];
		Random r = new Random();
		long time1 = System.currentTimeMillis();
		for(int i=0;i<numKmers;i++) {
			//Create kmer
			for(int j=0;j<kmerC.length;j++) {
				int bpI = r.nextInt(4);
				kmerC[j] = DNASequence.BASES_STRING.charAt(bpI);
			}
			String kmerStr = new String(kmerC);
			DNAShortKmer kmer = new DNAShortKmer(kmerStr);
			kmers.add(kmer);
			if((i+1)%100000==0) System.out.println("Created "+kmers.size()+" kmers");
		}
		long time2 = System.currentTimeMillis();
		System.out.println("Time creating k-mers: "+(time2-time1)+" Number of k-mers: "+kmers.size());
		time1 = time2;
		if(isList) {
			Collections.sort(((List<DNAShortKmer>)kmers));
			time2 = System.currentTimeMillis();
			System.out.println("Time sorting k-mers: "+(time2-time1));
			time1 = time2;
		}
		//Search experiment
		for(int i=0;i<1000000;i++) {
			//Create kmer
			for(int j=0;j<kmerC.length;j++) {
				int bpI = r.nextInt(4);
				kmerC[j] = DNASequence.BASES_STRING.charAt(bpI);
			}
			String kmerStr = new String(kmerC);
			DNAShortKmer kmer = new DNAShortKmer(kmerStr);
			if(isList) {
				Collections.binarySearch((List<DNAShortKmer>)kmers, kmer);
			} else {
				kmers.contains(kmer);
			}
		}
		time2 = System.currentTimeMillis();
		System.out.println("Time searching 1000000 k-mers: "+(time2-time1));
		time1 = time2;
	}
}
