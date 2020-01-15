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

import java.io.PrintStream;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Map.Entry;

import ngsep.math.Distribution;

/**
 * @author Jorge Gomez
 * @author Jorge Duitama
 */
public class DNAShortKmerClusterMap implements KmersMap {
	private int kmerLength;
	private int maxNumClusters;
	private int[][] table;
	private Map<DNAShortKmer, Integer> index = new HashMap<>();
	private int newIndex = 0;
	private int numClusters = 0;
	
	public DNAShortKmerClusterMap (int kmerLength, int maxNumClusters) {
		this.kmerLength = kmerLength;
		this.maxNumClusters = maxNumClusters;
		table=new int [kmerLength*maxNumClusters][DNASequence.BASES_ARRAY.length];
	}
	 
	
	public int getNumClusters() {
		return this.numClusters;
	}
	
	/**
	 * Adds the given k-mer to the cluster
	 * @param kmer New kmer
	 */
	public void addOcurrance(CharSequence seq) {
		if (!(seq instanceof DNAShortKmer)) throw new IllegalArgumentException("This class only can process objects of the DNAShortKmer class");
		DNAShortKmer kmer = (DNAShortKmer) seq;
		Integer k = inexactSearchKmerCluster(kmer);
		if(k != null) {
			append(kmer, k);
		} else {
			if (numClusters < maxNumClusters) {
				createCluster(kmer);
			}
		}
	}
	@Override
	public void setCount(CharSequence kmer, int count) {
		if(count>Integer.MAX_VALUE) count = Integer.MAX_VALUE;
		index.put(new DNAShortKmer(kmer), count);
	}
	/**
	 * Searches the hashmap for a matching kmer. If it is not
	 * found, it looks for a kmer that is one nucleotide apart
	 * (e.g ACATCCC[...] would match with ACGTCCC[...]).
	 * @param kmer
	 * @return int k or null if no kmer or neighboring kmer found.
	 */
	private Integer inexactSearchKmerCluster (DNAShortKmer kmer) {
		Integer k = index.get(kmer);
		if(k != null) return k;
		String kmerStr = kmer.toString();
		for(int i = 0; i < kmerStr.length(); i++) {
			int bpIdx = DNASequence.BASES_STRING.indexOf(kmerStr.charAt(i));
			for(int j = 0; j < DNASequence.BASES_STRING.length(); j++) {
				if(j==bpIdx) continue;
				DNAShortKmer neighbor = new DNAShortKmer(kmerStr.subSequence(0, i) + "" + DNASequence.BASES_STRING.charAt(j) + kmerStr.subSequence(i+1, kmerStr.length())); 
				k = index.get(neighbor);
				if(k != null) return k;
			}
		}
		return null;
	}
	
	/**
	 * Creates a new cluster with index newIndex++, associates it with the given
	 * kmer and adds the kmer to the table.
	 * @param kmer
	 */
	private void createCluster(DNAShortKmer kmer) {
		for(int i = 0; i < kmer.length(); i++) {
			char bp = Character.toUpperCase(kmer.charAt(i));
			if (!DNASequence.isInAlphabeth(bp)) {
				throw new RuntimeException("Non DNA kmer found " + kmer);
			}
			int ibp = DNASequence.BASES_STRING.indexOf(bp);
			table[newIndex*kmerLength + i][ibp]++;
		}
		index.put(kmer, newIndex++);
		numClusters++;
	}
	
	/**
	 * Finds the representative kmer for the given k, and updates the table
	 * by increasing the count of the correct nucleotide for each char of the
	 * kmer. (i.e if kmer is ACAT[...] the count for [0][0], [1][1], [2][0],
	 * [3][3], etc. are updated by one. The first index indicates the position
	 * on the kmer and the second index indicates the corresponding nucleotide based
	 * on BASE_ARRAY.)
	 *
	 * If after the update, the representative kmer has changed, the hashmap is updated
	 * to reflect this.
	 * @param kmer
	 * @param k
	 */
	private void append(DNAShortKmer kmer, int k) {
		DNAShortKmer oldKmer = getRepresentativeKmer(k);
		String oldKmerStr = oldKmer.toString();
		String kmerStr = kmer.toString();
		boolean alternative = !oldKmerStr.equals(kmerStr);
		for(int i = 0; i < kmer.length(); i++) {
			int j = DNASequence.BASES_STRING.indexOf(kmerStr.charAt(i));
			table[k*kmerLength + i][j]++;
		}
		if(!alternative) return;
		DNAShortKmer newKmer = getRepresentativeKmer(k);
		if(oldKmerStr.equals(newKmer.toString())) return;
		index.remove(oldKmer);
		index.put(newKmer, k);
	}
	
	/**
	 * Finds the kmer with the most likely sequence. (i.e. for each
	 * char of the kmer, it looks at the cluster table to find the
	 * nucleotide with most occurrences).
	 * @param k
	 * @return DNAShortKmer consensus
	 */
	private DNAShortKmer getRepresentativeKmer(int k) {
		char[] consensus = new char[kmerLength];
		
		for(int i = 0; i < consensus.length; i++) {
			int max = 0;
			for(int j = 0; j < DNASequence.BASES_STRING.length(); j++) {
				int next = table[k*kmerLength + i][j];
				if(max <= next) {
					consensus[i] = DNASequence.BASES_STRING.charAt(j);
					max = next;
				}
			}
		}
		return new DNAShortKmer(new String(consensus));
	}
	@Override
	public int size() {
		return index.size();
	}
	@Override
	public int getCount(CharSequence seq) {
		if (!(seq instanceof DNAShortKmer)) throw new IllegalArgumentException("This class only can process objects of the DNAShortKmer class");
		DNAShortKmer kmer = (DNAShortKmer) seq;
		Integer k = inexactSearchKmerCluster(kmer);
		if(k != null) {
			return getCount(k);
		}
		return 0;
	}
	
	/**
	 * Counts the number of kmers in a given cluster by adding the number of
	 * occurrences of each nucleotide
	 * @param k
	 * @return int count
	 */
	private int getCount(Integer k) {
		int count = 0;
		for(int j = 0; j < table[k].length; j++) {
			count += table[k*kmerLength][j];
		}
		return count;
	}
	@Override
	public void filterKmers(int minAbundance) {
		// TODO: Implement
		
	}
	@Override
	public Distribution calculateAbundancesDistribution() {
		Distribution kmerSpectrum = new Distribution(1, 200, 1);
		Iterator<Entry<DNAShortKmer, Integer>> it = index.entrySet().iterator();
		while (it.hasNext()) {
			Entry<DNAShortKmer, Integer> entry = it.next();
			int k = entry.getValue();
		    kmerSpectrum.processDatapoint(getCount(k));
		}
		return kmerSpectrum;
	}
	public Integer getCluster(DNAShortKmer kmer) {
		return inexactSearchKmerCluster(kmer);
	}
	/**
	 * Disposes memory resources associated with this table
	 */
	public void dispose () {
		table = null;
	}


	@Override
	public void save(PrintStream out) {
		Iterator<Entry<DNAShortKmer, Integer>> it = index.entrySet().iterator();
		while (it.hasNext()) {
			Entry<DNAShortKmer, Integer> entry = it.next();
			out.println(entry.getKey().toString()+"\t"+entry.getValue());
		}
	}
}
