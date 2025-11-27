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
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Queue;

import ngsep.math.Distribution;

/**
 * @author Jorge Gomez
 * @author Jorge Duitama
 */
public class DNAShortKmerClusterMap implements KmersMap {
	private int kmerLength;
	private int maxNumClusters;
	private short[][] table;
	private DNAShortKmer [] consensusByClusterId;
	private Map<DNAShortKmer, Integer> index = new HashMap<>();
	private int newIndex = 0;
	private Queue<Integer> clusterIdsToReuse = new LinkedList<Integer>();
	private Queue<Integer> clusterIdsToEvaluate = new LinkedList<Integer>();
	
	public DNAShortKmerClusterMap (int kmerLength, int maxNumClusters) {
		this.kmerLength = kmerLength;
		this.maxNumClusters = maxNumClusters;
		table=new short [kmerLength*maxNumClusters][DNASequence.BASES_ARRAY.length];
		consensusByClusterId = new DNAShortKmer[maxNumClusters];
		Arrays.fill(consensusByClusterId, null);
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
			createCluster(kmer);
		}
	}
	@Override
	public void setCount(CharSequence kmer, int count) {
		//TODO: Make a good implementation
		//index.put(new DNAShortKmer(kmer), count);
		throw new RuntimeException("Method not implemented");
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

	public void eliminateShallowClusters() {
		System.out.println("\tCurrent total clusters: "+size()+" Clusters to evaluate: " + clusterIdsToEvaluate.size());
		int removed = 0;
		while(clusterIdsToEvaluate.size()>0) {
			int nextToEvaluate = clusterIdsToEvaluate.remove();
			int depth = getCount(nextToEvaluate);
			if(depth <= 3) {
				clusterIdsToReuse.add(nextToEvaluate);
				removeCluster(nextToEvaluate);
				removed++;
			}
		}
		System.out.println("\tremoved: " + removed + " clusters. Remaining: "+size());
		System.out.println("\tCluster ids to reuse: " + clusterIdsToReuse.size());
	}
	
	/**
	 * Creates a new cluster with index newIndex++, associates it with the given
	 * kmer and adds the kmer to the table.
	 * @param kmer
	 */
	private boolean createCluster(DNAShortKmer kmer) {
		int clusterIndex;
		// Checks if there are any empty cluster positions to occupy
		if(clusterIdsToReuse.isEmpty()) {
			clusterIndex = newIndex;
			if(clusterIndex>=maxNumClusters) return false;
			newIndex++;
		} else {
			clusterIndex = clusterIdsToReuse.remove();
		}
		for(int i = 0; i < kmer.length(); i++) {
			char bp = Character.toUpperCase(kmer.charAt(i));
			if (!DNASequence.EMPTY_DNA_SEQUENCE.isInAlphabet(bp)) {
				throw new RuntimeException("Non DNA kmer found " + kmer);
			}
			int ibp = DNASequence.BASES_STRING.indexOf(bp);
			table[clusterIndex*kmerLength + i][ibp]=1;
		}
		clusterIdsToEvaluate.add(clusterIndex);
		consensusByClusterId[clusterIndex] = kmer;
		index.put(kmer, clusterIndex);
		return true;
	}
	
	private void checkClusterMem(int i, int k) {
		short count = 0;
		int row = k*kmerLength + i;
		for(int j=0;j < DNASequence.BASES_STRING.length(); j++) {
			if(table[row][j] == Short.MAX_VALUE) {
				count++;
			}
		}
		if(count >= 2){
			System.err.print("WARNING: counts for cluster "+k+" position "+i+" in K-mer table has surpassed the length of " +
					Short.MAX_VALUE+". The counts of the other possible nucleotides are");
			for(int m = 0; m < DNASequence.BASES_STRING.length(); m++) {
				System.err.print(" "+DNASequence.BASES_STRING.charAt(m)+": "+ table[row][m]);
			}
			System.err.println();
		}
	}
	
	private void removeCluster (int k) {
		DNAShortKmer consensus = consensusByClusterId[k];
		consensusByClusterId[k] = null;
		index.remove(consensus);
		for(int i = 0; i < consensus.length(); i++) {
			int row = k*kmerLength + i;
			Arrays.fill(table[row],(short)0); 
		}
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
		DNAShortKmer oldKmer = consensusByClusterId[k];
		String oldKmerStr = oldKmer.toString();
		String kmerStr = kmer.toString();
		boolean alternative = !oldKmerStr.equals(kmerStr);
		for(int i = 0; i < kmer.length(); i++) {
			int row = k*kmerLength + i;
			int j = DNASequence.BASES_STRING.indexOf(kmerStr.charAt(i));
			if(table[row][j] < Short.MAX_VALUE) {
				table[row][j]++;
			} else {
				checkClusterMem(i,k);
			}
		}
		if(!alternative) return;
		DNAShortKmer newKmer = calculateRepresentativeKmer(k);
		if(oldKmerStr.equals(newKmer.toString())) return;
		consensusByClusterId[k] = newKmer;
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
	private DNAShortKmer calculateRepresentativeKmer(int k) {
		char[] consensus = new char[kmerLength];
		
		for(int i = 0; i < consensus.length; i++) {
			int max = 0;
			for(int j = 0; j < DNASequence.BASES_STRING.length(); j++) {
				short next = table[k*kmerLength + i][j];
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


	@Override
	public List<CharSequence> getKmersWithCount(int count) {
		// TODO Auto-generated method stub
		return null;
	}
}
