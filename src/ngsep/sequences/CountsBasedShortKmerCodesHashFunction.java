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

import ngsep.math.PrimeNumbers;

/*
 * @author Jorge Duitama
 */
public class CountsBasedShortKmerCodesHashFunction implements ShortKmerCodesHashFunction {
	private int kmerLength;
	private KmersMap kmersMap;
	private KmersMapAnalyzer kmersAnalyzer;
	private PrimeNumbers primeNumbersHelper;
	
	
	public CountsBasedShortKmerCodesHashFunction(int kmerLength, KmersMapAnalyzer kmersAnalyzer) {
		super();
		this.kmerLength = kmerLength;
		this.kmersAnalyzer = kmersAnalyzer;
		this.kmersMap = kmersAnalyzer.getKmersMap();
		int mode = Math.max(1,kmersAnalyzer.getMode());
		//TODO: This is ok for assemblies but it is not ok for reference genomes
		int capacityPrimes = (int) Math.min(100000000, kmersAnalyzer.getNumKmers(mode));
		this.primeNumbersHelper = new PrimeNumbers(capacityPrimes);
	}


	@Override
	public int getHash(long dnaHash) {
		int prime = 1073676287;
		int count;
		if(kmersMap instanceof ShortArrayDNAKmersMapImpl) {
			//Select the 15-mer suffix
			count = ((ShortArrayDNAKmersMapImpl)kmersMap).getCount(dnaHash & 0x3FFFFFFF);
			//count = ((ShortArrayDNAKmersMapImpl)kmersMap).getCount(dnaHash);
			//count = 20;
		} else {
			String kmer = new String(AbstractLimitedSequence.getSequence(dnaHash, kmerLength, DNASequence.EMPTY_DNA_SEQUENCE));
			count = kmersMap.getCount(kmer);
		}
		long hash = Integer.MAX_VALUE;
		if(count > 0) {
			long rankingStart=kmersAnalyzer.getRanking(count);
			long kmersWithCount = kmersAnalyzer.getNumKmers(count);
			if(kmersWithCount<primeNumbersHelper.getCapacity()) prime = primeNumbersHelper.getNextPrime((int) kmersWithCount);
			hash = rankingStart+(dnaHash%prime);
			if(hash>=Integer.MAX_VALUE) hash = Integer.MAX_VALUE-1;
		}
		return (int)hash;
	}

}
