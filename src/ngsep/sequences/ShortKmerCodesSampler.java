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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * @author Jorge Duitama
 */
public class ShortKmerCodesSampler {
	public static final int DEF_KMER_LENGTH = 25;
	public static final int DEF_WINDOW_LENGTH = 40;
	private ShortKmerCodesSamplingAlgorithm algorithm;
	private ShortKmerCodesHashFunction hashFunction;
	private int kmerLength = DEF_KMER_LENGTH;
	private int windowLength = DEF_WINDOW_LENGTH;
	private LimitedSequence alphabetSequence = DNASequence.EMPTY_DNA_SEQUENCE;
	
	
	
	public int getKmerLength() {
		return kmerLength;
	}
	public void setKmerLength(int kmerLength) {
		this.kmerLength = kmerLength;
	}
	public int getWindowLength() {
		return windowLength;
	}
	public void setWindowLength(int windowLength) {
		this.windowLength = windowLength;
	}
	public ShortKmerCodesSampler () {
		this.algorithm = new MinimapShortKmerCodesSamplingAlgorithm();
		this.hashFunction = new MinimapShortKmerCodesHashFunction();
	}
	public ShortKmerCodesSampler(ShortKmerCodesSamplingAlgorithm algorithm, ShortKmerCodesHashFunction hashFunction) {
		this.algorithm = algorithm;
		this.hashFunction = hashFunction;
	}

	/**
	 * Calculates the selected codes of the given sequence following the same algorithm used for minimizers but saving the codes instead of the hashes
	 * @param sequenceId Id of the sequence to calculate
	 * @param sequence characters of the sequence to calculate
	 * @return List<MinimizersTableEntry> Codes selected for the given sequence.
	 */
	public List<KmerCodesTableEntry> computeSequenceCodes(int sequenceId, String sequence,int start,int end) {
		//new PrimeNumbers(1000000);
		long [] segmentCodes = KmersExtractor.extractKmerCodes(sequence, kmerLength, start, Math.min(sequence.length(), end+windowLength+kmerLength),alphabetSequence,false);
		List<KmerCodesTableEntry> selectedCodes = computeSequenceCodes(sequenceId, start, segmentCodes);
		//System.out.println("Selected codes for sequence "+sequenceId+" from "+start+" to "+end+" Number of codes: "+selectedCodes.size()+" pct: "+(100*selectedCodes.size()/segmentCodes.length));
		return selectedCodes;
	}
	public List<KmerCodesTableEntry> computeSequenceCodes(int sequenceId, int start, long[] segmentCodes) {
		int debugIdx = -2;
		List<KmerCodesTableEntry> answer = new ArrayList<KmerCodesTableEntry>();
		Integer [] hashcodes = new Integer [segmentCodes.length];
		Arrays.fill(hashcodes, null);
		for(int i=0;i<segmentCodes.length;i++) {
			long code = segmentCodes[i];
			if(code<0) continue;
			hashcodes[i]=hashFunction.getHash(code);
		}
		if(sequenceId==debugIdx) System.err.println("Calculated "+segmentCodes.length+" hash codes");
		boolean [] selected = algorithm.sample(hashcodes);
		for(int i=0;i<hashcodes.length;i++) {
			if(selected[i]) {
				long originalCode = segmentCodes[i];
				int globalStart = start+i;
				KmerCodesTableEntry entry = new KmerCodesTableEntry(originalCode, sequenceId, globalStart);
				answer.add(entry);
				//if(globalStart>56000 && globalStart<58000) System.err.println("New minimizer calculated Start: "+globalStart+" new code: "+originalCode+" kmer: "+new String (DNASequence.getDNASequence(originalCode, kmerLength)));
			}
		}
		if(sequenceId==debugIdx) System.err.println("Selected codes for sequence "+sequenceId+" from "+start+" Codes: "+answer.size());
		return answer;
	}
	/**
	 * Calculates the selected codes of the given sequence following the same algorithm used for minimizers but saving the codes instead of the hashes
	 * @param sequenceId Id of the sequence to calculate
	 * @param sequence characters of the sequence to calculate
	 * @return List<MinimizersTableEntry> Codes selected for the given sequence.
	 */
	public Map<Integer,Long> computeSequenceCodesAsMap(String sequence,int start,int end) {
		List<KmerCodesTableEntry> entries = computeSequenceCodes(-1, sequence, start , end);
		Map<Integer,Long> answer = new LinkedHashMap<>();
		for(KmerCodesTableEntry entry:entries) {
			answer.put(entry.getStart(), entry.getKmerCode());
		}
		return answer;
	}
	public Map<Integer,Long> computeSequenceCodesAsMap(String sequence,int start,int end, long [] kmerCodes) {
		List<KmerCodesTableEntry> entries = computeSequenceCodes(-1, start, kmerCodes);
		Map<Integer,Long> answer = new LinkedHashMap<>();
		for(KmerCodesTableEntry entry:entries) {
			answer.put(entry.getStart(), entry.getKmerCode());
		}
		return answer;
	}
}
