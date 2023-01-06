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
package ngsep.benchmark;

import java.io.IOException;
import java.io.PrintStream;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.KmersExtractor;
import ngsep.sequences.KmersMap;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.ShortArrayDNAKmersMapImpl;
import ngsep.sequences.ShortKmerCodesTable;
import ngsep.sequences.io.FastaFileReader;
import ngsep.sequences.io.KmersMapLoader;

/**
 * 
 * @author Jorge Duitama
 *
 */
public class KmerBasedSwitchErrorsFinder {

	private KmersMap kmersHap1;
	private KmersMap kmersHap2;
	private KmersMap kmersReads;
	private Set<Long> uniqueKmerCodesHap1;
	private Set<Long> uniqueKmerCodesHap2;
	
	private KmersMapLoader loader = new KmersMapLoader();
	private ShortKmerCodesTable table = new ShortKmerCodesTable(15, 40);
	
	public static void main(String[] args) throws IOException {
		KmerBasedSwitchErrorsFinder instance = new KmerBasedSwitchErrorsFinder();
		String kmersHap1File = args[0];
		String kmersHap2File = args[1];
		//String kmersReadsFile = args[2];
		String assemblyFile = args[2];
		instance.loadKmers(kmersHap1File,1);
		instance.loadKmers(kmersHap2File,2);
		//instance.loadKmers(kmersReadsFile,3);
		instance.selectPhaseInformativeKmers();
		instance.processAssembly(assemblyFile,System.out);
	}

	

	private void loadKmers(String kmersFile, int dataset) throws IOException {
		KmersMap kmers = loader.loadKmersMap(kmersFile, 15);
		if (dataset==1) kmersHap1 = kmers;
		if (dataset==2) kmersHap2 = kmers;
		if (dataset==3) kmersReads = kmers;
	}
	
	private void selectPhaseInformativeKmers() {
		//TODO: make it work for any KmersMap
		uniqueKmerCodesHap1 = new HashSet<>();
		uniqueKmerCodesHap2 = new HashSet<>();
		for(long code = 0; code<Integer.MAX_VALUE;code++) {
			int countH1 = ((ShortArrayDNAKmersMapImpl)kmersHap1).getCount(code);
			int countH2 = ((ShortArrayDNAKmersMapImpl)kmersHap2).getCount(code);
			int countReads = kmersReads!=null?((ShortArrayDNAKmersMapImpl)kmersHap1).getCount(code):-1;
			if(countH1==1 && countH2==0 && (countReads==-1 || countReads < 100)) uniqueKmerCodesHap1.add(code);
			if(countH1==0 && countH2==1 && (countReads==-1 || countReads < 100)) uniqueKmerCodesHap2.add(code);
		}
		System.out.println("Phase informative kmers: "+uniqueKmerCodesHap1.size()+" "+uniqueKmerCodesHap2.size());
	}
	private int calculateCountHaplotypes(Map<Integer, Long> codes) {
		int count = 0;
		for(long code:codes.values()) {
			//TODO: make it work for any KmersMap
			int countH1 = ((ShortArrayDNAKmersMapImpl)kmersHap1).getCount(code);
			int countH2 = ((ShortArrayDNAKmersMapImpl)kmersHap2).getCount(code);
			if(countH1 ==1) count++;
			if(countH2 ==1) count++;
		}
		return count;
	}

	private void processAssembly(String assemblyFile, PrintStream out) throws IOException {
		int switchErrors = 0;
		try (FastaFileReader reader = new FastaFileReader(assemblyFile)) {
			Iterator <QualifiedSequence> it = reader.iterator();
			while(it.hasNext()) {
				QualifiedSequence qseq = it.next();
				DNAMaskedSequence seq = (DNAMaskedSequence) qseq.getCharacters();
				String seqStr = seq.toString();
				Map<Integer, Long> kmerCodesF = KmersExtractor.extractDNAKmerCodes(seqStr, 15, 0, seqStr.length());	
				Map<Integer,Long> minimizersF = table.computeSequenceCodesAsMap(seqStr, 0, seqStr.length(),kmerCodesF);
				DNAMaskedSequence rcseq = seq.getReverseComplement();
				String rcseqStr = rcseq.toString();
				Map<Integer, Long> kmerCodesR = KmersExtractor.extractDNAKmerCodes(rcseqStr, 15, 0, rcseqStr.length());
				Map<Integer,Long> minimizersR = table.computeSequenceCodesAsMap(rcseqStr, 0, rcseq.length(),kmerCodesR);
				int strand = inferStrand (qseq.getName(),minimizersF,minimizersR);
				if(strand == 1) switchErrors+= processKmerCodes(qseq.getName(), seqStr.length(), kmerCodesF, out);
				else if(strand == 2) switchErrors+= processKmerCodes(qseq.getName(), rcseqStr.length(), kmerCodesR, out);
				
			}
		}
		out.println("Total switch errors: "+switchErrors);
	}

	private int inferStrand(String sequenceId, Map<Integer, Long> minimizersF, Map<Integer, Long> minimizersR) {
		int countF = calculateCountHaplotypes(minimizersF);
		int countR = calculateCountHaplotypes(minimizersR);
		if(countF>2*countR && countF > 100) return 1;
		if(countR>2*countF && countR > 100) return 2;
		System.out.println("Strand undefined for sequence: "+sequenceId+" counts: "+countF+" "+countR);
		return 0;
	}

	
	private int processKmerCodes(String seqId, int seqLen, Map<Integer, Long> kmerCodes, PrintStream out) {
		int nw = seqLen / 10000 + 1;
		int [][] windowCounts = new int [nw][2];  
		for(Map.Entry<Integer, Long> entry:kmerCodes.entrySet()) {
			int start = entry.getKey();
			long code = entry.getValue();
			int w = start/20000;
			if(uniqueKmerCodesHap1.contains(code)) {
				windowCounts[w][0]++;
				if(start%20000>=10000 && w<windowCounts.length-1) windowCounts[w+1][0]++;
				else if (start%20000<10000 && w>0) windowCounts[w-1][0]++;
			} else if (uniqueKmerCodesHap2.contains(code)) {
				windowCounts[w][1]++;
				if(start%20000>=10000 && w<windowCounts.length-1) windowCounts[w+1][1]++;
				else if (start%20000<10000 && w>0) windowCounts[w-1][1]++;
			}
		}
		int switchErrors = 0;
		int informativeWindows = 0;
		int countHap1=0;
		int countHap2=0;
		int countW1=0;
		int countW2=0;
		double lastProp=0;
		int haplotype = 0;
		for(int i=0;i<windowCounts.length;i++) {
			int sum = windowCounts[i][0] +windowCounts[i][1];
			if(sum<10) continue;
			double prop = 1.0*windowCounts[i][0] / sum;
			int nextHap = 0;
			countHap1+=windowCounts[i][0];
			countHap2+=windowCounts[i][1];
			
			if(prop>0.8) {
				nextHap = 1;
				countW1++;
			}
			if(prop<0.2) {
				nextHap = 2;
				countW2++;
			}
			if(nextHap>0) {
				informativeWindows++;
				if(haplotype == 0) haplotype = nextHap;
				else if (haplotype != nextHap) {
					System.out.println("Switch error at sequence "+seqId+" window "+i+" props "+lastProp +" "+prop+ " counts: "+countHap1+" "+countHap2+" current haplotype: "+haplotype+ " next "+nextHap);
					switchErrors++;
					haplotype = nextHap;
				}
			}
			lastProp = prop;
		}
		System.out.println("Total informative windows "+informativeWindows+" total counts: "+countHap1+" "+countHap2+" windows assigned: "+countW1+" "+countW2);
		return switchErrors;
	}
}
