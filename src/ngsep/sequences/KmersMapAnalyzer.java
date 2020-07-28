package ngsep.sequences;

import java.util.Arrays;
import java.util.List;

import ngsep.math.Distribution;

public class KmersMapAnalyzer {
	private KmersMap kmersMap;
	private boolean assembly;
	private int mode;
	private int localMinimum;
	private double average;
	private long expectedAssemblyLength;
	private int numKmersInUniqueZone;
	
	public KmersMapAnalyzer (KmersMap map, boolean assembly) {
		this.kmersMap = map;
		this.assembly = assembly;
		analyzeDistribution(kmersMap.calculateAbundancesDistribution());
	}
	private void analyzeDistribution(Distribution distribution) {
		long totalEntries = 0;
		int maxValueDist = (int)Math.round(distribution.getMaxValueDistribution());
		long [] completeCounts = new long [maxValueDist+1];
		for(int i=1;i<=maxValueDist;i++) {
			completeCounts[i] = Math.round(i*distribution.getDistributionCount(i));
			totalEntries+=completeCounts[i];
		}
		localMinimum=1;
		if(assembly) {
			mode = 1;
			expectedAssemblyLength = totalEntries/2;
			average = distribution.getAverage();
		} else {
			for(int k=10;k<=20;k++) {
				long maxValue = 0;
				for(int i=k;i<=maxValueDist;i++) {
					long value = completeCounts[i];
					if(value>maxValue) {
						maxValue = value;
						mode = i;
					}
				}
				if(mode>k) break;
			}
			long localMinValue = completeCounts[1];
			localMinimum = 1;
			for(int i=2;i<mode;i++) {
				long value = completeCounts[i];
				if(value<localMinValue) {
					localMinValue = value;
					localMinimum = i;
				}
			}
			localMinimum = Math.max(localMinimum, mode/2);
			long errorSum = 0;
			for(int i=1;i<localMinimum;i++) {
				errorSum+=completeCounts[i];
			}
			int diff = mode-localMinimum;
			expectedAssemblyLength = (totalEntries-errorSum) / (2*mode);
			average = 0;
			numKmersInUniqueZone = 0;
			double count = 0;
			for(int i=localMinimum;i<maxValueDist;i++) {
				long iCount = Math.round(distribution.getDistributionCount(i));
				if(numKmersInUniqueZone+iCount<Integer.MAX_VALUE && i<=mode+diff) numKmersInUniqueZone += iCount;
				average+=completeCounts[i];
				count+=iCount;
			}
			average/=count;
		}
		
		int maxDepthPrint = 5*Math.max(mode, (int)average);
		
		if(assembly) {
			distribution.printDistribution(System.out,true, maxDepthPrint);
		} else {
			for(int i=1;i<=maxDepthPrint;i++) {
				System.out.println(""+i+"\t"+distribution.getDistributionCount(i)+"\t"+completeCounts[i]);
			}
			long remainderCounts = 0;
			long remainderComplete = 0;
			for(int i=maxDepthPrint+1;i<=maxValueDist;i++) {
				remainderCounts+=distribution.getDistributionCount(i);
				remainderComplete+=completeCounts[i];
			}
			System.out.println("More\t"+remainderCounts+"\t"+remainderComplete);
			System.out.println("First minimum: "+localMinimum+" total value: "+completeCounts[localMinimum]);
			System.out.println("Local mode: "+mode+" total value: "+completeCounts[mode]);
			System.out.println("Average removing segment until first local minimum: "+average);
			System.out.println("Estimated assembly length: "+expectedAssemblyLength);
		}
	}
	public KmersMap getKmersMap() {
		return kmersMap;
	}
	public boolean isAssembly() {
		return assembly;
	}
	public int getMode() {
		return mode;
	}
	public double getAverage() {
		return average;
	}
	public long getExpectedAssemblyLength() {
		return expectedAssemblyLength;
	}
	public long[] extractKmerCodesInLocalSDZone() {
		long [] answer = new long [numKmersInUniqueZone];
		Arrays.fill(answer, -1);
		int localSD = getModeLocalSD();
		int idxAnswer = 0;
		for(int i=mode;i>=mode-localSD;i--) {
			idxAnswer = addKmerCodes(i,answer,idxAnswer);
			if(idxAnswer==answer.length) return answer;
			
		}
		for(int i=mode;i<=mode+localSD;i++) {
			idxAnswer = addKmerCodes(i,answer,idxAnswer);
			if(idxAnswer==answer.length) return answer;
		}
		return answer;
	}
	private int addKmerCodes(int depth, long[] answer, int idxAnswer) {
		if(kmersMap instanceof ShortArrayDNAKmersMapImpl) {
			List<Integer> codes = ((ShortArrayDNAKmersMapImpl)kmersMap).getKmerCodesWithCount(depth);
			for(int code:codes) {
				answer[idxAnswer] = code;
				idxAnswer++;
				if(idxAnswer==answer.length) return idxAnswer;
			}
		} else {
			List<CharSequence> kmersDepth = kmersMap.getKmersWithCount(depth);
			for(CharSequence kmer:kmersDepth) {
				long code = AbstractLimitedSequence.getHash(kmer, 0, kmer.length(), DNASequence.EMPTY_DNA_SEQUENCE);
				answer[idxAnswer] = code;
				idxAnswer++;
				if(idxAnswer==answer.length) return idxAnswer;
			}
		}
		return idxAnswer;
	}
	public int getModeLocalSD() {
		return mode-localMinimum;
	}
	
}
