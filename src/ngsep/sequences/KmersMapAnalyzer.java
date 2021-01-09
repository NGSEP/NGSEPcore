package ngsep.sequences;

import java.util.Arrays;
import java.util.List;
import java.util.Map;

import ngsep.math.Distribution;

public class KmersMapAnalyzer {
	private KmersMap kmersMap;
	private long [] kmerCounts;
	private long [] completeCounts;
	private long [] countRankings;
	private boolean assembly;
	private int mode;
	private int localMinimum;
	private double average;
	private long expectedAssemblyLength;
	
	public KmersMapAnalyzer (KmersMap map, boolean assembly) {
		this.kmersMap = map;
		this.assembly = assembly;
		analyzeDistribution(kmersMap.calculateAbundancesDistribution());
	}
	private void analyzeDistribution(Distribution distribution) {
		long totalEntries = 0;
		int maxValueDist = (int)Math.round(distribution.getMaxValueDistribution());
		completeCounts = new long [maxValueDist+1];
		kmerCounts = new long [maxValueDist+1];
		countRankings = new long [maxValueDist+1];
		for(int i=1;i<=maxValueDist;i++) {
			kmerCounts[i] = Math.round(distribution.getDistributionCount(i));
			completeCounts[i] = Math.round(i*distribution.getDistributionCount(i));
			totalEntries+=completeCounts[i];
		}
		localMinimum=1;
		if(assembly) {
			mode = 1;
			expectedAssemblyLength = totalEntries/2;
			average = distribution.getAverage();
		} else {
			mode = Math.max(10, 3*(int)distribution.getAverage());
			int limit = Math.min(10000, 100*(int)distribution.getAverage());
			for(int k=10;k<limit;k++) {
				long maxValue = completeCounts[k];
				int localMode = k;
				for(int i=k+1;i<=limit;i++) {
					long value = completeCounts[i];
					if(value>maxValue) {
						maxValue = value;
						localMode = i;
					}
				}
				if(localMode>k) {
					mode = localMode;
					break;
				}
			}	
			countRankings[mode]=0;
			long sum=kmerCounts[mode];
			for(int k=1;k<=maxValueDist-mode;k++) {
				int idx = mode+k;
				long count = kmerCounts[idx];
				countRankings[idx]=sum;
				sum+=count;
				idx = mode-k;
				if(idx>=5) {
					count = kmerCounts[idx];
					countRankings[idx]=sum;
					sum+=count;
				}
			}
			//Worst rankings for low abundance kmers
			for(int i=4;i>=1;i--) {
				long count = kmerCounts[i];
				countRankings[i]=sum;
				sum+=count;
			}
			long localMinValue = completeCounts[1];
			for(int i=2;i<mode;i++) {
				long value = completeCounts[i];
				if(value<localMinValue) {
					localMinValue = value;
					localMinimum = i;
				}
			}
			if(localMinimum==1 || localMinimum==mode-1 || localMinValue>completeCounts[1]/2 || localMinValue>completeCounts[mode]/2) {
				//Issues calculating local minimum. Set to half of the mode
				localMinimum = Math.max(5, mode/2);
			}
			
			long errorSum = 0;
			for(int i=1;i<localMinimum;i++) {
				errorSum+=completeCounts[i];
			}
			expectedAssemblyLength = (totalEntries-errorSum) / (2*mode);
			average = 0;
			double count = 0;
			for(int i=localMinimum;i<maxValueDist;i++) {
				average+=completeCounts[i];
				count+=kmerCounts[i];
			}
			average/=count;
		}
		
		int maxDepthPrint = 5*Math.max(mode, (int)average);
		
		if(assembly) {
			distribution.printDistribution(System.out,true, maxDepthPrint);
		} else {
			for(int i=1;i<=maxDepthPrint;i++) {
				System.out.println(""+i+"\t"+kmerCounts[i]+"\t"+completeCounts[i]+"\t"+countRankings[i]);
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
		int minValueKmers;
		int maxValueKmers;
		int numKmersToSort;
		if (assembly) {
			minValueKmers = 1;
			maxValueKmers = 1;
			numKmersToSort = (int)kmerCounts[1];
			
			for(int i=2;i<20;i++) {
				long newCount = numKmersToSort+kmerCounts[i];
				if(newCount>20000000) {
					break;
				}
				maxValueKmers = i;
				numKmersToSort=(int)newCount;
			}
		} else {
			int localSD = getModeLocalSD();
			minValueKmers = localMinimum;
			maxValueKmers = mode+localSD;
			numKmersToSort = (int)kmerCounts[mode];
			for(int i=1;i<=localSD;i++) {
				long newCount = numKmersToSort+kmerCounts[mode+i]+kmerCounts[mode-i];
				if(newCount>20000000 || newCount>2*expectedAssemblyLength) {
					minValueKmers=mode-i;
					maxValueKmers=mode+i;
					break;
				}
				numKmersToSort=(int)newCount;
			}
		}
		
		//System.out.println("Mode: "+mode+" limits: "+minValueKmers+" "+maxValueKmers+ " kmers to sort: "+numKmersToSort);
		long [] answer = new long [numKmersToSort];
		Arrays.fill(answer, -1);
		int idxAnswer = 0;
		if(kmersMap instanceof ShortArrayDNAKmersMapImpl) {
			Map<Integer,Short> codes = ((ShortArrayDNAKmersMapImpl)kmersMap).getKmerCodesWithCount(minValueKmers,maxValueKmers);
			for(int i=0;mode+i<=maxValueKmers;i++) {
				idxAnswer = addKmerCodes(codes,mode+i,answer,idxAnswer);
				if(idxAnswer==answer.length) return answer;
				if(i==0 || mode-i<minValueKmers) continue;
				idxAnswer = addKmerCodes(codes,mode-i,answer,idxAnswer);
				if(idxAnswer==answer.length) return answer;
			}
		} else {
			for(int i=0;mode+i<=maxValueKmers;i++) {
				idxAnswer = addKmerCodes(mode+i,answer,idxAnswer);
				if(idxAnswer==answer.length) return answer;
				if(i==0 || mode-i<minValueKmers) continue;
				idxAnswer = addKmerCodes(mode-i,answer,idxAnswer);
				if(idxAnswer==answer.length) return answer;
			}
		}
		
		return answer;
	}
	private int addKmerCodes(Map<Integer, Short> codes, int depth, long[] answer, int idxAnswer) {
		for(int code:codes.keySet()) {
			int value = codes.get(code);
			if(value == depth) {
				answer[idxAnswer]=code;
				idxAnswer++;
				if(idxAnswer==answer.length) return idxAnswer;
			}
		}
		return idxAnswer;
	}
	private int addKmerCodes(int depth, long[] answer, int idxAnswer) {
		List<CharSequence> kmersDepth = kmersMap.getKmersWithCount(depth);
		for(CharSequence kmer:kmersDepth) {
			long code = AbstractLimitedSequence.getHash(kmer, 0, kmer.length(), DNASequence.EMPTY_DNA_SEQUENCE);
			answer[idxAnswer] = code;
			idxAnswer++;
			if(idxAnswer==answer.length) return idxAnswer;
		}
		return idxAnswer;
	}
	public int getModeLocalSD() {
		return Math.max(1, mode-localMinimum);
	}
	public long getRanking(int count) {
		return countRankings[count];
	}
	public long getNumKmers(int count) {
		return kmerCounts[count];
	}
	
}
