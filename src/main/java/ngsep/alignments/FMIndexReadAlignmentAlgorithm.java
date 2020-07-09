package ngsep.alignments;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import ngsep.genome.GenomicRegion;
import ngsep.genome.GenomicRegionImpl;
import ngsep.genome.ReferenceGenomeFMIndex;
import ngsep.genome.io.SimpleGenomicRegionFileHandler;
import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.KmerHitsCluster;
import ngsep.sequences.KmersExtractor;
import ngsep.sequences.LimitedSequence;
import ngsep.sequences.PairwiseAlignmentAffineGap;
import ngsep.sequences.RawRead;
import ngsep.sequences.UngappedSearchHit;

public class FMIndexReadAlignmentAlgorithm implements ReadAlignmentAlgorithm {
	
	private int kmerLength;
	private int maxAlnsPerRead;
	private ReferenceGenomeFMIndex fMIndex;
	
	private Map<String, List<GenomicRegion>> knownSTRs;
	
	private Set<String> repetitiveKmers = new HashSet<String>();
	
	private boolean runFullAlignment = true;
	private boolean onlyPositiveStrand = false;
	
	private PairwiseAlignmentAffineGap alignerFullRead = new PairwiseAlignmentAffineGap(1000);
	private PairwiseAlignmentAffineGap alignerSTRsLeft = new PairwiseAlignmentAffineGap(500);
	private PairwiseAlignmentAffineGap alignerSTRsRight = new PairwiseAlignmentAffineGap(500);
	
	// Statistics
	private int fewMismatchesAlns = 0;
	private int completeAlns = 0;
	
	public FMIndexReadAlignmentAlgorithm(ReferenceGenomeFMIndex fMIndex, int kmerLength, int maxAlnsPerRead) {
		this.fMIndex = fMIndex;
		this.kmerLength = kmerLength;
		this.maxAlnsPerRead = maxAlnsPerRead;
		alignerSTRsLeft.setForceEnd1(false);
		alignerSTRsRight.setForceStart1(false);
		alignerFullRead.setForceStart2(false);
		alignerFullRead.setForceEnd2(false);
	}
	public ReferenceGenomeFMIndex getFMIndex() {
		return fMIndex;
	}
	public void setFMIndex(ReferenceGenomeFMIndex fMIndex) {
		this.fMIndex = fMIndex;
	}

	public Map<String, List<GenomicRegion>> getKnownSTRs() {
		return knownSTRs;
	}
	public void setKnownSTRs(Map<String, List<GenomicRegion>> knownSTRs) {
		this.knownSTRs = knownSTRs;
	}
	public int getCompleteAlns() {
		return completeAlns;
	}
	
	public int getFewMismatchesAlns() {
		return fewMismatchesAlns;
	}
	
	@Override
	public List<ReadAlignment> alignRead (RawRead read) {
		List<ReadAlignment> alignments = new ArrayList<>();
		String readSeq = read.getSequenceString();
		String qual = read.getQualityScores();
		String reverseQS = null;
		if(qual == null || qual.length()!=readSeq.length()) {
			qual = RawRead.generateFixedQSString('5', readSeq.length());
			reverseQS = qual;
		} else if (!onlyPositiveStrand) {
			reverseQS = new StringBuilder(qual).reverse().toString();
		}
		String reverseComplement = null;
		if(!onlyPositiveStrand) {
			reverseComplement = DNAMaskedSequence.getReverseComplement(readSeq).toString();
			
		}
		if(readSeq.length()<500) {
			int maxMismatches = 2;
			alignments.addAll(fewMismatchesSingleStrandSearch(readSeq,maxMismatches));
			//System.out.println("Read: "+read.getName()+" Forward exact alignments: "+alignments.size());
			if(reverseComplement!=null) {
				List<ReadAlignment> alnsR = fewMismatchesSingleStrandSearch(reverseComplement,maxMismatches);
				//System.out.println("Read: "+read.getName()+" Reverse exact alignments: "+alnsR.size());
				for (ReadAlignment aln:alnsR) aln.setNegativeStrand(true);
				alignments.addAll(alnsR);
			}
		}
		if(alignments.size()==0) {
			alignments.addAll(alignQueryToReference(readSeq));
			//System.out.println("Read: "+read.getName()+" Forward inexact alignments: "+alignments.size());
			if(reverseComplement!=null) {
				List<ReadAlignment> alnsR = alignQueryToReference(reverseComplement);
				//System.out.println("Read: "+read.getName()+" Reverse inexact alignments: "+alnsR.size());
				for (ReadAlignment aln:alnsR) aln.setNegativeStrand(true);
				alignments.addAll(alnsR);
			}
		} else fewMismatchesAlns++;
		
		//System.out.println("Read: "+read.getName()+" total alignments: "+alignments.size());
		for(ReadAlignment aln:alignments) {
			aln.setReadName(read.getName());
			if(!aln.isNegativeStrand()) aln.setQualityScores(qual);
			else aln.setQualityScores(reverseQS);
		}
		return alignments;
	}
	public List<ReadAlignment> alignQueryToReference (String query) {
		
		return kmerBasedSingleStrandInexactSearchAlgorithm(query);
	}
	
	/**
	 * Inexact search of kmers to an FM-index
	 * It iterates the Sequences of the genome and if there is at least MIN_ACCURACY percentage of the kmers
	 * it allow the alignment with the first an the last position of the kmers ocurrence
	 * Only tries to align the given quey in the positive strand
	 * @return List<ReadAlignment>
	 */
	private List<ReadAlignment> kmerBasedSingleStrandInexactSearchAlgorithm (String query) 
	{
		Map<Integer,CharSequence> kmersMap = KmersExtractor.extractKmersAsMap(query, kmerLength, kmerLength, true, true, true);
		List<ReadAlignment> finalAlignments =  new ArrayList<>();
		//System.out.println("Read: "+query+" length "+query.length()+" kmers: "+kmersMap.size());
		int kmersCount=kmersMap.size();
		if(kmersCount==0) return finalAlignments;
		List<UngappedSearchHit> initialKmerHits = searchKmers (kmersMap);
		List<KmerHitsCluster> clusteredKmerHits = clusterKmerHits(query, initialKmerHits);
		if(clusteredKmerHits.size()==0) return finalAlignments;
		//System.out.println("Initial kmer hits: "+initialKmerHits.size()+" Clusters: "+clusteredKmerHits.size());
		Collections.sort(clusteredKmerHits, (o1, o2) -> o2.getNumDifferentKmers()-o1.getNumDifferentKmers());
		//KmerHitsCluster cluster = clusteredKmerHits.get(0);
		//ReadAlignment readAln = createNewAlignmentFromConsistentKmers(cluster, query);
		//if(readAln!=null) finalAlignments.add(readAln);
		int kmersMaxCluster = 0;
		for (int i=0;i<clusteredKmerHits.size() && i<2*maxAlnsPerRead;i++) {
			KmerHitsCluster cluster = clusteredKmerHits.get(i);
			int numKmers = cluster.getNumDifferentKmers();
			//System.out.println("Processing cluster "+i+" spanning "+cluster.getSequenceName()+":"+cluster.getSubjectPredictedStart()+"-"+cluster.getSubjectPredictedEnd()+" Num kmers: "+cluster.getNumDifferentKmers()+" consistent: "+cluster.isAllConsistent());
			if(i==0) kmersMaxCluster = numKmers;
			else if (finalAlignments.size()>0 && (numKmers<2 || numKmers< 0.5*kmersMaxCluster)) break;
			ReadAlignment readAln = createNewAlignmentFromConsistentKmers(cluster, query);
			if(readAln!=null) finalAlignments.add(readAln);
		}
		//System.out.println("Found "+finalAlignments.size()+" alignments for query: "+query);
		return finalAlignments;
	}


	/**
	 * Searches the given kmers in the fmIndex 
	 * @param kmers to search
	 * @return List of alignments of each kmer. The read number of each alignment contains the kmer number.
	 */
	private List<UngappedSearchHit> searchKmers(Map<Integer,CharSequence> kmersMap) {
		List<UngappedSearchHit> answer = new ArrayList<>();
		for (int start:kmersMap.keySet()) {
			String kmer = kmersMap.get(start).toString();
			if(repetitiveKmers.contains(kmer)) continue;
			List<UngappedSearchHit> kmerHits=fMIndex.exactSearch(kmer);
			//System.out.println("Kmer: "+kmer+" hits: "+kmerHits.size());
			if(kmerHits.size()>50) {
				repetitiveKmers.add(kmer);
				continue;
			}
			for(UngappedSearchHit hit:kmerHits) {
				hit.setQueryIdx(start);
				answer.add(hit);
			}
		}
		return answer;
	}

	private List<KmerHitsCluster> clusterKmerHits(String query, List<UngappedSearchHit> initialKmerHits) {
		List<KmerHitsCluster> clusters = new ArrayList<>();
		Map<String,List<UngappedSearchHit>> hitsBySubjectName = new LinkedHashMap<String, List<UngappedSearchHit>>();
		for(UngappedSearchHit hit:initialKmerHits) {
			List<UngappedSearchHit> hitsSeq = hitsBySubjectName.computeIfAbsent(hit.getSequenceName(), k -> new ArrayList<>());
			hitsSeq.add(hit);
		}
		for(List<UngappedSearchHit> hitsSeq:hitsBySubjectName.values()) {
			Collections.sort(hitsSeq, (hit0,hit1)-> hit0.getStart()-hit1.getStart());
			clusters.addAll(clusterSequenceKmerAlns(query, hitsSeq));
		}
		return clusters;
	}
	
	private List<KmerHitsCluster> clusterSequenceKmerAlns(String query, List<UngappedSearchHit> sequenceHits) {
		List<KmerHitsCluster> answer = new ArrayList<>();
		//System.out.println("Alns to cluster: "+sequenceAlns.size());
		KmerHitsCluster cluster=null;
		for(UngappedSearchHit kmerHit:sequenceHits) {
			if(cluster==null || !cluster.addKmerHit(kmerHit, 0)) {
				cluster = new KmerHitsCluster(query, kmerHit);
				answer.add(cluster);
			}
		}
		return answer;
	}

	private ReadAlignment createNewAlignmentFromConsistentKmers(KmerHitsCluster cluster, String query) {
		String sequenceName = cluster.getSequenceName();
		int first = cluster.getSubjectPredictedStart()+1;
		int last = cluster.getSubjectPredictedEnd();
		int lastPerfect = first+query.length()-1;
		List<Integer> alignment = new ArrayList<Integer>(1);
		alignment.add(ReadAlignment.getAlnValue(query.length(), ReadAlignment.ALIGNMENT_MATCH));
		ReadAlignment aln = buildAln(query, sequenceName, first, lastPerfect, alignment);
		if(aln!=null) {
			//System.out.println("Built alignment at "+sequenceName+":"+aln.getFirst()+"-"+aln.getLast()+" CIGAR: "+aln.getCigarString());
			GenomicRegion region =findTandemRepeat(sequenceName,first,last);
			if(region!=null) {
				ReadAlignment newaln=verifyShortTandemRepeats(aln.getSequenceName(),aln.getFirst(), aln.getLast(),query,region);
				//System.out.println("Found overlapping tandem repeat at "+region.getSequenceName()+":"+region.getFirst()+"-"+region.getLast()+" new aln: "+newaln);
				if(newaln!=null) {
					return newaln;
				}
			}
			if(cluster.getNumDifferentKmers()>2 && cluster.isAllConsistent()) {
				int [] mismatches = countMismatches (query, aln);
				if(mismatches !=null && mismatches[0]<0.05*query.length() && mismatches[1]+mismatches[2] < 0.1*query.length()) {
					int ends = mismatches[1]+mismatches[2]; 
					//if (ends > mismatches[0]) System.err.println("Problem counting mismatches for "+sequenceName+":"+first+" read: "+query+" mismatches: "+mismatches[0]+" "+mismatches[1]+" "+mismatches[2]);
					if (ends>0) aln = buildAln(query, sequenceName, first+mismatches[1], lastPerfect-mismatches[2], encodeAlignment(query.length(),mismatches));
					if(aln!=null) {
						aln.setAlignmentQuality((byte) Math.round(100-5*mismatches[0]));
						aln.setNumMismatches((short) mismatches[0]);
						//System.out.println("Mismatches alignment at "+aln.getSequenceName()+":"+aln.getFirst()+"-"+aln.getLast()+": "+mismatches);
						return aln;
					}
				}
			}
		}
		if(!runFullAlignment) return null;
		//Perform smith waterman
		first = Math.max(1, first-3);
		last = Math.min(fMIndex.getReferenceLength(sequenceName), last+3);
		if(last-first+1>1.5*query.length()) return null;
		CharSequence refSeq = fMIndex.getSequence(sequenceName, first, last);
		if(refSeq == null) return null;
		
		//System.out.println("Aligning reference from "+first+" to "+last+ " to query. length: "+refSeq.length());
		completeAlns++;
		String [] rawAln;
		synchronized (alignerFullRead) {
			rawAln = alignerFullRead.getAlignment(query, refSeq.toString());
		}
		int mismatches = countMismatches(rawAln);
		if(mismatches>0.1*query.length()) return null;
		LinkedList<Integer> alnCodes = ReadAlignment.encodePairwiseAlignment(rawAln);
		aln = buildAln(query, sequenceName, first, last, alnCodes);
		if(aln==null) return null;
		if (!aln.clipBorders(kmerLength)) return null;
		//System.out.println("New genomic coordinates : "+first+"-"+last+" CIGAR:" +aln.getCigarString());
		aln.setAlignmentQuality((byte) Math.round(100-5*mismatches));
		aln.setNumMismatches((short)mismatches);
		
		return aln;
	}
	private List<Integer> encodeAlignment(int length, int[] mismatches) {
		List<Integer> answer = new LinkedList<Integer>();
		int l2 = length-mismatches[1]-mismatches[2];
		if(mismatches[1]>0) {
			answer.add(ReadAlignment.getAlnValue(mismatches[1], ReadAlignment.ALIGNMENT_SKIPFROMREAD));
		}
		answer.add(ReadAlignment.getAlnValue(l2, ReadAlignment.ALIGNMENT_MATCH));
		if(mismatches[2]>0) {
			answer.add(ReadAlignment.getAlnValue(mismatches[2], ReadAlignment.ALIGNMENT_SKIPFROMREAD));
		}
		return answer;
	}
	private int countMismatches(String[] alignedSequences) {
		int answer = 0;
		boolean lastIsGap = true;
		for(int i=0;i<alignedSequences[0].length();i++) {
			char c1 = alignedSequences[0].charAt(i);
			char c2 = alignedSequences[1].charAt(i);
			if(c1==LimitedSequence.GAP_CHARACTER || c2 == LimitedSequence.GAP_CHARACTER) {
				if(!lastIsGap) answer+=2;
				lastIsGap = true;
			} else {
				if(c1!=c2) answer++;
				lastIsGap = false;
			}
		}
		if(lastIsGap) answer-=2;
		return answer;
	}
	private int [] countMismatches(CharSequence query, ReadAlignment aln) {
		int [] answer = {0,0,0};
		CharSequence refS = fMIndex.getSequence(aln.getSequenceName(), aln.getFirst(), aln.getLast());
		if(refS==null) return null;
		String refSeq = refS.toString();
		int lastMismatch = -1;
		boolean startAssigned = false;
		for (int i=0;i<query.length() && i<refSeq.length();i++ ) {
			if(query.charAt(i)!=refSeq.charAt(i)) {
				answer[0]++;
				lastMismatch=i;
			} else if (startAssigned==false && answer[0]+3<i) {
				answer[1]=lastMismatch+1;
				startAssigned = true;
			}
		}
		if (query.length()!=refSeq.length()) {
			answer[0]+=Math.abs(query.length()-refSeq.length());
			answer[2]=Math.max(0, query.length()-refSeq.length());
		} else {
			lastMismatch=refSeq.length();
			int numM =0;
			for (int i=query.length()-1;i>=0;i-- ) {
				if (query.charAt(i)!=refSeq.charAt(i)) {
					lastMismatch = i;
					numM++;
				} else {
					int revIdx = refSeq.length()-1-i;
					if (numM+3<revIdx) {
						answer[2]=refSeq.length()-lastMismatch;
						break;
					}
				}
			}
		}
		
		return answer;
	}
	
	public GenomicRegion findTandemRepeat(String sequenceName, int first, int last) {
		GenomicRegionImpl region = new GenomicRegionImpl(sequenceName, first, last);
		if(knownSTRs==null) return null;
		List<GenomicRegion> l =knownSTRs.get(sequenceName);
		if(l==null) return null;
		return binaryContains(l, 0, l.size()-1, region);
	}

	private GenomicRegion binaryContains(List<GenomicRegion> l, int left, int rigth, GenomicRegion element) {
		if(rigth>=left) {
			int middle = left +(rigth-left)/2;
			GenomicRegion actual = l.get(middle);
			if(isOverlappgingSorted(actual, element)||isOverlappgingSorted(element, actual	)) {
				return actual;
			}
			if(l.get(middle).getFirst()>element.getFirst()) {
				return binaryContains(l, left, middle-1, element);
			}
			return binaryContains(l, middle+1,rigth, element);
		}
		return null;
	}
	
	/**
	 * Creates an alignment taking into account that the region overlaps with a tandem repeat
	 * @param aln
	 * @param read
	 * @param qualityScores
	 * @param region
	 * @return
	 */
	public ReadAlignment verifyShortTandemRepeats(String sequenceName, int first, int last, String read, GenomicRegion region) {
		int firstLeftPart = Math.max(first,1);
		int softClipLeft = 0;
		int softClipRight = 0;
		LinkedList<Integer> encodedLeftAln = null;
		LinkedList<Integer> encodedRightAln = null;
		int leftMismatches = 0;
		int rightMismatches = 0;
		if(first<region.getFirst()-5) {
			CharSequence refSeq = fMIndex.getSequence(sequenceName, firstLeftPart, region.getFirst()-1).toString();
			if(refSeq!=null) {
				int endReadSegment= Math.min(read.length(), region.getFirst()-first+5);
				String readSegment = read.substring(0,endReadSegment);
				//System.out.println(refSeq);
				//System.out.println(readSegment);
				String [] alignmentLeft;
				synchronized (alignerSTRsLeft) {
					alignmentLeft = alignerSTRsLeft.getAlignment(readSegment, refSeq.toString());
				}
				leftMismatches = countMismatches(alignmentLeft);
				encodedLeftAln = ReadAlignment.encodePairwiseAlignment(alignmentLeft);
				int lastCode = encodedLeftAln.getLast();
				if (leftMismatches<=readSegment.length()/10 && ReadAlignment.getOperator(lastCode)==ReadAlignment.ALIGNMENT_INSERTION) {
					softClipLeft = ReadAlignment.getOperationLength(lastCode);
					encodedLeftAln.removeLast();
				} else {
					encodedLeftAln = null;
				}
				softClipLeft+=(read.length()-endReadSegment);
			}	
		}
		if(last>region.getLast()+5) {
			CharSequence refSeq = fMIndex.getSequence(sequenceName, region.getLast()+1, last);
			if(refSeq!=null) {
				int startReadSegment= Math.max(0, read.length()-(last-region.getLast())-5);
				String readSegment = read.substring(startReadSegment);
				//System.out.println(refSeq);
				//System.out.println(readSegment);
				String [] alignmentRight;
				synchronized (alignerSTRsRight) {
					alignmentRight = alignerSTRsRight.getAlignment(readSegment, refSeq.toString());
				}
				rightMismatches = countMismatches(alignmentRight);
				encodedRightAln = ReadAlignment.encodePairwiseAlignment(alignmentRight);
				int firstCode = encodedRightAln.getFirst();
				if (rightMismatches<=readSegment.length()/10 && ReadAlignment.getOperator(firstCode)==ReadAlignment.ALIGNMENT_INSERTION) {
					softClipRight = ReadAlignment.getOperationLength(firstCode);
					encodedRightAln.removeFirst();
				} else {
					encodedRightAln=null;
				}
				softClipRight+=startReadSegment;
				
			}	
		}
		if(encodedLeftAln==null && encodedRightAln ==null) {
			return null;
		}
		if(encodedRightAln==null) {
			//Left alignment with right soft clip
			if(softClipLeft>0) encodedLeftAln.add(ReadAlignment.getAlnValue(softClipLeft, ReadAlignment.ALIGNMENT_SKIPFROMREAD));
			ReadAlignment aln = buildAln(read, sequenceName, firstLeftPart, region.getFirst()-1, encodedLeftAln);
			if(aln==null) return null;
			//System.out.println("Left alignment new genomic coordinates : "+first+"-"+last+" CIGAR:" +aln.getCigarString());
			if (!aln.clipBorders(kmerLength)) return null;
			aln.setAlignmentQuality((byte)(90-5*leftMismatches));
			aln.setNumMismatches((short) leftMismatches);
			return aln;
		}
		if(encodedLeftAln==null) {
			//Right alignment with left soft clip
			first = region.getLast()+1;
			if(softClipRight>0) encodedRightAln.addFirst(ReadAlignment.getAlnValue(softClipRight, ReadAlignment.ALIGNMENT_SKIPFROMREAD));
			ReadAlignment aln = buildAln(read, sequenceName, region.getLast()+1, last, encodedRightAln);
			if(aln==null) return null;
			//System.out.println("Right alignment new genomic coordinates : "+first+"-"+last+" CIGAR:" +aln.getCigarString());
			if (!aln.clipBorders(kmerLength)) return null;
			aln.setAlignmentQuality((byte)(90-5*rightMismatches));
			aln.setNumMismatches((short) rightMismatches);
			return aln;
		}
		first = firstLeftPart;
		int alignedLeft = read.length()-softClipLeft;
		int alignedRight = read.length()-softClipRight;
		int middleLength = read.length()-alignedLeft-alignedRight;
		if(middleLength<0) return null;
		int difference = region.length()-middleLength;
		//System.out.println("Aligned left: "+alignedLeft+" aligned right: "+alignedRight+" middle length: "+middleLength+" difference: "+difference);
		LinkedList<Integer> alignmentList = new LinkedList<Integer>();
		alignmentList.addAll(encodedLeftAln);
		if(difference>0) {
			// Region length > middle length. Add deletion
			alignmentList.add(ReadAlignment.getAlnValue(difference, ReadAlignment.ALIGNMENT_DELETION));
			if(middleLength>0) alignmentList.add(ReadAlignment.getAlnValue(middleLength, ReadAlignment.ALIGNMENT_MATCH));
		} else if (difference<0) {
			// Region length < middle length. Add insertion
			alignmentList.add(ReadAlignment.getAlnValue(-difference, ReadAlignment.ALIGNMENT_INSERTION));
			if(region.length()>0) alignmentList.add(ReadAlignment.getAlnValue(region.length(), ReadAlignment.ALIGNMENT_MATCH));
		} else {
			if(middleLength>0) alignmentList.add(ReadAlignment.getAlnValue(middleLength, ReadAlignment.ALIGNMENT_MATCH));
		}
		
		alignmentList.addAll(encodedRightAln);
		short mismatches = (short) (leftMismatches+rightMismatches);
		ReadAlignment aln = buildAln(read, sequenceName, first, last, alignmentList);
		if(aln==null) return null;
		//System.out.println("Building alignment from first "+first+" last: "+last+" softClipLeft: "+softClipLeft+" softClip right "+softClipRight+" cigar "+aln.getCigarString());
		if (!aln.clipBorders(kmerLength)) return null;
		aln.setAlignmentQuality((byte)(100-5*mismatches));
		aln.setNumMismatches(mismatches);
		return aln;
	}
	
	public List<ReadAlignment> fewMismatchesSingleStrandSearch(String query, int maxMismatches) {
		List<ReadAlignment> alns = new ArrayList<ReadAlignment>();
		List<Integer> alignment = new ArrayList<Integer>(1);
		alignment.add(ReadAlignment.getAlnValue(query.length(), ReadAlignment.ALIGNMENT_MATCH));
		// Whole read exact search
		List<UngappedSearchHit> readHits=fMIndex.exactSearch(query);
		for(int i=0;i<readHits.size()&& i<maxAlnsPerRead;i++) {
			UngappedSearchHit hit = readHits.get(i);
			ReadAlignment aln = buildAln (query, hit.getSequenceName(), hit.getStart()+1, hit.getStart()+query.length(), alignment);
			if(aln==null) continue;
			aln.setAlignmentQuality((byte) 100);
			alns.add(aln);
		}
		if (alns.size()>0) return alns;
		//One mismatch search
		int middle = query.length()/2;
		if(middle < 50) return alns;
		String firstPart = query.substring(0,middle);
		readHits=fMIndex.exactSearch(firstPart);
		for(UngappedSearchHit hit: readHits) {
			ReadAlignment aln = buildAln (query, hit.getSequenceName(), hit.getStart()+1, hit.getStart()+query.length(), alignment);
			if (aln==null) continue;
			int[] mismatches = countMismatches(query, aln);
			if(mismatches==null) continue;
			//System.out.println("first half. Next aln: "+aln.getSequenceName()+":"+aln.getFirst()+" mismatches: "+mismatches[0]+" CIGAR: "+aln.getCigarString()+" clip: "+mismatches[1]+" "+mismatches[2]);
			if(mismatches[0]<=maxMismatches) {
				if (mismatches[1]+mismatches[2]>0) {
					aln = buildAln(query, hit.getSequenceName(), hit.getStart()+1+mismatches[1], hit.getStart()+query.length()-mismatches[2], encodeAlignment(query.length(),mismatches));
				}
				if(aln==null) continue;
				aln.setAlignmentQuality((byte) (100-5*mismatches[0]));
				aln.setNumMismatches((short) mismatches[0]);
				
				alns.add(aln);
				//Best alignments selected during the filtering step
				if(alns.size()>=3*maxAlnsPerRead) break;
			}
		}
		if (alns.size()>0) return alns;
		String secondPart = query.substring(middle);
		readHits=fMIndex.exactSearch(secondPart);
		for(UngappedSearchHit hit: readHits) {
			int start = hit.getStart()-middle;
			if(start<0) continue;
			ReadAlignment aln = buildAln (query, hit.getSequenceName(), start+1, start+query.length(), alignment);
			if (aln==null) continue;
			int[] mismatches = countMismatches(query, aln);
			if(mismatches==null) continue;
			//System.out.println("Second half. Next aln: "+aln.getSequenceName()+":"+aln.getFirst()+" mismatches: "+mismatches[0]+" CIGAR: "+aln.getCigarString()+" clip: "+mismatches[1]+" "+mismatches[2]);
			if(mismatches[0]<=maxMismatches) {
				if (mismatches[1]+mismatches[2]>0) {
					aln = buildAln(query, hit.getSequenceName(), start+1+mismatches[1], start+query.length()-mismatches[2], encodeAlignment(query.length(),mismatches));
				}
				if(aln==null) continue;
				aln.setAlignmentQuality((byte) (100-5*mismatches[0]));
				aln.setNumMismatches((short) mismatches[0]);
				alns.add(aln);
				//Best alignments selected during the filtering step
				if(alns.size()>=3*maxAlnsPerRead) break;
			}
		}
		return alns;
	}
	private ReadAlignment buildAln(CharSequence query, String sequenceName, int first, int last, List<Integer> alignmentCodes) {
		if(first <=0) return null;
		ReadAlignment aln = new ReadAlignment(sequenceName, first, last, query.length(), 0);
		aln.setReadCharacters(query);
		if(alignmentCodes!=null)aln.setAlignment(alignmentCodes);
		//verify last exists
		if(!fMIndex.isValidPosition(sequenceName,aln.getLast())) return null;
		return aln;
	}
	public void loadSTRsFile(String strsFile) throws IOException {
		SimpleGenomicRegionFileHandler handler = new SimpleGenomicRegionFileHandler();
		knownSTRs=handler.loadRegionsAsMap(strsFile);
		flat();
	}

	private void flat() {
		Set<String>keys=knownSTRs.keySet();
		Iterator<String>it=keys.iterator();
		while(it.hasNext()) {
			String key = it.next();
			List<GenomicRegion> l =knownSTRs.get(key);
			List<GenomicRegion> newList = flat(l);
			while(isOverlappging(newList)) {
				newList=flat(newList);
			}
			knownSTRs.put(key, newList);
		}
	}
	private boolean isOverlappging(List<GenomicRegion>l) {
		for (int i = 0; i < l.size()-1; i++) {
			GenomicRegion current = l.get(i);
			GenomicRegion next = l.get(i+1);
			if(isOverlappgingSorted(current, next)) {
				return true;
			}
		}
		return false;
	}

	private boolean isOverlappgingSorted(GenomicRegion a,GenomicRegion b) {
		return b.getFirst()>=a.getFirst()&&b.getFirst()<=a.getLast();
	}

	private List<GenomicRegion> flat(List<GenomicRegion> l) {
		List<GenomicRegion> newList = new ArrayList<GenomicRegion>();
		for (int i = 0; i < l.size(); i++) {
			if(i+1==l.size())newList.add(l.get(i));
			else {
				GenomicRegion current =l.get(i);
				GenomicRegion next =l.get(i+1);
				if(isOverlappgingSorted(current, next)) {
					newList.add(new GenomicRegionImpl(current.getSequenceName(), current.getFirst(), Math.max(current.getLast(), next.getLast())));
					i++;
				}
				else {
					newList.add(current);
				}
			}	
		}

		return newList;
	}
}
