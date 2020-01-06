package ngsep.assembly;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Set;

import ngsep.alignments.ReadAlignment;
import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.FMIndexSingleSequence;
import ngsep.sequences.FMIndexUngappedSearchHit;
import ngsep.sequences.KmerHitsCluster;
import ngsep.sequences.KmerWithStart;

/**
 * Alignment of a main read with its embedded reads
 * @author David Guevara
 * @author Jorge Duitama
 */
public class EmbeddedReadsAlignerConstantGap
{
	AlignmentConstantGap aligner;
	
	public EmbeddedReadsAlignerConstantGap(int match, int gap, int mismatch)
	{
		this.aligner = new AlignmentConstantGap(match, gap, mismatch);
	}
		
	/**
	 * Sets the alignment of the reads with the main read
	 * @param mainRead Main read that contains every embedded read
	 * @param embeddedReads List of embedded reads
	 * @param match Score for a match
	 * @param gap Score for a gap
	 * @param mismatch Score for a mismatch
	 * @return alignments Aligned strings
	 */
	public List<ReadAlignment> alignEmbeddedReads(String mainRead, List<AssemblyEmbedded> embeddedReads, boolean reverse)
	{
		List<ReadAlignment> alignmentsEmbedded = new ArrayList<ReadAlignment>();
		FMIndexSingleSequence fm = new FMIndexSingleSequence(mainRead);
		for(int i = 0; i < embeddedReads.size(); i++)
		{
			AssemblyEmbedded embeddedRead = embeddedReads.get(i);
			String embedded = embeddedRead.getRead().toString();
			//If the main and the embedded read have different senses, the embedded read must be reversed
			if(reverse != embeddedRead.isReverse())
			{
				embedded = DNAMaskedSequence.getReverseComplement(embedded);
			}
			
			ReadAlignment embeddedAln = getAlignment(fm, mainRead, embedded);
			if(embeddedAln!=null) {
				alignmentsEmbedded.add(embeddedAln);
			}
		}
		return alignmentsEmbedded;
 	}

	private ReadAlignment getAlignment(FMIndexSingleSequence fmIndex, String subject, String query) {
		int kmerLength = 15;
		List<KmerWithStart> kmers = KmerWithStart.selectKmers(query, kmerLength, kmerLength);
		if(kmers==null) return null;
		List<FMIndexUngappedSearchHit> initialKmerHits = new ArrayList<>();
		for (KmerWithStart kmer:kmers) {
			Set<Integer> kmerAlns=fmIndex.exactSearch(kmer.getKmer().toString());
			//if(querySequenceId==1) System.out.println("Query: "+querySequenceId+" Found "+kmerAlns.size()+" alignments for kmer: "+kmer.getKmer().toString());
			for(int start:kmerAlns) {
				//if(querySequenceId==1) System.out.println("Kmer start: "+kmer.getStart()+" Next alignment: "+aln.getSequenceIndex()+": "+aln.getFirst()+"-"+aln.getLast()+" rc: "+aln.isNegativeStrand());
				FMIndexUngappedSearchHit kmerHit = new FMIndexUngappedSearchHit(query, "", start);
				kmerHit.setQueryIdx(kmer.getStart());
				initialKmerHits.add(kmerHit);
			}
		}
		if(initialKmerHits.size()==0) return null;
		List<KmerHitsCluster> clusters = GraphBuilderFMIndex.clusterSequenceKmerAlns(query, initialKmerHits);
		
		Collections.sort(clusters, new Comparator<KmerHitsCluster>() {
			@Override
			public int compare(KmerHitsCluster o1, KmerHitsCluster o2) {
				return o2.getNumDifferentKmers()-o1.getNumDifferentKmers();
			}
		});
		return buildCompleteAlignment(subject,query,clusters.get(0));
		
	}

	private ReadAlignment buildCompleteAlignment(String subject, String query, KmerHitsCluster kmerHitsCluster) {
		List<FMIndexUngappedSearchHit> kmerHits = kmerHitsCluster.getHitsByQueryIdx();
		
		
		int clusterFirst = kmerHitsCluster.getFirst();
		int clusterLast = kmerHitsCluster.getLast();
		int subjectNext = Math.max(0, clusterFirst);
		int queryNext = 0;
		int alnStart = -1;
		
		StringBuilder cigar = new StringBuilder();
		for(FMIndexUngappedSearchHit kmerHit:kmerHits) {
			int kmerLength = kmerHit.getQuery().length();
			if(kmerHit.getQueryIdx()==0) {
				alnStart = kmerHit.getStart();
				cigar.append(""+kmerLength+"M");
				subjectNext = alnStart+kmerLength;
				queryNext = kmerLength;
			} else if(kmerHit.getQueryIdx() >= queryNext && subjectNext<=kmerHit.getStart()) {
				//Possibly align left part
				String subjectStr = subject.substring(subjectNext,kmerHit.getStart());
				String queryStr = query.substring(queryNext,kmerHit.getQueryIdx());
				if(subjectStr.length()>0 && queryStr.length()>0) {
					String [] alignedFragments = aligner.getAlignment(subjectStr, queryStr);
					String cigarSegment = buildCigar(alignedFragments[0],alignedFragments[1]);
					//TODO: Softclip if it starts with indel
					cigar.append(cigarSegment);
					if(alnStart==-1) alnStart = kmerHit.getStart()-subjectStr.length();
				} else if (subjectStr.length()>0) {
					cigar.append(""+subjectStr.length()+"D");
				}  else if (queryStr.length()>0) {
					if(alnStart==-1) {
						cigar.append(""+queryStr.length()+"S");
						alnStart = kmerHit.getStart();
					} else {
						cigar.append(""+queryStr.length()+"I");
					}
					
				}
				cigar.append(""+kmerLength+"M");
				subjectNext = kmerHit.getStart()+kmerLength;
				queryNext = kmerHit.getQueryIdx()+kmerLength;
			}
		}
		int alnLast = subjectNext;
		if(queryNext<query.length()) {
			String queryStr = query.substring(queryNext);
			if(subjectNext<subject.length()) {
				String subjectStr = subject.substring(subjectNext,Math.min(clusterLast, subject.length()));
				//TODO: Softclip if it starts with indel
				String [] alignedFragments = aligner.getAlignment(subjectStr, queryStr);
				String cigarSegment = buildCigar(alignedFragments[0],alignedFragments[1]);
				
				cigar.append(cigarSegment);	
			} else {
				//Ignore last bp
				cigar.append(""+queryStr.length()+"S");
			}
		}
		ReadAlignment finalAlignment = new ReadAlignment(0, alnStart+1, alnLast, query.length(), 0);
		finalAlignment.setReadCharacters(query);
		finalAlignment.setCigarString(cigar.toString());
		return finalAlignment;
	}

	private String buildCigar(String subjectAln, String queryAln) {
		StringBuilder cigar = new StringBuilder();
		char nextOperator = 0;
		int nextLength = 0;
		for(int i=0;i<subjectAln.length();i++) {
			char subjectChar = subjectAln.charAt(i);
			char queryChar = queryAln.charAt(i);
			char op = 'M';
			if(subjectChar == AlignmentConstantGap.GAP_CHARACTER) {
				op = 'I';
			} else if(queryChar == AlignmentConstantGap.GAP_CHARACTER) {
				op = 'D';
			}
			if(op != nextOperator) {
				if(nextLength>0) {
					cigar.append(""+nextLength+""+nextOperator);
				}
				nextOperator = op;
				nextLength = 0;
			}
			nextLength++;
		}
		return cigar.toString();
	}
}
