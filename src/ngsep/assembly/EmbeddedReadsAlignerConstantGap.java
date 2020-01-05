package ngsep.assembly;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Set;

import ngsep.alignments.KmerAlignmentCluster;
import ngsep.alignments.KmerWithStart;
import ngsep.alignments.ReadAlignment;
import ngsep.genome.GenomicRegionPositionComparator;
import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.FMIndexSingleSequence;

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
		List<ReadAlignment> initialKmerAlns = new ArrayList<>();
		for (KmerWithStart kmer:kmers) {
			Set<Integer> kmerAlns=fmIndex.search(kmer.getKmer().toString());
			//if(querySequenceId==1) System.out.println("Query: "+querySequenceId+" Found "+kmerAlns.size()+" alignments for kmer: "+kmer.getKmer().toString());
			for(int start:kmerAlns) {
				//if(querySequenceId==1) System.out.println("Kmer start: "+kmer.getStart()+" Next alignment: "+aln.getSequenceIndex()+": "+aln.getFirst()+"-"+aln.getLast()+" rc: "+aln.isNegativeStrand());
				ReadAlignment kmerAln = new ReadAlignment("", start, start+kmerLength-1, kmerLength, 0);
				kmerAln.setReadNumber(kmer.getStart());
				initialKmerAlns.add(kmerAln);
			}
		}
		if(initialKmerAlns.size()==0) return null;
		Collections.sort(initialKmerAlns,GenomicRegionPositionComparator.getInstance());
		List<KmerAlignmentCluster> clusters = GraphBuilderFMIndex.clusterSequenceKmerAlns(query, initialKmerAlns);
		
		Collections.sort(clusters, new Comparator<KmerAlignmentCluster>() {
			@Override
			public int compare(KmerAlignmentCluster o1, KmerAlignmentCluster o2) {
				return o2.getNumDifferentKmers()-o1.getNumDifferentKmers();
			}
		});
		return buildCompleteAlignment(subject,query,clusters.get(0));
		
	}

	private ReadAlignment buildCompleteAlignment(String subject, String query, KmerAlignmentCluster kmerAlignmentCluster) {
		List<ReadAlignment> kmerAlns = kmerAlignmentCluster.getAlignmentsByReadNumber();
		
		
		int clusterFirst = kmerAlignmentCluster.getFirst();
		int clusterLast = kmerAlignmentCluster.getLast();
		int subjectNext = Math.max(0, clusterFirst);
		int queryNext = 0;
		int alnFirst = -1;
		
		StringBuilder cigar = new StringBuilder();
		for(ReadAlignment aln:kmerAlns) {
			int kmerLength = aln.getReadLength();
			if(aln.getReadNumber()==0) {
				alnFirst = aln.getFirst();
				cigar.append(""+kmerLength+"M");
				subjectNext = alnFirst+kmerLength;
				queryNext = kmerLength;
			} else if(subjectNext<=aln.getFirst()) {
				//Possibly align left part
				String subjectStr = subject.substring(subjectNext,aln.getFirst());
				String queryStr = query.substring(queryNext,aln.getReadNumber());
				if(subjectStr.length()>0 && queryStr.length()>0) {
					String [] alignedFragments = aligner.getAlignment(subjectStr, queryStr);
					cigar.append(buildCigar(alignedFragments[0],alignedFragments[1]));
					if(alnFirst==-1) alnFirst = aln.getFirst()-subjectStr.length();
				} else if (subjectStr.length()>0) {
					cigar.append(""+subjectStr.length()+"D");
				}  else if (queryStr.length()>0) {
					if(alnFirst==-1) {
						cigar.append(""+queryStr.length()+"S");
						alnFirst = aln.getFirst();
					} else {
						cigar.append(""+queryStr.length()+"I");
					}
					
				}
				cigar.append(""+kmerLength+"M");
				subjectNext = aln.getFirst()+kmerLength;
				queryNext = aln.getReadNumber()+kmerLength;
			}
		}
		int alnLast = subjectNext-1;
		if(queryNext<query.length()) {
			String queryStr = query.substring(queryNext);
			if(subjectNext<subject.length()) {
				String subjectStr = subject.substring(subjectNext,Math.min(clusterLast, subject.length()));
				
				String [] alignedFragments = aligner.getAlignment(subjectStr, queryStr);
				cigar.append(buildCigar(alignedFragments[0],alignedFragments[1]));	
			} else {
				//Ignore last bp
				cigar.append(""+queryStr.length()+"S");
			}
		}
		ReadAlignment finalAlignment = new ReadAlignment(0, alnFirst, alnLast, query.length(), 0);
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
