package ngsep.assembly;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import ngsep.alignments.ReadAlignment;
import ngsep.discovery.AlignmentsPileupGenerator;
import ngsep.discovery.IndelRealignerPileupListener;
import ngsep.discovery.VariantPileupListener;
import ngsep.genome.GenomicRegionPositionComparator;
import ngsep.genome.ReferenceGenome;
import ngsep.sequences.DNASequence;
import ngsep.sequences.FMIndexUngappedSearchHit;
import ngsep.sequences.KmerHitsCluster;
import ngsep.sequences.KmersExtractor;
import ngsep.sequences.PairwiseAlignmentAffineGap;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.RawRead;
import ngsep.variants.CalledGenomicVariant;

public class ConsensusBuilderBidirectionalWithPolishing implements ConsensusBuilder {
	private int tolerance;
	
	private static final String MOCK_REFERENCE_NAME = "Consensus";
	
	private PairwiseAlignmentAffineGap aligner = new PairwiseAlignmentAffineGap(1, 2, 1, 1);
	
	
	@Override
	public List<CharSequence> makeConsensus(AssemblyGraph graph) 
	{
		//List of final contigs
		List<CharSequence> consensusList = new ArrayList<CharSequence>();
		List<List<AssemblyEdge>> paths = graph.getPaths(); 
		for(int i = 0; i < paths.size(); i++)
		{
			List<AssemblyEdge> path = paths.get(i);
			CharSequence consensusPath = makeConsensus (graph, path);
			consensusList.add(consensusPath);
		}
		
		return consensusList;
	}
	
	private CharSequence makeConsensus(AssemblyGraph graph, List<AssemblyEdge> path) {
		StringBuilder rawConsensus = new StringBuilder();
		AssemblyVertex lastVertex = null;
		List<ReadAlignment> alignments = new ArrayList<ReadAlignment>();
		int totalReads = 0;
		int unalignedReads = 0;
		String pathS = "";
		if(path.size()==1) {
			rawConsensus.append(path.get(0).getVertex1().getRead());
			return rawConsensus;
		}
		for(int j = 0; j < path.size(); j++) {
			//Needed to find which is the origin vertex
			AssemblyEdge edge = path.get(j);
			AssemblyVertex vertexPreviousEdge;
			AssemblyVertex vertexNextEdge;
			//If the first edge is being checked, compare to the second edge to find the origin vertex
			if(j == 0) {
				AssemblyEdge nextEdge = path.get(j + 1);
				vertexNextEdge = edge.getSharedVertex(nextEdge);
				if(vertexNextEdge== null) throw new RuntimeException("Inconsistency found in first edge of path");
				vertexPreviousEdge = edge.getVertex1();
				if(vertexPreviousEdge == vertexNextEdge) vertexPreviousEdge = edge.getVertex2();
			}
			else if (lastVertex == edge.getVertex1()) {
				vertexPreviousEdge = edge.getVertex1();
				vertexNextEdge = edge.getVertex2();
			}
			else if (lastVertex == edge.getVertex2()) {
				vertexPreviousEdge = edge.getVertex2();
				vertexNextEdge = edge.getVertex1();
			}
			else {
				throw new RuntimeException("Inconsistency found in path");
			}
			if(j == 0) {
				pathS = pathS.concat(vertexPreviousEdge.getIndex() + ",");
				String seq = vertexPreviousEdge.getRead().toString();
				boolean reverse = !vertexPreviousEdge.isStart();
				if(reverse) seq = DNASequence.getReverseComplement(seq);
				rawConsensus.append(seq.toUpperCase());
			} 
			else if(vertexPreviousEdge.getRead()!=vertexNextEdge.getRead()) {
				//If the second string is not start, then the reverse complement is added to the consensus
				String seq = vertexNextEdge.getRead().toString();
				boolean reverse = !vertexNextEdge.isStart();
				if(reverse) seq = DNASequence.getReverseComplement(seq);
				int overlap = edge.getOverlap();
				if(seq.length() - overlap > tolerance) {
					pathS = pathS.concat(vertexNextEdge.getIndex() + ",");
					//String overlapSegment = nextSequence.substring(0, edge.getOverlap());
					String remainingSegment = seq.substring(edge.getOverlap());
					rawConsensus.append(remainingSegment.toUpperCase());
					
				} 
				else 
				{
					System.err.println("Non embedded edge has overlap: "+edge.getOverlap()+ " and length: "+seq.length());
				}
				
			}
			if(vertexPreviousEdge.getRead()==vertexNextEdge.getRead()) {
				CharSequence read = vertexPreviousEdge.getRead();
				boolean reverse = !vertexPreviousEdge.isStart();
				if(reverse) read = DNASequence.getReverseComplement(read.toString());
				int startConsensus = Math.max(0, rawConsensus.length() - read.length()-10);
				String subjectSequence = rawConsensus.substring(startConsensus);
				Map<CharSequence, Integer> uniqueKmersSubject = extractUniqueKmers(subjectSequence);
				totalReads++;
				ReadAlignment alnRead = alignRead(rawConsensus, read, startConsensus, uniqueKmersSubject, 0);
				if (alnRead!=null) alignments.add(alnRead);
				else unalignedReads++;
				
				List<AssemblyEmbedded> embeddedList = graph.getEmbedded(vertexPreviousEdge.getIndex());
				for(AssemblyEmbedded embedded:embeddedList) {
					CharSequence embeddedRead = embedded.getRead();
					boolean reverseE = (reverse!=embedded.isReverse());
					if(reverseE) embeddedRead = DNASequence.getReverseComplement(embeddedRead.toString());
					totalReads++;
					ReadAlignment alnEmbedded = alignRead(rawConsensus, embeddedRead, startConsensus, uniqueKmersSubject, embedded.getStartPosition());
					if(alnEmbedded!=null) alignments.add(alnEmbedded);
					else unalignedReads++;
				}
				System.out.println("Total reads: "+totalReads+" alignments: "+alignments.size()+" unaligned: "+unalignedReads);
			}
			lastVertex = vertexNextEdge;
		}
		System.out.println(pathS);
		String consensus = rawConsensus.toString();
		List<CalledGenomicVariant> variants = callVariants(consensus,alignments);
		System.out.println("Identified "+variants.size()+" total variants");
		return applyVariants(consensus, variants);
	}

	private Map<CharSequence, Integer> extractUniqueKmers(CharSequence sequence) {
		CharSequence [] rawKmers = KmersExtractor.extractKmers(sequence, 15, 1, true, true, true);
		Map<CharSequence, Integer> answer = new LinkedHashMap<CharSequence, Integer>();
		Map<CharSequence, Integer> reverseMap = new HashMap<CharSequence,Integer>();
		boolean [] unique = new boolean[rawKmers.length];
		Arrays.fill(unique, true);
		for(int i=0;i<rawKmers.length;i++) {
			CharSequence kmer = rawKmers[i];
			if(kmer == null) {
				unique[i] = false;
				continue;
			}
			Integer previousStart = reverseMap.get(kmer);
			if(previousStart!=null) {
				unique[previousStart] = false;
				unique[i] = false;
				continue;
			}
			reverseMap.put(kmer,i);
		}
		for(int i=0;i<rawKmers.length;i++) {
			CharSequence kmer = rawKmers[i];
			if(kmer!=null && unique[i]) answer.put(kmer,i);
		}
		return answer;
	}

	private ReadAlignment alignRead(StringBuilder rawConsensus, CharSequence read, int startConsensus, Map<CharSequence, Integer> uniqueKmersSubject, int offsetRead) {
		Map<CharSequence, Integer> uniqueKmersRead = extractUniqueKmers(read);
		List<FMIndexUngappedSearchHit> initialKmerHits = new ArrayList<FMIndexUngappedSearchHit>();
		for(CharSequence kmerRead:uniqueKmersRead.keySet()) {
			Integer subjectPos = uniqueKmersSubject.get(kmerRead);
			if(subjectPos==null) continue;
			subjectPos+=startConsensus;
			FMIndexUngappedSearchHit hit = new FMIndexUngappedSearchHit(kmerRead.toString(), "", subjectPos);
			hit.setQueryIdx(uniqueKmersRead.get(kmerRead));
			initialKmerHits.add(hit);
		}
		if(initialKmerHits.size()==0) return null;
		List<KmerHitsCluster> clusters = GraphBuilderFMIndex.clusterSequenceKmerAlns(read, initialKmerHits);
		Collections.sort(clusters, (o1,o2)->o2.getNumDifferentKmers()-o1.getNumDifferentKmers());
		if(clusters.size()>1) {
			int l0 = clusters.get(0).getNumDifferentKmers();
			int l1 = clusters.get(1).getNumDifferentKmers();
			if(l0<0.9*initialKmerHits.size()) {
				System.out.println("Number of clusters: "+clusters.size()+" best cluster length: "+l0+" start: "+clusters.get(0).getFirst()+" end: "+clusters.get(0).getLast());
				System.out.println("Second cluster length: "+l1+" start: "+clusters.get(1).getFirst()+" end: "+clusters.get(1).getLast());
				return null;
			}	
		}
		return buildCompleteAlignment(rawConsensus, read.toString(), clusters.get(0));
	}
	private ReadAlignment buildCompleteAlignment(StringBuilder subject, String query, KmerHitsCluster kmerHitsCluster) {
		List<FMIndexUngappedSearchHit> kmerHits = kmerHitsCluster.getHitsByQueryIdx();
		
		int clusterFirst = kmerHitsCluster.getFirst();
		int subjectNext = Math.max(0, clusterFirst-1);
		System.out.println("Consensus current length: "+subject.length()+". Next query length: "+query.length()+" kmer hits: "+kmerHits.size()+" subject next: "+subjectNext);
		int queryNext = 0;
		int alnStart = -1;
		char matchOp = ReadAlignment.ALIGNMENT_CHAR_CODES.charAt(ReadAlignment.ALIGNMENT_MATCH);
		char insertionOp = ReadAlignment.ALIGNMENT_CHAR_CODES.charAt(ReadAlignment.ALIGNMENT_INSERTION);
		char deletionOp = ReadAlignment.ALIGNMENT_CHAR_CODES.charAt(ReadAlignment.ALIGNMENT_DELETION);
		char softClipOp = ReadAlignment.ALIGNMENT_CHAR_CODES.charAt(ReadAlignment.ALIGNMENT_SKIPFROMREAD);
		StringBuilder cigar = new StringBuilder();
		int nextMatchLength = 0;
		for(FMIndexUngappedSearchHit kmerHit:kmerHits) {
			//System.out.println("Processing Kmer hit at pos: "+kmerHit.getQueryIdx()+" query next: "+queryNext+" subject next: "+subjectNext+" subject hit start: "+kmerHit.getStart()+" cigar length: "+cigar.length());
			int kmerLength = kmerHit.getQuery().length();
			if(alnStart==-1) {
				alnStart = kmerHit.getStart();
				if(kmerHit.getQueryIdx()>0) cigar.append(""+kmerHit.getQueryIdx()+""+softClipOp);
				nextMatchLength+=kmerLength;
				subjectNext = kmerHit.getStart()+kmerLength;
				queryNext = kmerHit.getQueryIdx()+kmerLength;
			} else if(kmerHit.getQueryIdx() >= queryNext && subjectNext<=kmerHit.getStart()) {
				//Kmer does not overlap with already aligned segments
				String subjectStr = subject.substring(subjectNext,kmerHit.getStart());
				String queryStr = query.substring(queryNext,kmerHit.getQueryIdx());
				if(subjectStr.length()==queryStr.length() && subjectStr.length()<10) {
					nextMatchLength+=subjectStr.length();
				} else {
					if(nextMatchLength>0 && (subjectStr.length()>0 || queryStr.length()>0)) {
						System.out.println("Found internal segment for possible alignment. Subject length "+subjectStr.length()+" query length "+queryStr.length()+" current match length: "+nextMatchLength);
						cigar.append(""+nextMatchLength+""+matchOp);
						nextMatchLength = 0;
					}
					if(subjectStr.length()>0 && queryStr.length()>0) {
						System.out.println("Aligning segment of length "+subjectStr.length()+" of subject with total length: "+subject.length()+" to segment with length "+queryStr.length()+" of query with total length: "+query.length());
						String [] alignedFragments = aligner.getAlignment(subjectStr, queryStr);
						String cigarSegment = buildCigar(alignedFragments[0],alignedFragments[1]);
						cigar.append(cigarSegment);
						
					} else if (subjectStr.length()>0) {
						cigar.append(""+subjectStr.length()+""+deletionOp);
					} else if (queryStr.length()>0) {
						cigar.append(""+queryStr.length()+""+insertionOp);
					}
				}
				nextMatchLength+=kmerLength;
				subjectNext = kmerHit.getStart()+kmerLength;
				queryNext = kmerHit.getQueryIdx()+kmerLength;
			}
			
			//System.out.println("Processed Kmer hit at pos: "+kmerHit.getQueryIdx()+" query next: "+queryNext+" subject next: "+subjectNext);
		}
		if(nextMatchLength>0) {
			cigar.append(""+nextMatchLength+""+matchOp);
			nextMatchLength = 0;
		}
		int alnFirst = alnStart+1;
		int alnLast = subjectNext;
		System.out.println("Aligned almost complete query. first: "+alnFirst+" last: "+alnLast+" CIGAR: "+cigar+" query next: "+queryNext+" query length: "+query.length());
		if(queryNext<query.length()) {
			String queryStr = query.substring(queryNext);
			//TODO: check if it is worth to align the sequence end
			/*if(subjectNext<subject.length()) {
				String subjectStr = subject.substring(subjectNext,Math.min(clusterLast, subject.length()));
				//TODO: Softclip if it starts with indel
				System.out.println("Aligning end of length "+subjectStr.length()+" of subject subsequence with total length: "+subject.length()+" to end with length "+queryStr.length()+" of query with total length: "+query.length());
				String [] alignedFragments = aligner.getAlignment(subjectStr, queryStr);
				String cigarSegment = buildCigar(alignedFragments[0],alignedFragments[1]);
				
				cigar.append(cigarSegment);	
			} else {*/
				//Ignore last bp
				cigar.append(""+queryStr.length()+""+softClipOp);
			//}
		}
		ReadAlignment finalAlignment = new ReadAlignment(MOCK_REFERENCE_NAME, alnStart+1, alnLast, query.length(), 0);
		finalAlignment.setReadCharacters(query);
		finalAlignment.setQualityScores(RawRead.generateFixedQSString('5', query.length()));
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
			char op = ReadAlignment.ALIGNMENT_CHAR_CODES.charAt(ReadAlignment.ALIGNMENT_MATCH);
			if(subjectChar == AlignmentConstantGap.GAP_CHARACTER) {
				op = ReadAlignment.ALIGNMENT_CHAR_CODES.charAt(ReadAlignment.ALIGNMENT_INSERTION);
			} else if(queryChar == AlignmentConstantGap.GAP_CHARACTER) {
				op = ReadAlignment.ALIGNMENT_CHAR_CODES.charAt(ReadAlignment.ALIGNMENT_DELETION);
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
		if (nextLength>0) cigar.append(""+nextLength+""+nextOperator);
		return cigar.toString();
	}
	private List<CalledGenomicVariant> callVariants(String consensus, List<ReadAlignment> alignments) {
		AlignmentsPileupGenerator generator = new AlignmentsPileupGenerator();
		ReferenceGenome genome = new ReferenceGenome(new QualifiedSequence(MOCK_REFERENCE_NAME, consensus));
		generator.setSequencesMetadata(genome.getSequencesMetadata());
		IndelRealignerPileupListener realignerListener = new IndelRealignerPileupListener();
		realignerListener.setGenome(genome);
		generator.addListener(realignerListener);
		VariantPileupListener varListener = new VariantPileupListener();
		varListener.setGenome(genome);
		generator.addListener(varListener);
		Collections.sort(alignments, GenomicRegionPositionComparator.getInstance());
		for(ReadAlignment aln:alignments) generator.processAlignment(aln);
		return varListener.getCalledVariants();
	}

	private CharSequence applyVariants(String consensus, List<CalledGenomicVariant> variants) {
		StringBuilder polishedConsensus = new StringBuilder();
		int l = consensus.length();
		int nextPos = 1;
		int appliedVariants = 0;
		for(CalledGenomicVariant call:variants) {
			if(call.isUndecided()|| call.isHeterozygous()) continue;
			appliedVariants++;
			String [] calledAlleles = call.getCalledAlleles();
			if(nextPos<call.getFirst()) {
				//Fill haplotypes with non variant segment
				String segment = consensus.substring(nextPos-1, call.getFirst()-1);
				polishedConsensus.append(segment);
			}
			polishedConsensus.append(calledAlleles[0]);
			nextPos = call.getLast()+1;
		}
		if(nextPos<=l) {
			//End of a chromosome
			CharSequence nonVarLast = consensus.substring(nextPos-1);
			polishedConsensus.append(nonVarLast);
		}
		System.out.println("Applied "+appliedVariants+" homozygous alternative variants");
		return new DNASequence(polishedConsensus.toString());
	}
}
