package ngsep.assembly;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import ngsep.alignments.ReadAlignment;
import ngsep.discovery.AlignmentsPileupGenerator;
import ngsep.discovery.IndelRealignerPileupListener;
import ngsep.discovery.VariantPileupListener;
import ngsep.genome.GenomicRegionPositionComparator;
import ngsep.genome.ReferenceGenome;
import ngsep.sequences.DNASequence;
import ngsep.sequences.QualifiedSequence;
import ngsep.variants.CalledGenomicVariant;

public class ConsensusBuilderBidirectionalWithPolishing implements ConsensusBuilder {
	int tolerance;
	
	boolean startConsensus = true;
	
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
				ReadAlignment alnRead = alignRead(rawConsensus, read, startConsensus);
				alignments.add(alnRead);
				List<AssemblyEmbedded> embeddedList = graph.getEmbedded(vertexPreviousEdge.getIndex());
				for(AssemblyEmbedded embedded:embeddedList) {
					CharSequence embeddedRead = embedded.getRead();
					boolean reverseE = (reverse!=embedded.isReverse());
					if(reverseE) embeddedRead = DNASequence.getReverseComplement(embeddedRead.toString());
					ReadAlignment alnEmbedded = alignRead(rawConsensus, read, startConsensus+embedded.getStartPosition());
					alignments.add(alnEmbedded);
				}
			}
			lastVertex = vertexNextEdge;
		}
		System.out.println(pathS);
		String consensus = rawConsensus.toString();
		List<CalledGenomicVariant> variants = callVariants(consensus,alignments);
		return applyVariants(consensus, variants);
	}

	private ReadAlignment alignRead(StringBuilder rawConsensus, CharSequence read, int startConsensus) {
		// TODO Auto-generated method stub
		return null;
	}
	private List<CalledGenomicVariant> callVariants(String consensus, List<ReadAlignment> alignments) {
		AlignmentsPileupGenerator generator = new AlignmentsPileupGenerator();
		String refId = "Consensus";
		ReferenceGenome genome = new ReferenceGenome(new QualifiedSequence(refId, consensus));
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
		for(CalledGenomicVariant call:variants) {
			if(call.isUndecided()|| call.isHeterozygous()) continue;
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
		return new DNASequence(polishedConsensus.toString());
	}
}
