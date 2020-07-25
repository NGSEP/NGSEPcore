package ngsep.assembly;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import ngsep.alignments.MinimizersTableReadAlignmentAlgorithm;
import ngsep.alignments.ReadAlignment;
import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.QualifiedSequence;

public class ConsensusBuilderBidirectionalSimple implements ConsensusBuilder {
	int tolerance;
	
	boolean startConsensus = true;
	
	@Override
	public List<QualifiedSequence> makeConsensus(AssemblyGraph graph) 
	{
		//List of final contigs
		List<QualifiedSequence> consensusList = new ArrayList<QualifiedSequence>();
		List<List<AssemblyEdge>> paths = graph.getPaths(); 
		for(int i = 0; i < paths.size(); i++)
		{
			List<AssemblyEdge> path = paths.get(i);
			String sequenceName = "Contig_"+(i+1);
			CharSequence consensusPath = makeConsensus (graph, path, i, sequenceName);
			consensusList.add(new QualifiedSequence(sequenceName,consensusPath));
		}
		
		return consensusList;
	}
	
	private CharSequence makeConsensus(AssemblyGraph graph, List<AssemblyEdge> path, int sequenceIdx, String sequenceName) 
	{
		StringBuilder consensus = new StringBuilder();
		AssemblyVertex lastVertex = null;
		MinimizersTableReadAlignmentAlgorithm aligner = new MinimizersTableReadAlignmentAlgorithm();
		
		String pathS = "";
		if(path.size()==1) {
			consensus.append(path.get(0).getVertex1().getRead());
			return consensus;
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
			} else if (lastVertex == edge.getVertex1()) {
				vertexPreviousEdge = edge.getVertex1();
				vertexNextEdge = edge.getVertex2();
			} else if (lastVertex == edge.getVertex2()) {
				vertexPreviousEdge = edge.getVertex2();
				vertexNextEdge = edge.getVertex1();
			} else {
				throw new RuntimeException("Inconsistency found in path");
			}
			if(j == 0) {
				pathS = pathS.concat(vertexPreviousEdge.getUniqueNumber() + ",");
				CharSequence seq = vertexPreviousEdge.getRead().getCharacters();
				boolean reverse = !vertexPreviousEdge.isStart();
				if(reverse) seq = DNAMaskedSequence.getReverseComplement(seq);
				consensus.append(seq);
			} else if(vertexPreviousEdge.getRead()!=vertexNextEdge.getRead()) {
				// Augment consensus with the next path read
				CharSequence nextPathSequence = vertexNextEdge.getRead().getCharacters();
				boolean reverse = !vertexNextEdge.isStart();
				if(reverse) nextPathSequence = DNAMaskedSequence.getReverseComplement(nextPathSequence);
				//if (rawConsensus.length()>490000 && rawConsensus.length()<530000) printAllOverlappingSeqs(graph,path,j,vertexPreviousEdge);
				
				//int startSuffix = edge.getOverlap();
				Map<CharSequence, Integer> uniqueKmersSubject = aligner.extractUniqueKmers(consensus,Math.max(0, consensus.length()-nextPathSequence.length()),consensus.length());
				ReadAlignment alnRead = aligner.alignRead(sequenceIdx, consensus, nextPathSequence, uniqueKmersSubject, 0.5);
				int startSuffix;
				if(alnRead!=null) {
					alnRead.setSequenceName(sequenceName);
					alnRead.setReadName(vertexNextEdge.getRead().getName());
					int posAlnRead = nextPathSequence.length()-1-alnRead.getSoftClipEnd();
					int lastPosSubject = alnRead.getReferencePositionAlignedRead(posAlnRead);
					//Just in case cycle but if the read aligns this should not enter
					while(posAlnRead>0 && lastPosSubject<0) {
						posAlnRead--;
						lastPosSubject = alnRead.getReferencePositionAlignedRead(posAlnRead);
					}
					if(lastPosSubject>=0) {
						startSuffix = posAlnRead + (consensus.length()-lastPosSubject+1);
					} else {
						startSuffix = edge.getOverlap();
					}
					//System.out.println("Calculated overlap from alignment: "+startSuffix+" alignment: "+alnRead+" edge: "+edge );
				} else {
					startSuffix = edge.getOverlap();
				}
				if(startSuffix<nextPathSequence.length()) {
					pathS = pathS.concat(vertexNextEdge.getUniqueNumber() + ",");
					String remainingSegment = nextPathSequence.subSequence(startSuffix, nextPathSequence.length()).toString();
					//if (consensus.length()>490000 && consensus.length()<510000) System.out.println("Consensus length: "+consensus.length()+" Vertex: "+vertexNextEdge.getUniqueNumber()+" read length: "+seq.length()+" overlap: "+edge.getOverlap()+" remaining: "+remainingSegment.length());
					consensus.append(remainingSegment.toUpperCase());
				}
			}
			lastVertex = vertexNextEdge;
		}
		System.out.println(pathS);
		return consensus;
	}
}
