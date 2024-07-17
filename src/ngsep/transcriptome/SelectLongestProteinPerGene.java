package ngsep.transcriptome;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;
import ngsep.sequences.io.FastaSequencesHandler;

public class SelectLongestProteinPerGene {

	public static void main(String[] args) throws Exception {
		FastaSequencesHandler handler = new FastaSequencesHandler();
		handler.setSequenceType(StringBuilder.class);
		List<QualifiedSequence> proteins = handler.loadSequences(args[0]);
		QualifiedSequenceList cds=null;
		if(args.length>1) {
			cds = new QualifiedSequenceList(handler.loadSequences(args[1]));
		}
		Map<String,QualifiedSequence> seqsByGene = new TreeMap<>();
		for(QualifiedSequence protein:proteins) {
			String geneId = protein.getComments();
			geneId = geneId.split(" |\t")[0];
			QualifiedSequence currentProtein = seqsByGene.get(geneId);
			if(currentProtein==null || currentProtein.getLength()<protein.getLength()) {
				CharSequence outSequence;
				if(cds!=null) outSequence = cds.get(protein.getName()).getCharacters();
				else outSequence = protein.getCharacters();
				seqsByGene.put(geneId, new QualifiedSequence(geneId,outSequence));
			}
		}
		handler.saveSequences(new ArrayList<>(seqsByGene.values()), System.out, 100);
	}

}
