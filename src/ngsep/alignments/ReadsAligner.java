package ngsep.alignments;

import java.util.List;

import ngsep.genome.GenomicRegion;
import ngsep.main.CommandsDescriptor;
import ngsep.sequences.FMIndex;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;
import ngsep.sequences.io.FastaSequencesHandler;

public class ReadsAligner {

	public static void main(String[] args) throws Exception 
	{
		FMIndex fMIndex = new FMIndex();
		int i = CommandsDescriptor.getInstance().loadOptions(fMIndex, args);
		if(i<0) return;
		String fMIndexFile = args[i++];
		String readsFile = args[i++];
		fMIndex = FMIndex.loadFromBinaries(fMIndexFile);
		FastaSequencesHandler fastaSequencesHandler = new FastaSequencesHandler();
		QualifiedSequenceList sequences = fastaSequencesHandler.loadSequences(readsFile);
		for (int j = 0; j < sequences.size(); j++) 
		{
			QualifiedSequence actual = sequences.get(i);
			List<GenomicRegion> r = fMIndex.search(actual.getCharacters().toString());
		}
	}

}
