package ngsep.alignments;

import java.util.List;

import ngsep.genome.GenomeIndexer;
import ngsep.genome.GenomicRegion;
import ngsep.main.CommandsDescriptor;
import ngsep.sequences.FMIndex;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;
import ngsep.sequences.io.FastaSequencesHandler;

public class ReadsAligner {

	public static void main(String[] args) throws Exception 
	{
		ReadsAligner instance = new ReadsAligner();
		int i = CommandsDescriptor.getInstance().loadOptions(instance, args);
		if(i<0) return;
		String fMIndexFile = args[i++];
		String readsFile = args[i++];
		
		FMIndex fMIndex = instance.loadIndex(fMIndexFile);
		FastaSequencesHandler fastaSequencesHandler = new FastaSequencesHandler();
		QualifiedSequenceList sequences = fastaSequencesHandler.loadSequences(readsFile);
		
		
		
		for (int j = 0; j < sequences.size(); j++) 
		{
			QualifiedSequence actual = sequences.get(j);
			List<GenomicRegion> r = fMIndex.search(actual.getCharacters().toString());
			System.out.println(actual.getName() +" found in:");
			for (int k = 0; k < r.size(); k++) 
			{
				System.out.println(r.get(k).getSequenceName()+" "+r.get(k).getFirst()+" "+r.get(k).getLast());
			}
		}
	}
	
	public FMIndex loadIndex(String fMIndexFile) throws Exception 
	{
		FMIndex f = FMIndex.loadFromBinaries(fMIndexFile);
		return f;
	}

}
