package ngsep.alignments;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.List;

import ngsep.genome.GenomicRegion;
import ngsep.main.CommandsDescriptor;
import ngsep.sequences.FMIndex;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;
import ngsep.sequences.RawRead;
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
		
		FileInputStream fis = new FileInputStream(readsFile);
		BufferedReader in = new BufferedReader(new InputStreamReader(fis));
		RawRead read = RawRead.load(in);
		while(read!=null) 
		{
			List<GenomicRegion> r = fMIndex.search(read.getCharacters().toString());
			System.out.println(read.getName() +" found in:");
			for (int k = 0; k < r.size(); k++) 
			{
				System.out.println(r.get(k).getSequenceName()+" "+r.get(k).getFirst()+" "+r.get(k).getLast());
			}
			read = RawRead.load(in);
		}
		
		fis.close();
		
			
	}
	
	public FMIndex loadIndex(String fMIndexFile) throws Exception 
	{
		FMIndex f = FMIndex.loadFromBinaries(fMIndexFile);
		return f;
	}

}
