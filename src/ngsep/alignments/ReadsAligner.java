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
			for (int k = 0; k < r.size(); k++) 
			{
				System.out.println(
						//1.query name
						read.getName()+"\t"+
								
						//2.Flag
						"1\t"+
						
						//3.reference sequence name
						r.get(k).getSequenceName()+"\t"+
						
						//4.POS
						r.get(k).getFirst()+"\t"+
						
						//5.MAPQ
						"255\t"+
						
						//6.CIGAR
						"*\t"+
						
						//7. RNEXT
						"*\t"+
						
						//8. PNEXT
						"0\t"+
						
						//9. TLEN
						"0\t"+
						
						//10. SEQ
						"*\t"+
						
						//11. QUAL
						"*\t"
						
						);
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
