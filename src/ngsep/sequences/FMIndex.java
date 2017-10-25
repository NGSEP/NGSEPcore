package ngsep.sequences;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

import ngsep.genome.GenomicRegion;
import ngsep.genome.ReferenceGenome;

public class FMIndex implements Serializable 
{
	private List<FMIndexSingleSequence> singleSequenceIndexes = new ArrayList<>();

	/**
	 * 
	 * @param filename
	 * @throws Exception 
	 */
	public static FMIndex create(String filename) throws Exception 
	{
		FMIndex fmIndex = new FMIndex();
		ReferenceGenome referenceGenome = new ReferenceGenome(filename);
		int n = referenceGenome.getNumSequences();
		for (int i = 0; i < n; i++) 
		{
			QualifiedSequence q = referenceGenome.getSequenceByIndex(i);
			FMIndexSingleSequence actual = new FMIndexSingleSequence(q);
			fmIndex.singleSequenceIndexes.add(actual);
		}
		return fmIndex;
	}
	
	public static FMIndex loadFromBinaries(String filename) throws Exception
	{
		FMIndex fmIndex = new FMIndex();
		ObjectInputStream ois = new ObjectInputStream(new FileInputStream(filename));
		fmIndex.singleSequenceIndexes= (List<FMIndexSingleSequence>) ois.readObject();
		return fmIndex;
	}
	
	public void save (String filename) throws Exception 
	{
		ObjectOutputStream oos = new ObjectOutputStream( new FileOutputStream( filename+".fmIndex" ) );
		oos.writeObject( singleSequenceIndexes );
		oos.close( ); 
	}
	
	public List<GenomicRegion> search (String searchSequence) 
	{
		List<GenomicRegion> alignments = new ArrayList<>();
		for (int i = 0; i < singleSequenceIndexes.size(); i++) 
		{
			List<GenomicRegion> actual = singleSequenceIndexes.get(i).search(searchSequence);
			alignments.addAll(actual);
		}
		return alignments;
	}
	public List<GenomicRegion> search (String sequenceName, String searchSequence) 
	{
		List<GenomicRegion> alignments = new ArrayList<>();
		for (int i = 0; i < singleSequenceIndexes.size(); i++) 
		{
			FMIndexSingleSequence actual = singleSequenceIndexes.get(i);
			if(actual.getSequenceName().equals(searchSequence))
			{
				return actual.search(searchSequence);
			}
		}
		return alignments;
	}
	
}
