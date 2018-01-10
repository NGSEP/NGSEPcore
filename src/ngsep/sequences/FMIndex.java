package ngsep.sequences;

import java.io.FileInputStream;
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
	/**
	 * Serial number
	 */
	private static final long serialVersionUID = 5577026857894649939L;
	private List<FMIndexSingleSequence> singleSequenceIndexes = new ArrayList<>();

	/**
	 * Loads the given genome in this FMIndex
	 * @param filename with the genome to load
	 * @throws IOException If the genome file can not be read
	 */
	public void loadGenome(String genomeFilename) throws IOException 
	{
		ReferenceGenome referenceGenome = new ReferenceGenome(genomeFilename);
		int n = referenceGenome.getNumSequences();
		
		for (int i = 0; i < n; i++) 
		{
			QualifiedSequence q = referenceGenome.getSequenceByIndex(i);
			DNAMaskedSequence seqChars = (DNAMaskedSequence)q.getCharacters();
			DNAMaskedSequence reverseComplement = seqChars.getReverseComplement();
			FMIndexSingleSequence seqForward = new FMIndexSingleSequence(q.getName()+"_F",seqChars);
			singleSequenceIndexes.add(seqForward);
			FMIndexSingleSequence seqReverse = new FMIndexSingleSequence(q.getName()+"_R",reverseComplement);
			singleSequenceIndexes.add(seqReverse);
		}
	}
	
	public void loadUnnamedSequences (String groupName, List<? extends CharSequence> sequences) {
		int n = sequences.size();
		for(int i=0;i<n;i++) {
			FMIndexSingleSequence seqForward = new FMIndexSingleSequence(groupName+"_"+i,sequences.get(i));
			singleSequenceIndexes.add(seqForward);
			if((i+1)%50 == 0) System.out.println("Built index for "+(i+1)+" sequences");
		}
	}
	
	/**
	 * Loads an instance of the FMIndex from a serialized binary file
	 * @param filename Binary file with the serialization of an FMIndex
	 * @return FMIndex serialized in the given file
	 * @throws IOException If there were errors reading the file
	 */
	public static FMIndex loadFromBinaries(String filename) throws IOException
	{
		FMIndex fmIndex = new FMIndex();
		try (FileInputStream fis = new FileInputStream(filename);
			 ObjectInputStream ois = new ObjectInputStream(fis);) {
			fmIndex = (FMIndex) ois.readObject();
		} catch (ClassNotFoundException e) {
			throw new RuntimeException("FMIndex class not found",e);
		}
		
		return fmIndex;
	}
	/**
	 * Saves this FM-Index as a serializable object
	 * @param filename
	 * @throws IOException
	 */
	public void save (String filename) throws IOException 
	{
		try (FileOutputStream fos = new FileOutputStream( filename );
			 ObjectOutputStream oos = new ObjectOutputStream( fos );) {
			oos.writeObject( this );
		}
	}
	
	public List<GenomicRegion> search (String searchSequence) {
		List<GenomicRegion> alignments = new ArrayList<>();
		for (FMIndexSingleSequence idxSeq:singleSequenceIndexes) {
			List<GenomicRegion> actual = idxSeq.search(searchSequence);
			alignments.addAll(actual);
		}
		return alignments;
	}
	public List<GenomicRegion> search (String sequenceName, String searchSequence) 
	{
		List<GenomicRegion> alignments = new ArrayList<>();
		for (FMIndexSingleSequence idxSeq:singleSequenceIndexes) 
		{
			if(idxSeq.getSequenceName().equals(searchSequence))
			{
				return idxSeq.search(searchSequence);
			}
		}
		return alignments;
	}

	
}
