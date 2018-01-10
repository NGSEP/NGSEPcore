package ngsep.sequences;

import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

import ngsep.alignments.ReadAlignment;
import ngsep.genome.ReferenceGenome;

public class FMIndex implements Serializable 
{
	/**
	 * Serial number
	 */
	private static final long serialVersionUID = 5577026857894649939L;
	private List<FMIndexSingleSequence> internalIndexes = new ArrayList<>();
	private List<List<SequenceMetadata>> realSequencesMap = new ArrayList<>();

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
			FMIndexSingleSequence seqForward = new FMIndexSingleSequence(seqChars.toString().toUpperCase());
			internalIndexes.add(seqForward);
			List<SequenceMetadata> metadata = new ArrayList<>();
			metadata.add(new SequenceMetadata(q.getName(), q.getLength(), false));
			realSequencesMap.add(metadata);
			seqChars = seqChars.getReverseComplement();
			FMIndexSingleSequence seqReverse = new FMIndexSingleSequence(seqChars.toString().toUpperCase());
			internalIndexes.add(seqReverse);
			metadata = new ArrayList<>();
			metadata.add(new SequenceMetadata(q.getName(), q.getLength(), true));
			realSequencesMap.add(metadata);
		}
	}
	
	public void loadUnnamedSequences (String groupName, List<? extends CharSequence> sequences) {
		int n = sequences.size();
		StringBuffer internalSequence = new StringBuffer();
		List<SequenceMetadata> internalSequenceMetadata = new ArrayList<>();
		for(int i=0;i<n;i++) {
			String next = sequences.get(i).toString();
			if(internalSequence.length() + next.length() > 1000000000) {
				System.out.println("Building index for "+internalSequenceMetadata.size()+" sequences. Total sequence length: "+internalSequence.length());
				long time = System.currentTimeMillis();
				FMIndexSingleSequence index = new FMIndexSingleSequence(internalSequence);
				System.out.println("Built index in "+(System.currentTimeMillis()-time)+" milliseconds");
				internalIndexes.add(index);
				realSequencesMap.add(internalSequenceMetadata);
				internalSequence = new StringBuffer();
				internalSequenceMetadata = new ArrayList<>();
			}
			internalSequence.append(next);
			internalSequenceMetadata.add(new SequenceMetadata(groupName+"_"+i, next.length(), false));
			
		}
		System.out.println("Building index for "+internalSequenceMetadata.size()+" sequences. Total sequence length: "+internalSequence.length());
		long time = System.currentTimeMillis();
		FMIndexSingleSequence index = new FMIndexSingleSequence(internalSequence);
		System.out.println("Built index in "+(System.currentTimeMillis()-time)+" milliseconds");
		internalIndexes.add(index);
		realSequencesMap.add(internalSequenceMetadata);
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
	
	public List<ReadAlignment> search (String searchSequence) {
		List<ReadAlignment> alignments = new ArrayList<>();
		for (int i=0;i<internalIndexes.size();i++) {
			FMIndexSingleSequence idxSeq = internalIndexes.get(i);
			List<Integer> matches = idxSeq.search(searchSequence);
			for (int internalPosMatch:matches) {
				ReadAlignment alignment = buildAlignmentFromMetadata(searchSequence, i,internalPosMatch);
				if(alignment!=null) alignments.add(alignment);
			}
		}
		return alignments;
	}
	private ReadAlignment buildAlignmentFromMetadata(String searchSequence, int i, int internalPosMatch) {
		List<SequenceMetadata> metadata = realSequencesMap.get(i);
		int internalPostStartSeq = 0;
		for(SequenceMetadata seqM:metadata) {
			if(internalPostStartSeq + seqM.getLength() > internalPosMatch) {
				String seqName = seqM.getSeqName();
				int first = internalPosMatch-internalPostStartSeq;
				int last = first + searchSequence.length() - 1;
				int flags = 0;
				if(seqM.isNegativeStrand()) {
					last = seqM.getLength()-1-first;
					first = last - searchSequence.length() + 1 ;
					flags+=ReadAlignment.FLAG_READ_REVERSE_STRAND;
				}
				if(first < 0) return null;
				if(last>=seqM.getLength()) return null;
				return new ReadAlignment(seqName, first, last,searchSequence.length(),flags);
			}
			internalPostStartSeq += seqM.getLength();
		}
		return null;
	}

	/*public List<GenomicRegion> search (String sequenceName, String searchSequence) 
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
	}*/

	
}
class SequenceMetadata {
	private String seqName;
	private int length;
	private boolean negativeStrand;
	public SequenceMetadata(String seqName, int length, boolean negativeStrand) {
		super();
		this.seqName = seqName;
		this.length = length;
		this.negativeStrand = negativeStrand;
	}
	/**
	 * @return the seqName
	 */
	public String getSeqName() {
		return seqName;
	}
	/**
	 * @return the length
	 */
	public int getLength() {
		return length;
	}
	/**
	 * @return the negativeStrand
	 */
	public boolean isNegativeStrand() {
		return negativeStrand;
	}
	
	
	
}
