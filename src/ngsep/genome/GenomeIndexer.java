package ngsep.genome;

import java.io.IOException;

import ngsep.main.CommandsDescriptor;
import ngsep.sequences.FMIndex;

public class GenomeIndexer 
{

	private int tallyDistance;
	public static void main(String[] args) throws Exception 
	{
		long t = System.currentTimeMillis();
		GenomeIndexer instance = new GenomeIndexer();
		int i = CommandsDescriptor.getInstance().loadOptions(instance, args);
		String genomeFile=args[i++];
		String outputFile=args[i++];
		instance.createIndex (genomeFile,outputFile);
		System.out.println(System.currentTimeMillis()-t);
	}

	public void createIndex(String genomeFile, String outputFile) throws IOException {
		// TODO Auto-generated method stub
		FMIndex fMIndex= new FMIndex();
		fMIndex.loadGenome(genomeFile);
		fMIndex.save(outputFile);
	}

	public int getTallyDistance() {
		return tallyDistance;
	}

	public void setTallyDistance(int tallyDistance) {
		this.tallyDistance = tallyDistance;
	}
	public void setTallyDistance(Integer tallyDistance) {
		this.setTallyDistance(tallyDistance.intValue());
	}
	

}
