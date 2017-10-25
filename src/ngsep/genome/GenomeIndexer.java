package ngsep.genome;

import ngsep.main.CommandsDescriptor;
import ngsep.sequences.FMIndex;

public class GenomeIndexer 
{

	private int tallyDistance;
	public static void main(String[] args) throws Exception 
	{
		GenomeIndexer instance = new GenomeIndexer();
		int i = CommandsDescriptor.getInstance().loadOptions(instance, args);
		String genomeFile=args[i++];
		String outPutFile=args[i++];
		instance.createIndex (genomeFile,outPutFile);
	}

	public void createIndex(String genomeFile, String outPutFile) throws Exception {
		// TODO Auto-generated method stub
		FMIndex fMIndex=FMIndex.create(genomeFile);
		fMIndex.save(outPutFile);
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
