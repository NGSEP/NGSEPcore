package ngsep.discovery;

import java.util.logging.Logger;

import ngsep.genome.ReferenceGenome;
import ngsep.main.CommandsDescriptor;
import ngsep.main.ProgressNotifier;
import ngsep.simulation.TillingPopulationSimulator;

public class TillingPoolsIndividualGenotyper {

	private Logger log = Logger.getLogger(TillingPopulationSimulator.class.getName());
	private ProgressNotifier progressNotifier=null;
	
	
	public static void main(String[] args) throws Exception {
		TillingPoolsIndividualGenotyper instance = new TillingPoolsIndividualGenotyper();
		int i = CommandsDescriptor.getInstance().loadOptions(instance, args);
		String poolsVCFFile = args[i++];
		String poolsDescriptor = args[i++];
		String outFilename = args[i++];
		instance.loadPools (poolsDescriptor);
		instance.run(poolsVCFFile, outFilename);

	}

	public void loadPools(String poolsDescriptor) {
		// TODO Auto-generated method stub
		
	}

	public void run(String poolsVCFFile, String outFilename) {
		// TODO Auto-generated method stub
		
	}

}
