package ngsep.vcf;

import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.util.Iterator;
import java.util.logging.Logger;

import ngsep.genome.ReferenceGenome;
import ngsep.main.CommandsDescriptor;
import ngsep.main.ProgressNotifier;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;

public class VCFVariantDensityCalculator {

	public static final int DEF_WINDOW_LENGTH=100000;
	private Logger log = Logger.getLogger(VCFVariantDensityCalculator.class.getName());
	private ProgressNotifier progressNotifier=null;
	private ReferenceGenome genome;
	private int windowLength = DEF_WINDOW_LENGTH;
	
	public static void main(String[] args) throws Exception {
		VCFVariantDensityCalculator instance = new VCFVariantDensityCalculator();
		int i=CommandsDescriptor.getInstance().loadOptions(instance, args);
		instance.genome = new ReferenceGenome(args[i++]);
		boolean systemInput = "-".equals(args[i]);
		if(systemInput) {
			instance.run(System.in, System.out);
		} else {
			String filename = args[i];
			instance.run(filename, System.out);
		}

	}
	
	/**
	 * @return the windowLength
	 */
	public int getWindowLength() {
		return windowLength;
	}



	/**
	 * @param windowLength the windowLength to set
	 */
	public void setWindowLength(int windowLength) {
		this.windowLength = windowLength;
	}



	public void run(String filename, PrintStream out) throws IOException {
		try (VCFFileReader in = new VCFFileReader(filename)) { 
			run(in, out);
		}
	}
	public void run(InputStream fis, PrintStream out) throws IOException {
		try (VCFFileReader in = new VCFFileReader(fis)) {
			run(in, out);
		}	
	}
	public void run(VCFFileReader in, PrintStream out) throws IOException {
		if(log!=null)in.setLog(log);
		in.setLoadMode(VCFFileReader.LOAD_MODE_MINIMAL);
		QualifiedSequenceList sequences = genome.getSequencesMetadata();
		in.setSequences(sequences);
		Iterator<VCFRecord> it = in.iterator();
		int n=0;
		
		String seqName = null;
		int endWindow = 0;
		int count = 0;
		int idxSeqNames = 0;
		while(it.hasNext()) {
			VCFRecord record = it.next();
			if(!record.getSequenceName().equals(seqName)) {
				if(seqName!=null) out.println(seqName+"\t"+(endWindow-windowLength+1)+"\t"+endWindow+"\t"+count);
				QualifiedSequence seq = sequences.get(idxSeqNames); 
				while (!seq.getName().equals(record.getSequenceName())) {
					processZeroes (seq,endWindow, out);
					idxSeqNames++;
					endWindow = 0;
					seq = sequences.get(idxSeqNames);
				}
				seqName = record.getSequenceName();
				endWindow = windowLength;
				count = 0;
			} else if (record.getFirst()>endWindow) {
				if(seqName!=null) out.println(seqName+"\t"+(endWindow-windowLength+1)+"\t"+endWindow+"\t"+count);
				endWindow+=windowLength;
				count = 0;
			}
			for(;endWindow<record.getFirst();endWindow+=windowLength) {
				out.println(seqName+"\t"+(endWindow-windowLength+1)+"\t"+endWindow+"\t0");
			}
			count++;
			n++;
			if (progressNotifier!=null && n%1000==0) {
				int progress = n/1000;
				if (!progressNotifier.keepRunning(progress)) {
					out.flush();
					return;
				}
			}
		}
		
	}

	private void processZeroes(QualifiedSequence seq, int endWindow, PrintStream out) {
		while(endWindow<seq.getLength()) {
			endWindow+=windowLength;
			out.println(seq.getName()+"\t"+(endWindow-windowLength+1)+"\t"+endWindow+"\t0");
			
		}
		
	}

}
