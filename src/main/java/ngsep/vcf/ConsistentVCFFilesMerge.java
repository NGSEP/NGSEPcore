/*******************************************************************************
 * NGSEP - Next Generation Sequencing Experience Platform
 * Copyright 2016 Jorge Duitama
 *
 * This file is part of NGSEP.
 *
 *     NGSEP is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     NGSEP is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with NGSEP.  If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/
package ngsep.vcf;

import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.logging.Logger;

import ngsep.genome.GenomicRegionComparator;
import ngsep.main.CommandsDescriptor;
import ngsep.main.ProgressNotifier;
import ngsep.sequences.QualifiedSequenceList;
import ngsep.sequences.io.SimpleSequenceListLoader;
import ngsep.variants.CalledGenomicVariant;
import ngsep.variants.CalledGenomicVariantImpl;
import ngsep.variants.DiversityStatistics;
import ngsep.variants.GenomicVariant;
import ngsep.variants.GenomicVariantAnnotation;
import ngsep.variants.GenomicVariantImpl;
import ngsep.variants.Sample;
/**
 * @author Jorge Duitama
 */
public class ConsistentVCFFilesMerge {
	// Constants for default values
	
	// Logging and progress
	private Logger log = Logger.getLogger(IndividualSampleVariantsMerge.class.getName());
	private ProgressNotifier progressNotifier=null;
	
	// Parameters
	private String sequenceNamesFile = null;
	private String outputFile = null;
	private List<String> vcfFiles = new ArrayList<String>();
	
	// Get and set methods
	public Logger getLog() {
		return log;
	}
	public void setLog(Logger log) {
		this.log = log;
	}
	
	public ProgressNotifier getProgressNotifier() {
		return progressNotifier;
	}
	public void setProgressNotifier(ProgressNotifier progressNotifier) { 
		this.progressNotifier = progressNotifier;
	}
	
	
	public String getSequenceNamesFile() {
		return sequenceNamesFile;
	}
	public void setSequenceNamesFile(String sequenceNamesFile) {
		this.sequenceNamesFile = sequenceNamesFile;
	}
	
	public String getOutputFile() {
		return outputFile;
	}
	public void setOutputFile(String outputFile) {
		this.outputFile = outputFile;
	}

	public static void main(String[] args) throws Exception {
		ConsistentVCFFilesMerge instance = new ConsistentVCFFilesMerge();
		int i = CommandsDescriptor.getInstance().loadOptions(instance, args);
		instance.vcfFiles = new ArrayList<String>();
		for(;i<args.length;i++) {
			instance.vcfFiles.add(args[i]);
		}
		instance.run();
	}
	
	public void run() throws IOException {
		if (sequenceNamesFile==null || sequenceNamesFile.isEmpty()) throw new IOException("The file with sequence names is a required parameter");
		if (outputFile==null || outputFile.isEmpty()) throw new IOException("The output file is a required parameter");
		SimpleSequenceListLoader listHandler = new SimpleSequenceListLoader();
		QualifiedSequenceList sequenceNames = listHandler.loadSequences(sequenceNamesFile);
		log.info("Loaded "+sequenceNames.size()+" sequence names from file "+sequenceNamesFile);
		mergeFiles(sequenceNames, vcfFiles, outputFile);
		log.info("Process finished");
	}
	public void mergeFiles(QualifiedSequenceList sequenceNames, List<String> vcfFiles, String outputFile) throws IOException {
		try (PrintStream out = new PrintStream(outputFile)) {
			mergeFiles(sequenceNames, vcfFiles, out);
		}
	}

	public void mergeFiles(QualifiedSequenceList sequenceNames, List<String> vcfFiles, PrintStream out) throws IOException {
		VCFFileWriter writer = new VCFFileWriter();
		List<VCFFileReader> readers = new ArrayList<VCFFileReader>();
		List<Iterator<VCFRecord>> iterators = new ArrayList<Iterator<VCFRecord>>();
		VCFFileHeader outHeader = VCFFileHeader.makeDefaultEmptyHeader();
		VCFRecord [] calls = new VCFRecord [vcfFiles.size()];
		int [] samplesPerFile = new int[vcfFiles.size()];
		for(int i=0;i<vcfFiles.size();i++) {
			VCFFileReader reader = new VCFFileReader(vcfFiles.get(i));
			reader.setLog(log);
			readers.add(reader);
			addFilesamples(reader,outHeader);
			samplesPerFile[i] = reader.getHeader().getSamples().size();
			Iterator<VCFRecord> it = reader.iterator();
			if(it.hasNext()) calls[i] = it.next();
			else calls[i] = null;
			iterators.add(it);
		}
		GenomicVariantStrictComparator comparator = new GenomicVariantStrictComparator(new GenomicRegionComparator(sequenceNames));
		comparator.setLog(log);
		writer.printHeader(outHeader, out);
		int n=0;
		while (true){
			VCFRecord record = createNextRecord(iterators,calls,samplesPerFile,comparator,outHeader);
			if(record == null) break;
			boolean save = false;
			for(CalledGenomicVariant cv:record.getCalls()) {
				if(!cv.isUndecided()) {
					save = true;
					break;
				}
			}
			if(save) {
				writer.printVCFRecord(record, out);
			}
			n++;
			if (progressNotifier!=null && n%1000==0) {
				int progress = n/1000;
				if (!progressNotifier.keepRunning(progress)) {
					for(VCFFileReader reader:readers) reader.close();
					out.flush();
					return;
				}
			}
		}
		for(VCFFileReader reader:readers) reader.close();
		
	}
	private void addFilesamples(VCFFileReader reader, VCFFileHeader outHeader) {
		VCFFileHeader inHeader = reader.getHeader();
		List<Sample> samples = inHeader.getSamples();
		Set<String> sampleIdsHeaderLine = inHeader.getSamplesWithHeaderLine().keySet();
		for(Sample s:samples) {
			outHeader.addSample(s, sampleIdsHeaderLine.contains(s.getId()));
		}
	}


	private VCFRecord createNextRecord(List<Iterator<VCFRecord>> iterators,VCFRecord[] currentCalls,int [] samplesPerFile,GenomicVariantStrictComparator comparator, VCFFileHeader outHeader) {
		GenomicVariant nextVariant = null;
		int [] nextFormat = null;
		for(int i=0;i<currentCalls.length;i++) {
			VCFRecord nextCall = currentCalls[i];
			if(nextCall!=null) {
				if(nextVariant==null || comparator.compare(nextCall.getVariant(), nextVariant)<0) {
					nextVariant = nextCall.getVariant();
					nextFormat = nextCall.getFieldsFormat();
				}
			}
		}
		if(nextVariant == null) return null;
		CalledGenomicVariantImpl undecidedCall = new CalledGenomicVariantImpl(nextVariant, new byte[0]);
		List<CalledGenomicVariant> calls = new ArrayList<CalledGenomicVariant>();
		short maxQS = nextVariant.getVariantQS();
		byte type = nextVariant.getType();
		int countCNVs = 0;
		for(int i=0;i<currentCalls.length;i++) {
			VCFRecord record = currentCalls[i];
			if(record == null) {
				for(int j=0;j<samplesPerFile[i];j++) calls.add(undecidedCall);
			} else {
				GenomicVariant var2 = record.getVariant();
				if(comparator.compare(nextVariant, var2)==0) {
					for(CalledGenomicVariant call:record.getCalls()) {
						calls.add(call);
						short genotypeQS = call.getGenotypeQuality();
						if(!call.isUndecided() && !call.isHomozygousReference() && genotypeQS > maxQS) {
							maxQS = genotypeQS;
						}
					}
					GenomicVariantAnnotation ann = record.getInfoField(GenomicVariantAnnotation.ATTRIBUTE_IN_CNV);
					//Support for old boolean format kept for compatibility with old CNV files
					if(ann!=null && ann.getValue()!=null) {
						if((ann.getValue() instanceof Boolean) && ((Boolean)ann.getValue()).booleanValue()) countCNVs++;
						else if ((ann.getValue() instanceof Integer)) countCNVs+= (Integer)ann.getValue();
						else if ((ann.getValue() instanceof String)) countCNVs+= Integer.parseInt((String)ann.getValue());
					}
					type = GenomicVariantImpl.mergeType(type, var2.getType());
					Iterator<VCFRecord> it = iterators.get(i);
					if(it.hasNext()) currentCalls[i] = it.next();
					else currentCalls[i] = null;
				} else {
					for(int j=0;j<samplesPerFile[i];j++) calls.add(undecidedCall);
				}
			}
		}
		nextVariant.setVariantQS(maxQS);
		nextVariant.setType(type);
		VCFRecord answer = new VCFRecord(nextVariant, new ArrayList<String>(), new ArrayList<GenomicVariantAnnotation>(), nextFormat, calls,outHeader);
		DiversityStatistics divStats = DiversityStatistics.calculateDiversityStatistics(calls, false);
		if(divStats!=null) {
			answer.addAnnotation(new GenomicVariantAnnotation(nextVariant, GenomicVariantAnnotation.ATTRIBUTE_SAMPLES_GENOTYPED, divStats.getNumSamplesGenotyped()));
			answer.addAnnotation(new GenomicVariantAnnotation(nextVariant, GenomicVariantAnnotation.ATTRIBUTE_NUMBER_ALLELES, divStats.getNumCalledAlleles()));
			if(divStats.getNumCalledAlleles()==2) answer.addAnnotation(new GenomicVariantAnnotation(nextVariant, GenomicVariantAnnotation.ATTRIBUTE_MAF, divStats.getMaf()));
		}
		if(countCNVs>0) answer.addAnnotation(new GenomicVariantAnnotation(nextVariant, GenomicVariantAnnotation.ATTRIBUTE_IN_CNV, countCNVs));
		return answer;
		
	}
}

class GenomicVariantStrictComparator implements Comparator<GenomicVariant> {

	private GenomicRegionComparator regComp;
	private Logger log;
	
	
	
	/**
	 * @param regComp
	 */
	public GenomicVariantStrictComparator(GenomicRegionComparator regComp) {
		super();
		this.regComp = regComp;
	}
	
	public Logger getLog() {
		return log;
	}

	public void setLog(Logger log) {
		this.log = log;
	}

	@Override
	public int compare(GenomicVariant v1, GenomicVariant v2) {
		int cmp = regComp.compare(v1, v2);
		if(cmp!=0) return cmp;
		String [] all1 = v1.getAlleles();
		String [] all2 = v2.getAlleles();
		//Compare reference alleles and choose first that spanning the smallest region
		cmp = all1[0].compareTo(all2[0]);
		if(cmp!=0) {
			log.warning("Inconsistent variants found at "+v1.getSequenceName()+":"+v1.getFirst()+". References differ. Reference 1: "+all1[0]+" Reference 2: "+all2[0]);
			return cmp;
		}
		//Choose first the variant with less alleles
		cmp = all1.length - all2.length;
		if(cmp!=0) {
			log.warning("Inconsistent variants found at "+v1.getSequenceName()+":"+v1.getFirst()+". Number of alleles differ. Var 1 has "+all1.length+" alleles. Var 2 has "+all2.length+" alleles");
			return cmp;
		}
		//Compare alternative alleles
		for(int i=1;i<all1.length;i++) {
			cmp = all1[i].compareTo(all2[i]);
			if(cmp!=0) {
				log.warning("Inconsistent variants found at "+v1.getSequenceName()+":"+v1.getFirst()+". Alternative alleles differ. Inconsistent alternative alleles are: "+all1[i]+" and "+all2[i]);
				return cmp;
			}
		}
		return 0;
	}	
}