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
package ngsep.variants.imputation;

import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.logging.Logger;

import ngsep.main.CommandsDescriptor;
import ngsep.main.ProgressNotifier;
import ngsep.variants.CalledGenomicVariant;
import ngsep.variants.CalledSNV;
import ngsep.variants.GenomicVariant;
import ngsep.variants.SNV;
import ngsep.variants.Sample;
import ngsep.vcf.VCFFileHeader;
import ngsep.vcf.VCFFileReader;
import ngsep.vcf.VCFFileWriter;
import ngsep.vcf.VCFRecord;


public class GenotypeImputer {
	public static final int DEF_K = 20;
	private Logger log = Logger.getLogger(GenotypeImputer.class.getName());
	private ProgressNotifier progressNotifier=null;
	private int progress = 0;
	private Double avgCMPerKbp = null;
	private boolean fixedTransitions = false;
	private List<String> parentIds = new ArrayList<String>() ;
	private int k = DEF_K;
	//private PrintStream outFounders;
	private PrintStream outAssignments;
	private PrintStream outGenotypes;
	 
	
	/**
	 * @param args
	 */
	public static void main(String[] args) throws Exception {
		if (args.length == 0 || args[0].equals("-h") || args[0].equals("--help")){
			CommandsDescriptor.getInstance().printHelp(GenotypeImputer.class);
			return;
		}
		GenotypeImputer instance = new GenotypeImputer();
		int i=0;
		while(i<args.length && args[i].charAt(0)=='-') {
			if("-c".equals(args[i])) {
				i++;
				instance.setAvgCMPerKbp(Double.parseDouble(args[i]));
			} else if("-p".equals(args[i])) {
				i++;
				instance.setParentIds(Arrays.asList(args[i].split(",")));
			} else if("-k".equals(args[i])) {
				i++;
				instance.setK(Integer.parseInt(args[i]));
			} else if("-t".equals(args[i])) {
				instance.setFixedTransitions(true);
			} else {
				System.err.println("Unrecognized option: "+args[i]);
				CommandsDescriptor.getInstance().printHelp(GenotypeImputer.class);
				return;
			}
			i++;
		}
		String vcfFile = args[i++];
		String outPrefix = args[i++];
		//String fileFounders = outPrefix+"_founders.txt";
		String fileAssignments = outPrefix+"_assignments.txt";
		String fileGenotypes = outPrefix+"_imputed.vcf";
		
		//instance.outFounders = null;
		instance.outAssignments = null;
		instance.outGenotypes = null;
		try {
			//instance.outFounders = new PrintStream(fileFounders);
			instance.outAssignments = new PrintStream(fileAssignments);
			instance.outGenotypes = new PrintStream(fileGenotypes);
			instance.impute(vcfFile);
		} finally {
			//if(instance.outFounders!=null) instance.outFounders.close();
			if(instance.outAssignments!=null) instance.outAssignments.close();
		}

	}
	
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

	public double getAvgCMPerKbp() {
		return avgCMPerKbp;
	}

	public void setAvgCMPerKbp(double avgCMPerKbp) {
		this.avgCMPerKbp = avgCMPerKbp;
	}
	
	public List<String> getParentIds() {
		return parentIds;
	}
	
	public void setParentIds(List<String> parentIds) {
		this.parentIds = parentIds;
	}
	
	public int getK() {
		return k;
	}

	public void setK(int k) {
		this.k = k;
	}
	
	public boolean isFixedTransitions() {
		return fixedTransitions;
	}

	public void setFixedTransitions(boolean fixedTransitions) {
		this.fixedTransitions = fixedTransitions;
	}

	public PrintStream getOutAssignments() {
		return outAssignments;
	}



	public void setOutAssignments(PrintStream outAssignments) {
		this.outAssignments = outAssignments;
	}



	public PrintStream getOutGenotypes() {
		return outGenotypes;
	}



	public void setOutGenotypes(PrintStream outGenotypes) {
		this.outGenotypes = outGenotypes;
	}



	public void impute(String filename) throws IOException {
		logParameters(filename);
		Map<String, List<CalledSNV>> genotypes = new TreeMap<String, List<CalledSNV>>();
		List<VCFRecord> records = new ArrayList<VCFRecord>();
		VCFFileReader reader = null;
		VCFFileWriter out = new VCFFileWriter();
		try {
			reader = new VCFFileReader(filename);
			if(log!=null) reader.setLog(log);
			reader.setLoadMode(VCFFileReader.LOAD_MODE_QUALITY);
			//Ids of the samples
			VCFFileHeader header = reader.getHeader();
			out.printHeader(header,outGenotypes);
			List<Sample> samples = header.getSamples();
			for(Sample sample:samples) genotypes.put(sample.getId(), new ArrayList<CalledSNV>());
			String lastSeqName = null;
			Iterator<VCFRecord> it = reader.iterator();
			while(it.hasNext()) {
				VCFRecord record = it.next();
				GenomicVariant var = record.getVariant();
				if(!(var instanceof SNV)) continue;
				
				SNV snv = (SNV) var;
				if(!var.getSequenceName().equals(lastSeqName)) {
					if(lastSeqName!=null) {
						imputeGenotypes(genotypes);
						for(VCFRecord r:records) out.printVCFRecord(r, outGenotypes);
						progress++;
						if(progressNotifier!=null && !progressNotifier.keepRunning(progress)) return;
					}
					lastSeqName = var.getSequenceName();
					for(List<CalledSNV> genotypesSample:genotypes.values()) genotypesSample.clear();
					records.clear();
				}
				records.add(record);
				List<CalledGenomicVariant> genotypeCalls = record.getCalls();
				for(int i=0;i<genotypeCalls.size();i++) {
					String sampleId = samples.get(i).getId();
					CalledGenomicVariant genotypeCall = genotypeCalls.get(i);
					CalledSNV csnv;
					if(genotypeCall instanceof CalledSNV) csnv = (CalledSNV)genotypeCall;
					else csnv = new CalledSNV(snv, CalledSNV.GENOTYPE_UNDECIDED);
					//LightweightCalledSNV lcsnv = new LightweightCalledSNV(csnv);
					genotypes.get(sampleId).add(csnv);
				}
				
			}
			if(lastSeqName!=null) {
				imputeGenotypes(genotypes);
				for(VCFRecord r:records) out.printVCFRecord(r, outGenotypes);
			}
		} finally {
			if (reader!=null) reader.close();
		}
	}
	
	private void logParameters(String filename) {
		log.info("Running imputation for VCF file: "+filename);
		if(avgCMPerKbp!=null) log.info("Average cm per kbp: "+avgCMPerKbp);
		if(parentIds!=null) log.info("Number of parents: "+parentIds.size());
	}

	public void imputeGenotypes(Map<String,List<CalledSNV>> genotypes) {
		progress = 0;
		imputeGenotypesHMM(genotypes);
	}

	public void imputeGenotypesHMM(Map<String, List<CalledSNV>> genotypes) {
		
		GenotypeImputationHMM  hmm = GenotypeImputationHMM.createHMM(genotypes, parentIds, k);
		if(avgCMPerKbp!=null) hmm.setAvgCMPerKbp(avgCMPerKbp);
		hmm.setFixedTransitions(fixedTransitions);
		hmm.setLog(log);
		
		if(progressNotifier!=null) {
			progress++;
			if(!progressNotifier.keepRunning(progress)) return;
		}
		List<CalledSNV> snvs = genotypes.values().iterator().next();
		Map<String, List<Integer>> assignments = hmm.imputeGenotypes(genotypes); 
		printAssignments(snvs,assignments,hmm);
	}
	
	private void printAssignments(List<CalledSNV> snvs, Map<String, List<Integer>> assignments,GenotypeImputationHMM hmm) {
		Set<String> sampleIds = assignments.keySet();
		int m = snvs.size();
		outAssignments.print("Chr\tPos");
		for(String sampleId:sampleIds) outAssignments.print("\t"+sampleId);
		outAssignments.println();
		for (int i=0;i<m;i++) {
			CalledSNV snv = snvs.get(i);
			outAssignments.print(snv.getSequenceName()+"\t"+snv.getPosition());
			for(String sampleId: sampleIds) {
				int assignment = assignments.get(sampleId).get(i);
				String stateId = hmm.getState(assignment).getId();
				if(stateId==null) outAssignments.print("\t"+assignment);
				else outAssignments.print("\t"+stateId);
			}
			outAssignments.println();
		}
		
	}
	
}
