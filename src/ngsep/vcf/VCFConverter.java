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

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.logging.Logger;

import ngsep.main.CommandsDescriptor;
import ngsep.main.ProgressNotifier;
import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.io.FastaSequencesHandler;
import ngsep.variants.CalledGenomicVariant;
import ngsep.variants.GenomicVariant;
import ngsep.variants.SNV;
import ngsep.variants.Sample;
import ngsep.variants.io.SimpleSamplesFileHandler;

public class VCFConverter {
	
	// Constants for default values
	
	// Logging and progress
	private Logger log = Logger.getLogger(VCFConverter.class.getName());
	private ProgressNotifier progressNotifier=null;
	
	// Parameters
	private String inputFile = null;
	private String outputPrefix = null;
	private boolean printStructure = false;
	private boolean printFasta = false;
	private boolean printrrBLUP = false;
	private boolean printMatrix = false;
	private boolean printHapmap = false;
	private boolean printGWASPoly = false;
	private boolean printSpagedi = false;
	private boolean printPlink = false;
	private boolean printHaploview = false;
	private boolean printEmma = false;
	private boolean printPowerMarker = false;
	private boolean printEigensoft = false;
	private boolean printFlapjack = false;
	private boolean printDarwin = false;
	private boolean printTreeMix= false;
	private boolean printJoinMap= false;
	private boolean printPhase = false;
	private boolean printFineStructure = false;
	private boolean printGenePop = false;
	private String sequenceName = null;
	private String idParent1 = null;
	private String idParent2 = null;
	private String populationFile=null;
	
	
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
	
	public String getInputFile() {
		return inputFile;
	}
	public void setInputFile(String inputFile) {
		this.inputFile = inputFile;
	}
	
	public String getOutputPrefix() {
		return outputPrefix;
	}
	public void setOutputPrefix(String outputPrefix) {
		this.outputPrefix = outputPrefix;
	}
	
	public boolean isPrintStructure() {
		return printStructure;
	}
	public void setPrintStructure(boolean printStructure) {
		this.printStructure = printStructure;
	}
	public void setPrintStructure(Boolean printStructure) {
		this.setPrintStructure(printStructure.booleanValue());
	}
	
	public boolean isPrintFasta() {
		return printFasta;
	}
	public void setPrintFasta(boolean printFasta) {
		this.printFasta = printFasta;
	}
	public void setPrintFasta(Boolean printFasta) {
		this.setPrintFasta(printFasta.booleanValue());
	}

	public boolean isPrintrrBLUP() {
		return printrrBLUP;
	}
	public void setPrintrrBLUP(boolean printrrBLUP) {
		this.printrrBLUP = printrrBLUP;
	}
	public void setPrintrrBLUP(Boolean printrrBLUP) {
		this.setPrintrrBLUP(printrrBLUP.booleanValue());
	}

	public boolean isPrintMatrix() {
		return printMatrix;
	}
	public void setPrintMatrix(boolean printMatrix) {
		this.printMatrix = printMatrix;
	}
	public void setPrintMatrix(Boolean printMatrix) {
		this.setPrintMatrix(printMatrix.booleanValue());
	}

	public boolean isPrintHapmap() {
		return printHapmap;
	}
	public void setPrintHapmap(boolean printHapmap) {
		this.printHapmap = printHapmap;
	}
	public void setPrintHapmap(Boolean printHapmap) {
		this.setPrintHapmap(printHapmap.booleanValue());
	}
	
	public boolean isPrintGWASPoly() {
		return printGWASPoly;
	}
	public void setPrintGWASPoly(boolean printGWASPoly) {
		this.printGWASPoly = printGWASPoly;
	}
	public void setPrintGWASPoly(Boolean printGWASPoly) {
		this.setPrintGWASPoly(printGWASPoly.booleanValue());
	}
	
	public boolean isPrintSpagedi() {
		return printSpagedi;
	}
	public void setPrintSpagedi(boolean printSpagedi) {
		this.printSpagedi = printSpagedi;
	}
	public void setPrintSpagedi(Boolean printSpagedi) {
		this.setPrintSpagedi(printSpagedi.booleanValue());
	}
	public boolean isPrintPlink() {
		return printPlink;
	}

	public void setPrintPlink(boolean printPlink) {
		this.printPlink = printPlink;
	}	
	public void setPrintPlink(Boolean printPlink) {
		this.setPrintPlink(printPlink.booleanValue());
	}
	public boolean isPrintHaploview() {
		return printHaploview;
	}
	public void setPrintHaploview(boolean printHaploview) {
		this.printHaploview = printHaploview;
	}
	public void setPrintHaploview(Boolean printHaploview) {
		this.setPrintHaploview(printHaploview.booleanValue());
	}

	public boolean isPrintEmma() {
		return printEmma;
	}
	public void setPrintEmma(boolean printEmma) {
		this.printEmma = printEmma;
	}
	public void setPrintEmma(Boolean printEmma) {
		this.setPrintEmma(printEmma.booleanValue());
	}

	public boolean isPrintPowerMarker() {
		return printPowerMarker;
	}
	public void setPrintPowerMarker(boolean printPowerMarker) {
		this.printPowerMarker = printPowerMarker;
	}
	public void setPrintPowerMarker(Boolean printPowerMarker) {
		this.setPrintPowerMarker(printPowerMarker.booleanValue());
	}
	
	public boolean isPrintEigensoft() {
		return printEigensoft;
	}
	public void setPrintEigensoft(boolean printEigensoft) {
		this.printEigensoft = printEigensoft;
	}
	public void setPrintEigensoft(Boolean printEigensoft) {
		this.setPrintEigensoft(printEigensoft.booleanValue());
	}

	public boolean isPrintFlapjack() {
		return printFlapjack;
	}
	public void setPrintFlapjack(boolean printFlapjack) {
		this.printFlapjack = printFlapjack;
	}
	public void setPrintFlapjack(Boolean printFlapjack) {
		this.setPrintFlapjack(printFlapjack.booleanValue());
	}

	public boolean isPrintDarwin() {
		return printDarwin;
	}
	public void setPrintDarwin(boolean printDarwin) {
		this.printDarwin = printDarwin;
	}
	public void setPrintDarwin(Boolean printDarwin) {
		this.setPrintDarwin(printDarwin.booleanValue());
	}
	
	public boolean isPrintTreeMix() {
		return printTreeMix;
	}
	public void setPrintTreeMix(boolean printTreeMix) {
		this.printTreeMix = printTreeMix;
	}
	public void setPrintTreeMix(Boolean printTreeMix) {
		this.setPrintTreeMix(printTreeMix.booleanValue());
	}

	public boolean isPrintJoinMap() {
		return printJoinMap;
	}
	public void setPrintJoinMap(boolean printJoinMap) {
		this.printJoinMap = printJoinMap;
	}
	public void setPrintJoinMap(Boolean printJoinMap) {
		this.setPrintJoinMap(printJoinMap.booleanValue());
	}
	
	public boolean isPrintPhase() {
		return printPhase;
	}
	public void setPrintPhase(boolean printPhase) {
		this.printPhase = printPhase;
	}
	public void setPrintPhase(Boolean printPhase) {
		this.setPrintPhase(printPhase.booleanValue());
	}

	public boolean isPrintFineStructure() {
		return printFineStructure;
	}
	public void setPrintFineStructure(boolean printFineStructure) {
		this.printFineStructure = printFineStructure;
	}
	public void setPrintFineStructure(Boolean printFineStructure) {
		this.setPrintFineStructure(printFineStructure.booleanValue());
	}
	
	public boolean isPrintGenePop() {
		return printGenePop;
	}
	public void setPrintGenePop(boolean printGenePop) {
		this.printGenePop = printGenePop;
	}
	public void setPrintGenePop(Boolean printGenePop) {
		this.setPrintGenePop(printGenePop.booleanValue());
	}
	
	public String getSequenceName() {
		return sequenceName;
	}
	public void setSequenceName(String sequenceName) {
		this.sequenceName = sequenceName;
	}
	
	public String getIdParent1() {
		return idParent1;
	}
	public void setIdParent1(String idParent1) {
		this.idParent1 = idParent1;
	}

	public String getIdParent2() {
		return idParent2;
	}
	public void setIdParent2(String idParent2) {
		this.idParent2 = idParent2;
	}
	
	public String getPopulationFile() {
		return populationFile;
	}
	public void setPopulationFile(String populationFile) {
		this.populationFile = populationFile;
	}
	
	
	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		VCFConverter instance = new VCFConverter();
		CommandsDescriptor.getInstance().loadOptions(instance, args);
		instance.run();
	}
	
	public void run() throws IOException {
		logParameters();
		if(outputPrefix == null) throw new IOException("The prefix of the output files is a required parameter");
		if(printTreeMix && populationFile==null) throw new IOException("The file with the description of the populations is required for conversion to TreeMix");
		if(printPhase && sequenceName==null) throw new IOException("The sequence name is required for conversion to phase");
		if(printFineStructure && sequenceName==null) throw new IOException("The sequence name is required for conversion to fineStructure");
		if(printJoinMap && (idParent1==null || idParent2==null)) throw new IOException("Parent ids are required for conversion to JoinMap");
		if(inputFile==null) {
			process(System.in,outputPrefix);
		} else {
			process(inputFile,outputPrefix);
		}
		log.info("Process finished");
	}
	
	private void logParameters() {
		ByteArrayOutputStream os = new ByteArrayOutputStream();
		PrintStream out = new PrintStream(os);
		if(inputFile != null) out.println("Input file: "+inputFile);
		else out.println("System standard input");
		if(outputPrefix != null) out.println("Prefix for output files: "+outputPrefix);
		out.println("Formats:");
		if (printDarwin) out.print(" darwin");
		if (printEigensoft) out.print(" eigensoft");
		if (printEmma) out.print(" emma");
		if (printFasta) out.print(" fasta");
		if (printFineStructure) out.print(" fineStructure");
		if (printFlapjack) out.print(" flapjack");
		if (printGWASPoly) out.print(" GWASPoly");
		if (printHaploview) out.print(" haploview");
		if (printHapmap) out.print(" hapmap");
		if (printJoinMap) out.print(" joinMap");
		if (printMatrix) out.print(" matrix");
		if (printPhase) out.print(" phase");
		
		if (printPlink) out.print(" plink");
		if (printPowerMarker) out.print(" powerMarker");
		if (printrrBLUP) out.print(" rrBLUP");
		if (printSpagedi) out.print(" spagedi");
		if (printStructure) out.print(" structure");
		if (printTreeMix) out.print(" treeMix");
		if (printGenePop) out.print(" GenePop");
		out.println();
		if (sequenceName!=null) out.println("Sequence name for phase or fineStructure: "+sequenceName);
		if (populationFile!=null) out.println("File with population assignments: "+populationFile);
		if (idParent1!=null) out.println("First parent: "+idParent1);
		if (idParent2!=null) out.println("Second parent: "+idParent2);
	   
	    log.info(""+os.toString());
	}
 
	public void process(String vcfFile, String prefix) throws IOException {
		try (VCFFileReader reader = new VCFFileReader(vcfFile)) {
			process(reader,prefix);
		}
	}
	public void process(InputStream is, String prefix) throws IOException {
		try (VCFFileReader reader = new VCFFileReader(is)) {
			process(reader,prefix);
		}
	}
	public void process(VCFFileReader reader, String prefix) throws IOException {
		
		PrintStream outMatrix = null;
		PrintStream outHapmap = null;
		PrintStream outGWASPoly = null;
		PrintStream outJoinMap = null;
		PrintStream outTreemix = null;
		//Load the matrix if at least one format need the matrix to be transposed
		boolean loadMatrix = printFasta || printStructure || printrrBLUP || printSpagedi || printEmma || printPlink || printHaploview || printPowerMarker || printFlapjack || printEigensoft || printDarwin || printFineStructure || printGenePop;
		boolean loadMatrixSeqName = printPhase;
		List<List<CalledGenomicVariant>> callsPerVariant = new ArrayList<List<CalledGenomicVariant>>();
		List<String> sampleIds = null;
		try {
			if(log!=null)reader.setLog(log);
			if(!printGWASPoly) reader.setLoadMode(VCFFileReader.LOAD_MODE_MINIMAL);
			VCFFileHeader header = reader.getHeader();
			sampleIds = header.getSampleIds();
			Iterator<VCFRecord> it = reader.iterator();
			//Print header for formats that do not need transposition
			if(printMatrix) {
				outMatrix = new PrintStream(prefix+"_genotypes.txt");
				printMatrixHeader (sampleIds,outMatrix);
			}
			if(printHapmap) {
				outHapmap = new PrintStream(prefix+"_hmp.txt");
				printHapmapHeader (sampleIds,outHapmap);
			}
			if(printGWASPoly) {
				outGWASPoly = new PrintStream(prefix+"_GWASPoly.txt");
				printGWASPolyHeader (sampleIds,outGWASPoly);
			}
			Map<String, List<Integer>> groupsWithSampleIdxs = null;
			if(printTreeMix) {
				outTreemix = new PrintStream(prefix+"_treemix.txt");
				SimpleSamplesFileHandler samplesFH = new SimpleSamplesFileHandler();
				Map<String,Sample> samplesMap = samplesFH.loadSamplesAsMap(populationFile);
				if(samplesMap.size()==0)throw new IOException("The populations file does not contain any sample");
				groupsWithSampleIdxs = Sample.getGroupsWithSampleIdxs(samplesMap, sampleIds);
				
				if(groupsWithSampleIdxs.size()==0)throw new IOException("The populations file does not contain a group id for any sample");
				printTreeMixHeader (groupsWithSampleIdxs.keySet(),outTreemix);
			}

			int ip1=-1;
			int ip2=-1;
			if(printJoinMap){
				if(idParent1==null) throw new IOException("Id parent 1 is mandatory for conversion to JoinMap");
				if(idParent2==null) throw new IOException("Id parent 2 is mandatory for conversion to JoinMap");
				ip1=header.getIndexSampleId(idParent1);
				if(ip1==-1) throw new IOException("Id parent "+idParent1+" not found in vcf file");
				ip2=header.getIndexSampleId(idParent2);
				if(ip2==-1) throw new IOException("Id parent "+idParent2+" not found in vcf file");
				if(ip1 == ip2) throw new IOException("The two parents must be different");
				outJoinMap=new PrintStream(prefix+"_joinmap.txt");
				printJoinMapHeader(sampleIds, outJoinMap, ip1, ip2);		
			}
					
			int n=0;
			while(it.hasNext()) {
				VCFRecord record = it.next();
				List<CalledGenomicVariant> calls = record.getCalls();
				//Adds the calls for the new variant if needed
				if(loadMatrix) callsPerVariant.add(calls);
				//PHASE converter. Loads the variants only for the sequence name requested by the user
				else if (loadMatrixSeqName && record.getVariant().getSequenceName().equals(sequenceName)) callsPerVariant.add(calls); 
				//Print genotypes for the current variant for formats that do not need transposition
				if(outMatrix!=null)printGenotypesMatrix(record,outMatrix);
				if(outHapmap!=null)printHapmap(record,outHapmap);
				if(outGWASPoly!=null)printGWASPoly(record,outGWASPoly,(n+1));
				if(outTreemix!=null && record.getVariant().isBiallelic()) printTreeMix(groupsWithSampleIdxs,record.getCalls(),outTreemix);
				if(outJoinMap!=null)printJoinMap(record, outJoinMap, ip1, ip2);
				n++;
				if (progressNotifier!=null && n%1000==0) {
					int progress = n/1000;
					if (!progressNotifier.keepRunning(progress)) {
						return;
					}
				}
			}
		} finally {
			if(outMatrix!=null) {
				outMatrix.flush();
				outMatrix.close();
			}
			if(outHapmap!=null) {
				outHapmap.flush();
				outHapmap.close();
			}
			if(outGWASPoly!=null) {
				outGWASPoly.flush();
				outGWASPoly.close();
			}
			if (outTreemix!=null) {
				outTreemix.flush();
				outTreemix.close();
			}
			if (outJoinMap!=null) {
				outJoinMap.flush();
				outJoinMap.close();
			}
		}
		
		
		
		if(!loadMatrix && !loadMatrixSeqName) return;
		System.out.println("Loaded "+callsPerVariant.size()+" calls");
		//Print formats that require transposition
		if(printFasta) printFasta(sampleIds,callsPerVariant,prefix+"_aln.fa");
		if(printDarwin)printDarwin(sampleIds, callsPerVariant, prefix);
		if(printStructure) printStructure(sampleIds,callsPerVariant,prefix+"_structure.in");
		if(printrrBLUP) printrrBLUP(sampleIds,callsPerVariant,prefix);
		if(printSpagedi) printSpagedi(sampleIds, callsPerVariant, prefix+"_spagedi.in");
		if(printEmma) printEmma(callsPerVariant, prefix+"_emma.in");
		if(printPlink) printPlink(sampleIds,callsPerVariant,prefix+"_plink",true);
		if(printHaploview) printPlink(sampleIds,callsPerVariant,prefix+"_haploview",false);
		if(printPowerMarker) printPowerMarker(sampleIds,callsPerVariant,prefix);
		if(printFlapjack) printFlapjack(sampleIds,callsPerVariant,prefix);
		if(printEigensoft) printEigensoft(sampleIds,callsPerVariant,prefix);
		if(printGenePop) printGenePop(sampleIds,callsPerVariant,prefix+"_genePop.in");
		if(printPhase) printPhase(sampleIds,callsPerVariant,prefix+"_"+sequenceName+"_phase.inp");
		if(printFineStructure) printFineStructure(sampleIds.size(),callsPerVariant,prefix+"_fs_"+sequenceName+".phase");
		
	}
	
	private void printFlapjack(List<String> sampleIds,List<List<CalledGenomicVariant>> calls, String outPrefix) throws IOException {
		StringBuilder [] sequences = new StringBuilder[sampleIds.size()];
		for(int i=0;i<sequences.length;i++) sequences[i] = new StringBuilder(sampleIds.get(i));
		
		PrintStream outMap = new PrintStream(outPrefix+"_flapjack.map");
		PrintStream outGenotypes = new PrintStream(outPrefix+"_flapjack.gen");
		outMap.println("# fjFile = MAP");
		outGenotypes.println("# fjFile = GENOTYPE");
		int snpId=1;
		for(List<CalledGenomicVariant> callsVariant:calls) {
			if(callsVariant.size()>0 && !callsVariant.get(0).isBiallelic()) continue;
			for(int i=0;i<callsVariant.size();i++) {
				CalledGenomicVariant calledVar = callsVariant.get(i);
				
				String [] alleles = calledVar.getCalledAlleles(); 
				String allele = "-";
				if(!calledVar.isUndecided()) {
					if(alleles.length>=2) {
						allele = ""+alleles[0]+"/"+alleles[1];
					} else {
						allele = ""+alleles[0];
					}
				}
				sequences[i].append("\t"+allele);
				if(i==0) {
					String markerId = "Marker"+snpId;
					outMap.println(markerId+"\t"+calledVar.getSequenceName()+"\t"+calledVar.getFirst());
					outGenotypes.print("\t"+markerId);
					snpId++;
				}
			}
		}
		outGenotypes.println();
		for(int i=0;i<sequences.length;i++) {
			outGenotypes.println(sequences[i].toString());
		}
		outMap.flush();
		outMap.close();
		outGenotypes.flush();
		outGenotypes.close();	
	}
	
	private void printPowerMarker(List<String> sampleIds,List<List<CalledGenomicVariant>> calls, String outPrefix) throws IOException {
		StringBuilder [] sequences = new StringBuilder[sampleIds.size()];
		for(int i=0;i<sequences.length;i++) sequences[i] = new StringBuilder(sampleIds.get(i));
		
		PrintStream out = new PrintStream(outPrefix+"_powermarker.in");
		PrintStream outSNPs = new PrintStream(outPrefix+"_powermarker.snp");
		out.print("Sample");
		int snpId=1;
		for(List<CalledGenomicVariant> callsVariant:calls) {
			if(callsVariant.size()>0 && !callsVariant.get(0).isBiallelic()) continue;
			for(int i=0;i<callsVariant.size();i++) {
				CalledGenomicVariant calledVar = callsVariant.get(i);
				String allele = "-9\t-9";
				if(!calledVar.isUndecided()) {
					if(calledVar.isHeterozygous()) {
						allele = "0\t1";
					} else if (calledVar.isHomozygousReference()) {
						allele = "0\t0";
					} else {
						allele = "1\t1";
					}
				}
				sequences[i].append("\t"+allele);
				if(i==0) {
					out.print("\t"+snpId+"\t"+snpId);
					outSNPs.println(snpId+"\t"+calledVar.getSequenceName()+"\t"+calledVar.getFirst());
					snpId++;
				}
			}
		}
		out.println();
		for(int i=0;i<sequences.length;i++) {
			out.println(sequences[i].toString());
		}
		out.flush();
		out.close();
		outSNPs.flush();
		outSNPs.close();
		
	}
	private void printPlink(List<String> sampleIds,List<List<CalledGenomicVariant>> calls, String prefix, boolean includeChrInMap) throws IOException {
		StringBuilder [] sequences = new StringBuilder[sampleIds.size()];
		for(int i=0;i<sequences.length;i++) sequences[i] = new StringBuilder(sampleIds.get(i)+"\t"+sampleIds.get(i)+"\t0\t0\t0\t0");
		PrintStream outMap;
		if(includeChrInMap) outMap = new PrintStream(prefix+".map");
		else outMap = new PrintStream(prefix+".info");
		String currentChr = null;
		int chrN = 0;
		for(List<CalledGenomicVariant> callsVariant:calls) {
			if(callsVariant.size()==0) continue;
			CalledGenomicVariant firstCall = callsVariant.get(0);
			if(!firstCall.isBiallelic() || !firstCall.isSNV() ) continue;
			for(int i=0;i<callsVariant.size();i++) {
				CalledGenomicVariant calledVar = callsVariant.get(i);
				if(i==0) {
					String seqName = calledVar.getSequenceName();
					if(!seqName.equals(currentChr)) {
						chrN++;
						currentChr = seqName;
					}
					if(includeChrInMap) outMap.print(""+chrN+"\t");
					outMap.println(calledVar.getSequenceName()+"_"+calledVar.getFirst()+"\t"+calledVar.getFirst());
				}
				String [] calledAlleles = calledVar.getCalledAlleles(); 
				String allele = "0 0";
				if(!calledVar.isUndecided()) {
					if(calledAlleles.length>=2) {
						allele = ""+calledAlleles[0]+" "+calledAlleles[1];
					} else {
						allele = ""+calledAlleles[0]+" "+calledAlleles[0];
					}
				}
				sequences[i].append("\t"+allele);
			}
		}
		outMap.flush();
		outMap.close();
		PrintStream out = new PrintStream(prefix+".ped");
		for(int i=0;i<sequences.length;i++) {
			out.println(sequences[i].toString());
		}
		out.flush();
		out.close();
		
	}
	private void printEigensoft(List<String> sampleIds, List<List<CalledGenomicVariant>> selectedCalls, String outPrefix) throws IOException {
		PrintStream outGenotype = new PrintStream(outPrefix+"_eigensoft.geno");
		PrintStream outSNP = new PrintStream(outPrefix+"_eigensoft.snp");
		int chrNum = 0;
		String lastChr = "";
		for(List<CalledGenomicVariant> varCalls:selectedCalls) {
			if(varCalls.size()>0 && !varCalls.get(0).isBiallelic()) continue;
			for(int i=0;i<varCalls.size();i++) {
				CalledGenomicVariant calledVar =varCalls.get(i);
				if(calledVar.isUndecided()) outGenotype.print("9");
				else if (calledVar.isHomozygousReference()) outGenotype.print("2");
				else if (calledVar.isHeterozygous()) outGenotype.print("1");
				else outGenotype.print("0");
				if(i==0) {
					if(!lastChr.equals(calledVar.getSequenceName())) {
						chrNum++;
						lastChr = calledVar.getSequenceName();
					}
					outSNP.printf("%30s %20d      0.0%20d %1s %1s\n",calledVar.getSequenceName()+"_"+calledVar.getFirst(),chrNum,calledVar.getFirst(),calledVar.getReference(),calledVar.getAlleles()[1]);
				}
			}
			outGenotype.println();
		}
		outGenotype.flush();
		outGenotype.close();
		outSNP.flush();
		outSNP.close();
		PrintStream outInd = new PrintStream(outPrefix+"_eigensoft.ind");
		for(String sampleId:sampleIds) {
			outInd.printf("%39s U Control\n", sampleId);
		}
		outInd.flush();
		outInd.close();
	}
	private void printEmma(List<List<CalledGenomicVariant>> selectedCalls, String outFile) throws IOException {
		PrintStream out = new PrintStream(outFile);
		for(List<CalledGenomicVariant> varCalls:selectedCalls) {
			if(varCalls.size()>0 && !varCalls.get(0).isBiallelic()) continue;
			for(int i=0;i<varCalls.size();i++) {
				if(i>0) out.print("\t");
				CalledGenomicVariant calledVar =varCalls.get(i);
				if(calledVar.isUndecided()) out.print("NA");
				else if (calledVar.isHomozygousReference()) out.print("0");
				else if (calledVar.isHeterozygous()) out.print("0.5");
				else out.print("1");
			}
			out.println();
		}
		out.flush();
		out.close();
		
	}
	private void printSpagedi(List<String> sampleIds,List<List<CalledGenomicVariant>> calls, String outFile) throws IOException {
		StringBuilder [] sequences = new StringBuilder[sampleIds.size()];
		for(int i=0;i<sequences.length;i++) sequences[i] = new StringBuilder(sampleIds.get(i));
		
		for(List<CalledGenomicVariant> callsVariant:calls) {
			if(callsVariant.size()>0 && !callsVariant.get(0).isBiallelic()) continue;
			for(int i=0;i<callsVariant.size();i++) {
				CalledGenomicVariant calledVar = callsVariant.get(i);
				String allele = "0,0";
				if(!calledVar.isUndecided()) {
					if(calledVar.isHeterozygous()) {
						allele = "1,2";
					} else if (calledVar.isHomozygousReference()) {
						allele = "1,1";
					} else {
						allele = "2,2";
					}
				}
				sequences[i].append("\t"+allele);
			}
		}
		PrintStream out = new PrintStream(outFile);
		out.println(""+sampleIds.size()+"\t0\t0\t"+calls.size()+"\t1\t2");
		out.println("0");
		out.print("Ind");
		for(int i=0;i<calls.size();i++) {
			out.print("\tSNP_"+(i+1));
		}
		out.println();
		for(int i=0;i<sequences.length;i++) {
			out.println(sequences[i].toString());
		}
		out.println("END");
		out.flush();
		out.close();
	}
	
	private void printMatrixHeader(List<String> sampleIds, PrintStream out) {
		out.print("Chr\tPosition\tId\tReference\tAlternative");
		for(String sampleId:sampleIds) {
			out.print("\t"+sampleId);
		}
		out.println();
		
	}
	
	private void printGenotypesMatrix(VCFRecord record, PrintStream out) throws IOException {
		GenomicVariant var = record.getVariant();
		out.print(var.getSequenceName());
		out.print("\t"+var.getFirst());
		out.print("\t");
		if(var.getId()!=null) out.print(var.getId());
		else out.print("-");
		String [] alleles = var.getAlleles();
		out.print("\t"+alleles[0]);
		out.print("\t");
		for(int j=1;j<alleles.length;j++) {
			if(j>1) out.print(",");
			out.print(alleles[j]);
		}
		for(int i=0;i<record.getCalls().size();i++) {
			CalledGenomicVariant calledVar =record.getCalls().get(i);
			String [] calledAlleles = calledVar.getCalledAlleles();
			Set<String> printedAlleles = new TreeSet<String>();
			out.print("\t");
			if(calledAlleles.length==0) out.print("N");
			for(int j=0;j<calledAlleles.length;j++) {
				if(printedAlleles.contains(calledAlleles[j])) continue;
				if(printedAlleles.size()>0) out.print("/");
				out.print(calledAlleles[j]);
				printedAlleles.add(calledAlleles[j]);
			}
		}
		out.println();
	}
	private void printHapmapHeader(List<String> sampleIds, PrintStream out) {
		out.print("rs#\talleles\tchrom\tpos\tstrand\tassembly#\tcenter\tprotLSID\tassayLSID\tpanel\tQCcode");
		for(String sampleId:sampleIds) {
			out.print("\t"+sampleId);
		}
		out.println();
	}
	private void printHapmap(VCFRecord record, PrintStream out) throws IOException{
		GenomicVariant var = record.getVariant();
		if(!var.isBiallelic()) return;
		String [] alleles = var.getAlleles();
		if(var.getId()!=null) out.print(var.getId());
		else out.print(var.getSequenceName()+"pos"+var.getFirst()+".1");
		out.print("\t"+alleles[0]+"/"+alleles[1]);
		out.print("\t"+var.getSequenceName());
		out.print("\t"+var.getFirst());
		
		out.print("\t+\tNA\tNA\tNA\tNA\tNA\tNA");
		for(int i=0;i<record.getCalls().size();i++) {
			CalledGenomicVariant calledVar =record.getCalls().get(i);
			if(calledVar.isUndecided()) out.print("\tNN");
			else if (calledVar.isHomozygousReference()) out.print("\t"+alleles[0]+""+alleles[0]);
			else if (calledVar.isHeterozygous()) out.print("\t"+alleles[0]+""+alleles[1]);
			else out.print("\t"+alleles[1]+""+alleles[1]);
		}
		out.println();	
	}
	
	private void printGWASPolyHeader(List<String> sampleIds, PrintStream out) {
		out.print("Marker,Chrom,Position");
		for(String sampleId:sampleIds) {
			out.print(","+sampleId);
		}
		out.println();
	}

	private void printGWASPoly(VCFRecord record, PrintStream out, int index) {
		GenomicVariant var = record.getVariant();
		if(!(var instanceof SNV)) return;
		
		String [] alleles = record.getVariant().getAlleles();
		out.print(""+index+","+record.getSequenceName()+","+record.getFirst());
		List<CalledGenomicVariant> calls = record.getCalls();
		for(int i=0;i<calls.size();i++) {
			CalledGenomicVariant calledVar = calls.get(i);
			short [] allelesCN = calledVar.getAllelesCopyNumber();
			if(calledVar.isUndecided()) out.print(",NA");
			else {
				out.print(",");
				for(int j=0;j<allelesCN.length;j++) {
					for(int k=0;k<allelesCN[j];k++) {
						out.print(alleles[j]);
					}
				}
			}
		}
		out.println();
	}
		
	private void printStructure(List<String> sampleIds,List<List<CalledGenomicVariant>> calls, String outFile) throws IOException {
		StringBuilder [] sequences = new StringBuilder[sampleIds.size()];
		for(int i=0;i<sequences.length;i++) sequences[i] = new StringBuilder(sampleIds.get(i));
		
		for(List<CalledGenomicVariant> callsVariant:calls) {
			if(callsVariant.size()>0 && !callsVariant.get(0).isBiallelic()) continue;
			for(int i=0;i<callsVariant.size();i++) {
				CalledGenomicVariant calledVar = callsVariant.get(i);
				String allele = "-9 -9";
				if(!calledVar.isUndecided()) {
					if(calledVar.isHeterozygous()) {
						allele = "0 1";
					} else if (calledVar.isHomozygousReference()) {
						allele = "0 0";
					} else {
						allele = "1 1";
					}
				}
				sequences[i].append(" "+allele);
			}
		}
		PrintStream out = new PrintStream(outFile);
		for(int i=0;i<sequences.length;i++) {
			out.println(sequences[i].toString());
		}
		out.flush();
		out.close();
	}
	
	private void printGenePop(List<String> sampleIds,List<List<CalledGenomicVariant>> calls, String outFile) throws IOException {
		StringBuilder [] sequences = new StringBuilder[sampleIds.size()];
		for(int i=0;i<sequences.length;i++) sequences[i] = new StringBuilder(sampleIds.get(i)+",");
		String [] formattedAlleles = new String[100];
		for(int a=0;a<10;a++) formattedAlleles[a] = "0"+a;
		for(int a=10;a<100;a++) formattedAlleles[a] = ""+a;
		Map<String, List<Integer>> groupsWithSampleIdxs;
		if(populationFile!=null) {
			SimpleSamplesFileHandler samplesFH = new SimpleSamplesFileHandler();
			Map<String,Sample> samplesMap = samplesFH.loadSamplesAsMap(populationFile);
			if(samplesMap.size()==0)throw new IOException("The populations file does not contain any sample");
			groupsWithSampleIdxs = Sample.getGroupsWithSampleIdxs(samplesMap, sampleIds);
			if(groupsWithSampleIdxs.size()==0)throw new IOException("The populations file does not contain any group");
		} else {
			groupsWithSampleIdxs = new HashMap<>();
			List<Integer> allIds = new ArrayList<>();
			for(int i=0;i<sampleIds.size();i++) allIds.add(i);
			groupsWithSampleIdxs.put("All", allIds);
		}
		
		try (PrintStream out = new PrintStream(outFile)) {
			out.println("GenePopFile");
			for(List<CalledGenomicVariant> callsVariant:calls) {
				if(callsVariant.size()>0 && callsVariant.get(0).getAlleles().length>99) continue;
				
				for(int i=0;i<callsVariant.size();i++) {
					CalledGenomicVariant calledVar = callsVariant.get(i);
					if(i==0) out.println(calledVar.getSequenceName()+":"+calledVar.getFirst());
					String allele = "0000";
					if(!calledVar.isUndecided()) {
						byte [] idxsAlleles = calledVar.getIndexesCalledAlleles();
						if(idxsAlleles.length==1) {
							String a0 = formattedAlleles[idxsAlleles[0]+1];
							allele = a0 + a0;
						} else if(idxsAlleles.length>1) {
							String a0 = formattedAlleles[idxsAlleles[0]+1];
							String a1 = formattedAlleles[idxsAlleles[1]+1];
							allele = a0 + a1;
						}
					}
					sequences[i].append(" "+allele);
				}
			}
			for(List<Integer> ids:groupsWithSampleIdxs.values()) {
				out.println("POP");
				for(int i:ids) {
					out.println(sequences[i].toString());
				}
			}
			
		}
	}
	
	private void printrrBLUP(List<String> sampleIds,List<List<CalledGenomicVariant>> calls, String prefix) throws IOException {
		String outFile = prefix+"_rrBLUP.in";
		String outSamples = prefix+"_rrBLUP_samples.txt";
		PrintStream outS = new PrintStream(outSamples);
		for(String sampleId:sampleIds) {
			outS.println(sampleId);
		}
		outS.flush();
		outS.close();
		StringBuilder [] sequences = new StringBuilder[sampleIds.size()];
		for(int i=0;i<sequences.length;i++) sequences[i] = new StringBuilder();
		
		for(List<CalledGenomicVariant> callsVariant:calls) {
			if(callsVariant.size()>0 && !callsVariant.get(0).isBiallelic()) continue;
			for(int i=0;i<callsVariant.size();i++) {
				CalledGenomicVariant calledVar = callsVariant.get(i);
				String allele = "NA";
				if(!calledVar.isUndecided()) {
					if(calledVar.isHeterozygous()) {
						allele = "0";
					} else if (calledVar.isHomozygousReference()) {
						allele = "1";
					} else {
						allele = "-1";
					}
				}
				if(i>=0) sequences[i].append(" ");
				sequences[i].append(allele);
			}
		}
		PrintStream out = new PrintStream(outFile);
		for(int i=0;i<sequences.length;i++) {
			out.println(sequences[i].toString());
		}
		out.flush();
		out.close();
	}
	
	private void printFasta(List<String> sampleIds, List<List<CalledGenomicVariant>> calls, String outFile) throws IOException {
		StringBuilder [] sequences = new StringBuilder[sampleIds.size()];
		for(int i=0;i<sequences.length;i++) sequences[i] = new StringBuilder();
		for(List<CalledGenomicVariant> callsVariant:calls) {
			if(callsVariant.size()>0 && !callsVariant.get(0).isBiallelic()) continue;
			for(int i=0;i<callsVariant.size();i++) {
				CalledGenomicVariant calledVar = callsVariant.get(i);
				String [] alleles = calledVar.getAlleles();
				char call = 'N';
				if(calledVar.isHomozygousReference()) {
					call = alleles[0].charAt(0);
				} else if (calledVar.isHomozygous()) {
					call = alleles[1].charAt(0);
				} 
				sequences[i].append(call);
			}
		}
		FastaSequencesHandler handler = new FastaSequencesHandler();
		List<QualifiedSequence> seqsList = new ArrayList<QualifiedSequence>();
		for(int i=0;i<sequences.length;i++) {
			String seqId = sampleIds.get(i);
			DNAMaskedSequence seq = new DNAMaskedSequence(sequences[i].toString());
			seqsList.add(new QualifiedSequence(seqId, seq));
		}
		PrintStream out = new PrintStream(outFile);
		handler.saveSequences(seqsList, out, 100);
		out.flush();
		out.close();
	}
	
	private void printPhase(List<String> sampleIds,List<List<CalledGenomicVariant>> calls, String outFile) throws IOException {
		StringBuilder [] sequences = new StringBuilder[2*sampleIds.size()];
		for(int i=0;i<sequences.length;i++) sequences[i] = new StringBuilder();
		StringBuilder linePositions = new StringBuilder("P");
		int nSites =0;
		for(List<CalledGenomicVariant> callsVariant:calls) {
			if(callsVariant.size()>0 && !callsVariant.get(0).isBiallelic()) continue;
			GenomicVariant variant = callsVariant.get(0);
			if(!variant.getSequenceName().equals(sequenceName)) continue;
			nSites++;
			linePositions.append(" "+variant.getFirst());
			for(int i=0;i<callsVariant.size();i++) {
				CalledGenomicVariant calledVar = callsVariant.get(i);
				String allele1 = "?";
				String allele2 = "?";
				if(!calledVar.isUndecided()) {
					if(calledVar.isHeterozygous()) {
						allele1 = "0";
						allele2 = "1";
					} else if (calledVar.isHomozygousReference()) {
						allele1 = "0";
						allele2 = "0";
					} else {
						allele1 = "1";
						allele2 = "1";
					}
				}
				sequences[2*i].append(allele1);
				sequences[2*i+1].append(allele2);
			}
		}
		if(nSites==0) throw new IOException("No biallelic variants found for the given sequence name");
		try (PrintStream out = new PrintStream(outFile)) {
			out.println(sampleIds.size());
			out.println(nSites);
			out.println(linePositions);
			for(int i=0;i<nSites;i++) out.print("S");
			out.println();
			for(int i=0;i<sampleIds.size();i++) {
				out.println("#"+sampleIds.get(i));
				out.println(sequences[2*i].toString());
				out.println(sequences[2*i+1].toString());
			}
		}
	}
	
	private void printFineStructure(int numSamples, List<List<CalledGenomicVariant>> calls, String outFile) throws IOException {
		//log.info("Converting to fineStructure. Variants: "+calls.size()+" samples: "+numSamples);
		StringBuilder [] sequences = new StringBuilder[2*numSamples];
		for(int i=0;i<sequences.length;i++) sequences[i] = new StringBuilder();
		StringBuilder linePositions = new StringBuilder("P");
		int nSites =0;
		for(List<CalledGenomicVariant> callsVariant:calls) {
			if(callsVariant.size()==0 || !callsVariant.get(0).isBiallelic()) continue;
			GenomicVariant variant = callsVariant.get(0);
			if(!variant.getSequenceName().equals(sequenceName)) continue;
			if(callsVariant.size()!=numSamples) {
				log.warning("Number of samples for variant "+variant.getSequenceName()+":"+variant.getFirst()+" does not match expected number: "+numSamples);
				continue;
			}
			linePositions.append(" "+variant.getFirst());
			for(int i=0;i<callsVariant.size();i++) {
				CalledGenomicVariant calledVar = callsVariant.get(i);
				byte [] idxsPhasedAlleles = calledVar.getIndexesPhasedAlleles();
				String allele1 = "0";
				String allele2 = "0";
				if(!calledVar.isUndecided()) {
					if(calledVar.isHeterozygous()) {
						if(idxsPhasedAlleles!=null && idxsPhasedAlleles.length==2) {
							allele1 = ""+idxsPhasedAlleles[0];
							allele2 = ""+idxsPhasedAlleles[1];
						} else {
							log.warning("Unphased heterozygous call at variant at "+variant.getSequenceName()+":"+variant.getFirst()+" indivudual: "+i);
							allele1 = "0";
							allele2 = "1";
						}
		
					} else if (!calledVar.isHomozygousReference()) {
						allele1 = "1";
						allele2 = "1";
					}
				} else {
					log.warning("Imputing reference genotype for undecided call at variant at "+variant.getSequenceName()+":"+variant.getFirst()+" indivudual: "+i);
				}
				sequences[2*i].append(allele1);
				sequences[2*i+1].append(allele2);
			}
			
			nSites++;
		}
		if(nSites==0) throw new IOException("No biallelic variants found for the given sequence name");
		try (PrintStream out = new PrintStream(outFile)) {
			out.println(2*numSamples);
			out.println(nSites);
			out.println(linePositions);
			for(int i=0;i<numSamples;i++) {
				out.println(sequences[2*i].toString());
				out.println(sequences[2*i+1].toString());
			}
		}
	}
	
	private void printTreeMixHeader(Set<String> groups, PrintStream out) throws IOException {
		int i=0;
		for (String group:groups) {
			if(i>0) out.print("\t");
			out.print(group);
			i++;
		}
		out.println();
	}
	
	private void printTreeMix(Map<String, List<Integer>> groups, List<CalledGenomicVariant> calls, PrintStream out) throws IOException {
		int j=0;
		for(String group:groups.keySet()) {
			List<Integer> groupIdxs = groups.get(group);
			int count0 = 0;
			int count1 = 0;
			for(int i:groupIdxs) {
				CalledGenomicVariant calledVar = calls.get(i);
				if(calledVar.isUndecided()) continue;
				if(calledVar.isHeterozygous()) {
					count0 ++;
					count1 ++;
				} else if (calledVar.isHomozygousReference()) {
					count0 += 2 ;
				} else {
					count1 += 2;
				}
			}
			if(j>0) out.print("\t");
			out.print(count0+","+count1);
			j++;
		}
		out.println();
	}

	
	
	private void printDarwin(List<String> sampleIds,List<List<CalledGenomicVariant>> calls, String outPrefix) throws IOException {
		StringBuilder [] sequences = new StringBuilder[sampleIds.size()];
		for(int i=0;i<sequences.length;i++) sequences[i] = new StringBuilder(sampleIds.get(i));
		PrintStream out = new PrintStream(outPrefix+"_Darwin.in");
		PrintStream outSNPs = new PrintStream(outPrefix+"_Darwin.snp");
		out.println(sampleIds.size()+"\t"+(2*calls.size()));
		out.print("N");
		int snpId=1;
		for(List<CalledGenomicVariant> callsVariant:calls) {
			for(int i=0;i<callsVariant.size();i++) {
				CalledGenomicVariant calledVar = callsVariant.get(i);
				String allele = "-9\t-9";
				if(!calledVar.isUndecided()) {
					if(calledVar.isHeterozygous()) {
						allele = "0\t1";
					} else if (calledVar.isHomozygousReference()) {
						allele = "0\t0";
					} else {
						allele = "1\t1";
					}
				}
				sequences[i].append("\t"+allele);
				if(i==0) {
					out.print("\t"+snpId+"\t"+snpId);
					outSNPs.println(snpId+"\t"+calledVar.getSequenceName()+"\t"+calledVar.getFirst());
					snpId++;
				}
			}
		}
		out.println();
		for(int i=0;i<sequences.length;i++) {
			out.println(sequences[i].toString());
		}
		out.flush();
		out.close();
		outSNPs.flush();
		outSNPs.close();
	}
	
	
	private void printJoinMapHeader(List<String> sampleIds, PrintStream out, int parent1, int parent2){
		out.print("SNPID\tSegregation\tClasification\t"+sampleIds.get(parent1)+"\t"+sampleIds.get(parent2));
		for (int i=0;i<sampleIds.size();i++) {
			if(i!=parent1 && i != parent2) {
				out.print("\t"+sampleIds.get(i));
			}
		}
		out.println();
	}
	
	private void printJoinMap(VCFRecord record,PrintStream out, int iPJ1, int iPJ2) throws IOException{
		GenomicVariant var = record.getVariant();
		if(!(var instanceof SNV)) return;
		String unknownAllele = "--";
		String errorAllele = "E";
		List<CalledGenomicVariant> genotypeCalls = record.getCalls();
		CalledGenomicVariant genP1 = genotypeCalls.get(iPJ1);
		CalledGenomicVariant genP2 = genotypeCalls.get(iPJ2);
		if(genP1.isUndecided() || genP2.isUndecided()) return;
		if(genP1.isHomozygous() && genP2.isHomozygous()) return;
		out.print(""+var.getSequenceName()+"_"+var.getFirst());
		String hetAllele;
		String homoAllele;
		String homoAllele2;
		if(genP1.isHeterozygous() && genP2.isHomozygous()) {
			if(genP2.isHomozygousReference()){
				hetAllele = "lm";
				homoAllele = "ll";
				homoAllele2= errorAllele;
				out.print("\t<lmxll>\t(ll,lm)\t"+hetAllele+"\t"+homoAllele);
			}else{
				hetAllele = "lm";
				homoAllele2 = "ll";
				homoAllele = errorAllele;
				out.print("\t<lmxll>\t(ll,lm)\t"+hetAllele+"\t"+homoAllele2);
			}
		} else if (genP1.isHomozygous() && genP2.isHeterozygous()){
			if(genP1.isHomozygousReference()){
				hetAllele = "np";
				homoAllele = "nn";
				homoAllele2 = errorAllele;
				out.print("\t<nnxnp>\t(nn,np)\t"+homoAllele+"\t"+hetAllele);
			}else{
				hetAllele = "np";
				homoAllele2 = "nn";
				homoAllele = errorAllele;
				out.print("\t<nnxnp>\t(nn,np)\t"+homoAllele2+"\t"+hetAllele);
			}		
		} else {
			hetAllele = "hk";
			homoAllele = "hh";
			homoAllele2 = "kk";
			out.print("\t<hkxhk>\t(hh,hk,kk)\t"+hetAllele+"\t"+hetAllele);
		}
		for(int i=0;i<genotypeCalls.size();i++) {
			if(i==iPJ1 || i == iPJ2) continue;
			CalledGenomicVariant genotypeCall = genotypeCalls.get(i);
			out.print("\t");
			if(genotypeCall.isUndecided()) {
				out.print(unknownAllele);
			} else if (genotypeCall.isHeterozygous()) {
				out.print(hetAllele);
			} else if (genotypeCall.isHomozygous()){
				String allele = homoAllele;
				if (!genotypeCall.isHomozygousReference()) allele = homoAllele2;
				if(allele!=errorAllele) out.print(allele);
				else {
					log.warning("Inconsistent homozygous genotype call found at "+genotypeCall.getSequenceName()+":"+genotypeCall.getFirst()+"-"+genotypeCall.getLast()+" sample: "+genotypeCall.getSampleId());
					out.print(unknownAllele);
				}
			} 
		}
		out.println();
	}
}
