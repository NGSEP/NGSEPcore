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
import java.io.InputStream;
import java.io.PrintStream;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.logging.Logger;

import ngsep.main.CommandsDescriptor;
import ngsep.variants.CalledGenomicVariant;
import ngsep.variants.DiversityStatistics;
import ngsep.variants.GenomicVariant;
import ngsep.variants.Sample;
import ngsep.variants.io.SimpleSamplesFileHandler;


public class VCFDiversityCalculator {
	
	public static final String COMMAND_DIVERSITYSTATS = "DiversityStats";
	private Logger log = Logger.getLogger(VCFDiversityCalculator.class.getName());
	private Map<String,Sample> samplesMap=null;
	private boolean assumeAlwaysDiploid = false;
	private DecimalFormat fmt = new DecimalFormat("#0.0000",DecimalFormatSymbols.getInstance(Locale.ENGLISH));
	
	public static void main(String[] args) throws Exception {
		VCFDiversityCalculator instance = new VCFDiversityCalculator();
		int i = CommandsDescriptor.getInstance().loadOptions(instance, args);
		String vcfFile = args[i++];
		String samplesFile = null;
		if(i<args.length) {
			samplesFile = args[i++];
		}
		
		if(samplesFile!=null) instance.loadSamplesFile(samplesFile);
		if("-".equals(vcfFile)) instance.processFile(System.in,System.out);
		else instance.processFile(vcfFile,System.out);
	}
	
	public Logger getLog() {
		return log;
	}


	public void setLog(Logger log) {
		this.log = log;
	}


	public void loadSamplesFile(String samplesFile) throws IOException {
		SimpleSamplesFileHandler handler = new SimpleSamplesFileHandler();
		samplesMap = handler.loadSamplesAsMap(samplesFile);
	}
	public void processFile(String vcfFile, PrintStream out) throws IOException {
		VCFFileReader in = null;
		try {
			in = new VCFFileReader(vcfFile);
			processFile(in, out);
		} finally {
			if(in!=null) in.close(); 
		}
	}
	public void processFile(InputStream input, PrintStream out) throws IOException {
		VCFFileReader in = null;
		try {
			in = new VCFFileReader(input);
			processFile(in, out);
		} finally {
			if(in!=null) in.close(); 
		}
	}
	public void processFile(VCFFileReader in, PrintStream out) throws IOException {
		
		if(log!=null)in.setLog(log);
		in.setLoadMode(VCFFileReader.LOAD_MODE_COPY_NUMBER);
		List<String> sampleIdsVCF = in.getSampleIds();
		Map<String, List<Integer>> groupSampleIdxs = Sample.getGroupsWithSampleIdxs(samplesMap, sampleIdsVCF);
		printDefaultHeader(out);
		out.print("Chr\tfirst\tLast\tRef\tAll");
		for(String groupId:groupSampleIdxs.keySet()) {
			out.print("\t"+groupId);
		}
		out.println();
		Iterator<VCFRecord> it = in.iterator();
		while(it.hasNext()) {
			VCFRecord record = it.next();
			List<CalledGenomicVariant> allCalls = record.getCalls();
			List<DiversityStatistics> groupsStats = new ArrayList<DiversityStatistics>();
			groupsStats.add(DiversityStatistics.calculateDiversityStatistics(allCalls, assumeAlwaysDiploid));
			for(String groupId:groupSampleIdxs.keySet()) {
					List<Integer> indexes = groupSampleIdxs.get(groupId);
					List<CalledGenomicVariant> groupCalls = selectCalls(allCalls,indexes);
					groupsStats.add(DiversityStatistics.calculateDiversityStatistics(groupCalls, assumeAlwaysDiploid));
			}
			printStats (record.getVariant(),groupsStats,out);
		}
		
		out.flush();
	}

	private void printDefaultHeader(PrintStream out) {
		out.println("#Diversity statisitcs. Fields calculated per population and separated by semicolon:");
		out.println("#1. Number of samples genotyped");
		out.println("#2. Expected heterozygosity (under HWE)");
		out.println("#3. Observed heterozygosity");
		out.println("#4. F-statistic (1-OH/EH)");
		out.println("#5. Minor allele frequency (MAF)");
		out.println("#6. Reference allele frequency");
		out.println("#7. Chi-square value of departure from HWE");
		out.println("#8. Uncorrected p-value of the Chi-square test for departure from HWE");
	}

	private void printStats(GenomicVariant v,List<DiversityStatistics> groupsStats, PrintStream out) {
		out.print(v.getSequenceName()+"\t"+v.getFirst()+"\t"+v.getLast()+"\t"+v.getReference());
		for(DiversityStatistics stats:groupsStats) {
			if(stats == null) {
				out.print("\t0:0:0:0:0:0:0:0");
			} else {
				out.print("\t"+stats.getNumSamplesGenotyped());
				out.print(":"+fmt.format(stats.getExpectedHeterozygosity()));
				out.print(":"+fmt.format(stats.getObservedHeterozygosity()));
				out.print(":"+fmt.format(stats.getfStatistic()));
				out.print(":"+fmt.format(stats.getMaf()));
				out.print(":"+fmt.format(stats.getAlleleFrequencies()[0]));
				out.print(":"+fmt.format(stats.getChiSquareValue()));
				out.print(":"+fmt.format(stats.getChiSquarePValue()));
			}
		}
		out.println();
		
	}

	private List<CalledGenomicVariant> selectCalls(List<CalledGenomicVariant> allCalls, List<Integer> indexes) {
		List<CalledGenomicVariant> answer = new ArrayList<CalledGenomicVariant>();
		for(int i:indexes) answer.add(allCalls.get(i));
		return answer;
	}

	
}
