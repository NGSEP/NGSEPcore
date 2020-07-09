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
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.logging.Logger;

import ngsep.main.CommandsDescriptor;
import ngsep.main.ProgressNotifier;
import ngsep.main.io.ParseUtils;
import ngsep.variants.CalledGenomicVariant;

/**
 * Program to calculate LD statistics
 * @author Jorge Duitama
 *
 */
public class VCFLDCalculator {
	
	
	public static final int MODE_WINDOW = 0;
	public static final int MODE_SEQUENCE_NAMES = 1;
	public static final int MODE_ALL_PAIRS = 2;
	
	
	
	private Logger log = Logger.getLogger(VCFLDCalculator.class.getName());
	private ProgressNotifier progressNotifier=null;
	private int mode = MODE_WINDOW;
	
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
	
	public static void main(String[] args) throws Exception {
		VCFLDCalculator instance = new VCFLDCalculator();
		int i=CommandsDescriptor.getInstance().loadOptions(instance, args);
		boolean systemInput = "-".equals(args[i]);
		if(systemInput) {
			instance.run(System.in, System.out);
		} else {
			String filename = args[i];
			instance.run(filename, System.out);
		}

	}

	/**
	 * @return the mode
	 */
	public int getMode() {
		return mode;
	}
	/**
	 * @param mode the mode to set
	 */
	public void setMode(int mode) {
		this.mode = mode;
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
	
	public void run(VCFFileReader in, PrintStream out) {
		if(log!=null)in.setLog(log);
		
		List<VCFRecord> recordsInMemory = new LinkedList<>();
		in.setLoadMode(VCFFileReader.LOAD_MODE_MINIMAL);
		//TODO: Implement modes
		Iterator<VCFRecord> it = in.iterator();
		String lastSeqName = null;
		int n=0;
		while(it.hasNext()) {
			VCFRecord record = it.next();
			if(!record.getVariant().isBiallelic()) continue;
			if(!record.getSequenceName().equals(lastSeqName)) {
				if(recordsInMemory.size()>0) {
					calculateLDStatistics(recordsInMemory,out);
					recordsInMemory.clear();
				}
				lastSeqName = record.getSequenceName();
			}
			recordsInMemory.add(record);
			
			n++;
			if (progressNotifier!=null && n%1000==0) {
				int progress = n/1000;
				if (!progressNotifier.keepRunning(progress)) {
					out.flush();
					return;
				}
			}
		}
		//System.out.println("Loaded: "+recordsInMemory.size()+" variants");
		if(recordsInMemory.size()>0) calculateLDStatistics(recordsInMemory,out);
	}
	/**
	 * Calculates LD statistics for all pairs of records within the given list
	 * @param records to process
	 * @param out stream to write results
	 */
	public void calculateLDStatistics(List<VCFRecord> records, PrintStream out) {
		DecimalFormat fmt = ParseUtils.ENGLISHFMT_PROBABILITIES;
		int n = records.size();
		//Array to ensure constant lookup time
		VCFRecord [] recordsArray = records.toArray(new VCFRecord[0]);
		for(int i=0;i<n;i++) {
			for(int j=i+1;j<n;j++) {
				VCFRecord r1 = recordsArray[i];
				VCFRecord r2 = recordsArray[j];
				LDStatistics stats = calculateLDStatistics (r1, r2);
				System.out.print(r1.getSequenceName()+"\t"+r1.getFirst()+"\t"+r1.getLast()+"\t"+r2.getSequenceName()+"\t"+r2.getFirst()+"\t"+r2.getLast());
				System.out.println("\t"+(r2.getFirst()-r1.getFirst())+"\t"+stats.getSharedVariants()+"\t"+fmt.format(stats.getD())+"\t"+fmt.format(stats.getDPrime())+"\t"+fmt.format(stats.getR2()));
			}
		}
		
	}
	public LDStatistics calculateLDStatistics(VCFRecord record1, VCFRecord record2) {
		List<CalledGenomicVariant> calls1= record1.getCalls();
		List<CalledGenomicVariant> calls2= record2.getCalls();
		int n = calls1.size();
		//Frequency of alleles together
		double n00=0;
		//Individual frequencies of allele zero in shared sites
		double n01=0;
		double n02=0;
		
		int shared = 0;
		for(int i=0;i<n;i++) {
			CalledGenomicVariant call1 = calls1.get(i);
			CalledGenomicVariant call2 = calls2.get(i);
			if(call1.isUndecided() || call1.isHeterozygous()) continue;
			if(call2.isUndecided() || call2.isHeterozygous()) continue;
			shared++;
			if(call1.isHomozygousReference()) {
				n01++;
				if(call2.isHomozygousReference()) n00++;
			}
			if(call2.isHomozygousReference()) n02++;
		}
		if(shared == 0) return new LDStatistics(0, 0, 0, shared);
		double p00 = n00/shared;
		double p01 = n01/shared;
		double p02 = n02/shared;
		double d = p00-p01*p02;
		double dPrime;
		if(p01==0 || p02==0 || p01==1 || p02==1) dPrime = 0; 
		else if(d<0) dPrime = d/Math.min(p01*p02, (1-p01)*(1-p02));
		else dPrime = d/Math.min(p01*(1-p02), (1-p01)*p02);
		double r2 = d*d;
		if(p01==0 || p02==0 || p01==1 || p02==1) r2 = 0; 
		else r2/= (p01*p02*(1-p01)*(1-p02));
		return new LDStatistics(d, dPrime, r2, shared);
	}

}
class LDStatistics {
	private double d;
	private double dPrime;
	private double r2;
	private int sharedVariants;
	public LDStatistics(double d, double dPrime, double r2, int shared) {
		super();
		this.d = d;
		this.dPrime = dPrime;
		this.r2 = r2;
		this.sharedVariants = shared;
	}
	/**
	 * @return the d
	 */
	public double getD() {
		return d;
	}
	/**
	 * @return the dPrime
	 */
	public double getDPrime() {
		return dPrime;
	}
	/**
	 * @return the r2
	 */
	public double getR2() {
		return r2;
	}
	/**
	 * @return the sharedVariants
	 */
	public int getSharedVariants() {
		return sharedVariants;
	}
}
