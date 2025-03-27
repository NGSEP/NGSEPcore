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
package ngsep.transposons;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.logging.Logger;

import ngsep.main.CommandsDescriptor;
import ngsep.main.ProgressNotifier;
import ngsep.sequences.DNAMaskedSequence;
import ngsep.transposons.io.TransposableElementLibraryHandler;

/**
 * @author Jorge Duitama
 * @author Ana Castellanos
 */
public class TransposableElementLibraryFilter {

	// Logging and progress
	private Logger log = Logger.getLogger(TransposableElementLibraryFilter.class.getName());
	private ProgressNotifier progressNotifier = null;

	public static final int DEF_NUM_THREADS = 1;

	// Parameters
	private String inputFile = null;
	private String outputFile = null;

	private int numThreads = DEF_NUM_THREADS;

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

	public String getOutputFile() {
		return outputFile;
	}

	public void setOutputFile(String outputFile) {
		this.outputFile = outputFile;
	}

	public int getNumThreads() {
		return numThreads;
	}

	public void setNumThreads(int numThreads) {
		this.numThreads = numThreads;
	}

	public void setNumThreads(String value) {
		this.setNumThreads(Integer.parseInt(value));
	}

	public static void main(String[] args) throws Exception {
		TransposableElementLibraryFilter instance = new TransposableElementLibraryFilter();
		CommandsDescriptor.getInstance().loadOptions(instance, args);
		instance.run();
	}

	public void run() throws IOException, InterruptedException {
		long time = System.currentTimeMillis();
		logParameters();
		if (inputFile == null)
			throw new IOException("The input library is a required parameter");
		if (outputFile == null)
			throw new IOException("The output file is a required parameter");
		TransposableElementLibraryHandler handler = new TransposableElementLibraryHandler();
		List<TransposableElement> tes = handler.load(inputFile);
		List<TransposableElement> filteredTEs = filterTEs(tes);
		handler.save(filteredTEs, outputFile);
		double seconds = (System.currentTimeMillis() - time);
		seconds /= 1000;
		log.info("Process finished in " + seconds + " seconds");
	}

	public void logParameters() {
		ByteArrayOutputStream os = new ByteArrayOutputStream();
		PrintStream out = new PrintStream(os);
		out.println("Input file:" + inputFile);
		out.println("Output file:" + outputFile);
		out.println("Number of threads: " + numThreads);
		log.info(os.toString());
	}

	public List<TransposableElement> filterTEs(List<TransposableElement> tes) {
		int n = tes.size();
		List<TransposableElement> filteredTEsWithNulls = new ArrayList<TransposableElement>(n);
		for (int i = 0; i < n; i++)
			filteredTEsWithNulls.add(null);
		try (ThreadPoolExecutor pool = new ThreadPoolExecutor(numThreads, numThreads, n, TimeUnit.SECONDS,
				new LinkedBlockingQueue<Runnable>())) {
			for (int i = 0; i < n; i++) {
				TransposableElement te = tes.get(i);
				final int x = i;
				pool.execute(() -> filteredTEsWithNulls.set(x, verifyTransposon(te)));
			}
			pool.shutdown();
			if (!pool.isShutdown()) {
				throw new RuntimeException("The ThreadPoolExecutor was not shutdown after an await termination call");
			}
		}
		List<TransposableElement> answer = new ArrayList<TransposableElement>();
		for (TransposableElement te : filteredTEsWithNulls) {
			if (te != null)
				answer.add(te);
		}
		return answer;
	}

	private TransposableElement verifyTransposon(TransposableElement te) {
		if(!te.verifyEnds()) {
			log.info("Terminal repeats at ends of TE: "+te.getId()+" could not be found.");
			return null;
		}
		HMMTransposonDomainsFinder finder = new HMMTransposonDomainsFinder();
		finder.loadHMMsFromClasspath();
		TransposableElementFamily givenFamily = te.getFamily();
		finder.assignFamily(te);
		TransposableElementFamily newFamily = te.getFamily();
		if(newFamily==null) {
			log.info("Family of TE: "+te.getId()+" could not be inferred.");
			return null;
		}
		if(givenFamily!=null) {
			if(!givenFamily.getOrder().equals(newFamily.getOrder())) {
				log.info("Calculated order of TE: "+te.getId()+" does not match given order. Given: "+givenFamily.getOrder()+" calculated: "+newFamily.getOrder());
				return null;
			}
			if(givenFamily.getId()!=TransposableElementFamily.UNKNOWN && newFamily.getId()!=TransposableElementFamily.UNKNOWN && !givenFamily.getId().equals(newFamily.getId())) {
				log.info("Calculated family of TE: "+te.getId()+" does not match given family. Given: "+givenFamily.getId()+" calculated: "+newFamily.getId());
				return null;
			}
		} else {
			te.modifyIdFromFamily();
		}
		TransposableElement answer = te;
		if(te.getDomainAlns()!=null && te.getDomainAlns().get(0).isReverse()) {
			TransposableElement revElem = new TransposableElement(te.getId(), DNAMaskedSequence.getReverseComplement(te.getSequence()));
			revElem.setFamily(newFamily);
			revElem.setDomainAlns(te.getDomainAlns());
			answer = revElem;
		}
		
		return answer;
	}
}
