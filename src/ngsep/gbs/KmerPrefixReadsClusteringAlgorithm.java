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
package ngsep.gbs;

import java.io.IOException;
import java.util.Iterator;
import java.util.List;

import ngsep.main.CommandsDescriptor;
import ngsep.main.ProgressNotifier;
import ngsep.sequences.DNASequence;
import ngsep.sequences.DNAShortKmer;
import ngsep.sequences.DNAShortKmerClusterMap;
import ngsep.sequences.RawRead;
import ngsep.sequences.io.FastqFileReader;
import ngsep.sequencing.ReadsDemultiplex;
import ngsep.sequencing.SequencingLane;

/**
 * @author Jorge Gomez
 * @author Jorge Duitama
 */
public class KmerPrefixReadsClusteringAlgorithm {

	private static int KMER_LENGTH = 31;
	private static int PREFIX = 10;
	private ReadsDemultiplex demultiplex = new ReadsDemultiplex();
	
	
	public static void main(String[] args) throws Exception {
		KmerPrefixReadsClusteringAlgorithm instance = new KmerPrefixReadsClusteringAlgorithm();
		int i = CommandsDescriptor.getInstance().loadOptions(instance, args);
		
		String indexFile = args[i++];
		instance.run(indexFile);
	}

	
	/**
	 * @param progressNotifier
	 * @see ngsep.sequencing.ReadsDemultiplex#setProgressNotifier(ngsep.main.ProgressNotifier)
	 */
	public void setProgressNotifier(ProgressNotifier progressNotifier) {
		demultiplex.setProgressNotifier(progressNotifier);
	}


	/**
	 * @return
	 * @see ngsep.sequencing.ReadsDemultiplex#getProgressNotifier()
	 */
	public ProgressNotifier getProgressNotifier() {
		return demultiplex.getProgressNotifier();
	}


	/**
	 * @return
	 * @see ngsep.sequencing.ReadsDemultiplex#getOutDirectory()
	 */
	public String getOutDirectory() {
		return demultiplex.getOutDirectory();
	}


	/**
	 * @param outDirectory
	 * @see ngsep.sequencing.ReadsDemultiplex#setOutDirectory(java.lang.String)
	 */
	public void setOutDirectory(String outDirectory) {
		demultiplex.setOutDirectory(outDirectory);
	}


	/**
	 * @return
	 * @see ngsep.sequencing.ReadsDemultiplex#getLaneFilesDescriptor()
	 */
	public String getLaneFilesDescriptor() {
		return demultiplex.getLaneFilesDescriptor();
	}


	/**
	 * @param laneFilesDescriptor
	 * @see ngsep.sequencing.ReadsDemultiplex#setLaneFilesDescriptor(java.lang.String)
	 */
	public void setLaneFilesDescriptor(String laneFilesDescriptor) {
		demultiplex.setLaneFilesDescriptor(laneFilesDescriptor);
	}


	/**
	 * @return
	 * @see ngsep.sequencing.ReadsDemultiplex#getFlowcell()
	 */
	public String getFlowcell() {
		return demultiplex.getFlowcell();
	}


	/**
	 * @param flowcell
	 * @see ngsep.sequencing.ReadsDemultiplex#setFlowcell(java.lang.String)
	 */
	public void setFlowcell(String flowcell) {
		demultiplex.setFlowcell(flowcell);
	}


	/**
	 * @return
	 * @see ngsep.sequencing.ReadsDemultiplex#getLane()
	 */
	public String getLane() {
		return demultiplex.getLane();
	}


	/**
	 * @param lane
	 * @see ngsep.sequencing.ReadsDemultiplex#setLane(java.lang.String)
	 */
	public void setLane(String lane) {
		demultiplex.setLane(lane);
	}


	/**
	 * @param prefix
	 * @see ngsep.sequencing.ReadsDemultiplex#setPrefix(java.lang.String)
	 */
	public void setPrefix(String prefix) {
		demultiplex.setPrefix(prefix);
	}


	/**
	 * @return
	 * @see ngsep.sequencing.ReadsDemultiplex#getMinReadLength()
	 */
	public int getMinReadLength() {
		return demultiplex.getMinReadLength();
	}


	/**
	 * @param minReadLength
	 * @see ngsep.sequencing.ReadsDemultiplex#setMinReadLength(int)
	 */
	public void setMinReadLength(int minReadLength) {
		demultiplex.setMinReadLength(minReadLength);
	}


	/**
	 * @param minReadLength
	 * @see ngsep.sequencing.ReadsDemultiplex#setMinReadLength(java.lang.Integer)
	 */
	public void setMinReadLength(Integer minReadLength) {
		demultiplex.setMinReadLength(minReadLength);
	}


	/**
	 * @return
	 * @see ngsep.sequencing.ReadsDemultiplex#isDualBarcode()
	 */
	public boolean isDualBarcode() {
		return demultiplex.isDualBarcode();
	}


	/**
	 * @param dualBarcode
	 * @see ngsep.sequencing.ReadsDemultiplex#setDualBarcode(boolean)
	 */
	public void setDualBarcode(boolean dualBarcode) {
		demultiplex.setDualBarcode(dualBarcode);
	}


	/**
	 * @param dualBarcode
	 * @see ngsep.sequencing.ReadsDemultiplex#setDualBarcode(java.lang.Boolean)
	 */
	public void setDualBarcode(Boolean dualBarcode) {
		demultiplex.setDualBarcode(dualBarcode);
	}


	/**
	 * @return
	 * @see ngsep.sequencing.ReadsDemultiplex#getTrimSequence()
	 */
	public String getTrimSequence() {
		return demultiplex.getTrimSequence();
	}


	/**
	 * @param trimSequence
	 * @see ngsep.sequencing.ReadsDemultiplex#setTrimSequence(java.lang.String)
	 */
	public void setTrimSequence(String trimSequence) {
		demultiplex.setTrimSequence(trimSequence);
	}


	/**
	 * @return
	 * @see ngsep.sequencing.ReadsDemultiplex#isUncompressedOutput()
	 */
	public boolean isUncompressedOutput() {
		return demultiplex.isUncompressedOutput();
	}


	/**
	 * @param uncompressedOutput
	 * @see ngsep.sequencing.ReadsDemultiplex#setUncompressedOutput(boolean)
	 */
	public void setUncompressedOutput(boolean uncompressedOutput) {
		demultiplex.setUncompressedOutput(uncompressedOutput);
	}


	/**
	 * @param uncompressedOutput
	 * @see ngsep.sequencing.ReadsDemultiplex#setUncompressedOutput(java.lang.Boolean)
	 */
	public void setUncompressedOutput(Boolean uncompressedOutput) {
		demultiplex.setUncompressedOutput(uncompressedOutput);
	}


	public void run(String indexFile) throws IOException {
		demultiplex.loadIndexAndLaneFiles(indexFile);
		clusterReads();
		
	}

	public void clusterReads() throws IOException {
		List<SequencingLane> lanes = demultiplex.getLanes();
		DNAShortKmerClusterMap map = new DNAShortKmerClusterMap();
		for(SequencingLane lane:lanes) {
			List<String> filesForward = lane.getFilesForward();
			for(int i=0;i<filesForward.size();i++) {
				String filename1 = filesForward.get(i);
				processFile(map, filename1);
			}
		}
		
	}

	private void processFile(DNAShortKmerClusterMap kmerMap, String filename) throws IOException {
		int readCount = 0;
		try (FastqFileReader openFile = new FastqFileReader(filename);) {
			Iterator<RawRead> reader = openFile.iterator();
			while(reader.hasNext()) {
				RawRead read = reader.next();
				String s = read.getSequenceString();
				//System.out.println(s);
				if(DNASequence.isDNA(s)) {
					kmerMap.addOcurrance(new DNAShortKmer(s.substring(PREFIX,PREFIX + KMER_LENGTH)));
					readCount++;
				}
			}
		}
		System.out.println("Processed a total of " + readCount + " reads for file: "+filename);
		
	}
	
}
