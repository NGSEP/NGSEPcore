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
package ngsep.sequencing.io;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

import ngsep.sequences.DNASequence;
import ngsep.sequencing.SequencingLane;

/**
 * @author Jorge Duitama
 */
public class BarcodesFileLoader {

	private Logger log = Logger.getLogger(BarcodesFileLoader.class.getName());
	public List<SequencingLane> loadMultipleLanes(String indexFile, boolean dualBarcode) throws IOException {
		List<SequencingLane> answer = new ArrayList<>();
		try (FileInputStream fis = new FileInputStream(indexFile);
			BufferedReader in = new BufferedReader(new InputStreamReader(fis));) {
			String line = in.readLine();
			//Skip header
			log.info("Loading index file: "+indexFile);
			line = in.readLine();
			SequencingLane lane = null;
			for (int i=1;line != null;i++) {
				String[] items = line.split("\t");
				if(lane==null || !lane.equals(items[0],items[1])) {
					lane = new SequencingLane(items[0], items[1]);
					answer.add(lane);
					log.info("Loading barcode map for flowcell: "+lane.getFlowcell()+" lane: "+lane.getNumber());
				}
				processLine(lane, i, items, dualBarcode);
				line = in.readLine();
			}
		}
		return answer;
	}
	public void loadSingleLane(String indexFile, boolean dualBarcode, SequencingLane lane) throws IOException {
		try (FileInputStream fis = new FileInputStream(indexFile);
			BufferedReader in = new BufferedReader(new InputStreamReader(fis));) {
			String line = in.readLine();
			//Skip header
			log.info("Loading index file: "+indexFile);
			line = in.readLine();
			for (int i=1;line != null;i++) {
				String[] items = line.split("\t");
				if(lane.equals(items[0],items[1])) {
					processLine(lane, i, items, dualBarcode);
				}
				line = in.readLine();
			}
		}
	}
	private void processLine(SequencingLane lane, int i, String[] items, boolean dualBarcode) throws IOException {
		if (dualBarcode) {
			log.info("Barcode 1: "+items[2]+" barcode 2: "+items[3]+" sample: "+items[4]);
			if(!DNASequence.isDNA(items[2])) throw new IOException("Barcode: "+items[2]+" at line "+i+" is not a DNA sequence");
			if(!DNASequence.isDNA(items[3])) throw new IOException("Barcode: "+items[3]+" at line "+i+" is not a DNA sequence");
			lane.addDualBarcode(items[2],items[3], items[4]);
		} else {
			log.info("Barcode: "+items[2]+" sample: "+items[3]);
			if(!DNASequence.isDNA(items[2])) throw new IOException("Barcode: "+items[2]+" at line "+i+" is not a DNA sequence");
			lane.addSingleBarcode(items[2], items[3]);
		}
	}
}
