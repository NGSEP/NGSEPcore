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
package ngsep.sequencing;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * 
 * @author Jorge Duitama
 *
 */
public class SequencingLane {
	private String flowcell;
	private String number;
	private BarcodeMap barcodeMap = new BarcodeMap();
	private List<String> filesForward = new ArrayList<>();
	private List<String> filesReverse = new ArrayList<>();
	
	public SequencingLane(String flowcell, String number) {
		super();
		this.flowcell = flowcell;
		this.number = number;
	}
	/**
	 * @return the flowcell
	 */
	public String getFlowcell() {
		return flowcell;
	}
	
	/**
	 * @return the number
	 */
	public String getNumber() {
		return number;
	}
	/**
	 * 
	 * @return String a single unique key for this lane
	 */
	public String getUniqueKey() {
		return flowcell+"_"+number;
	}
	
	/**
	 * @return the barcodeMap
	 */
	public BarcodeMap getBarcodeMap() {
		return barcodeMap;
	}
	/**
	 * @param barcodeMap the barcodeMap to set
	 */
	public void setBarcodeMap(BarcodeMap barcodeMap) {
		this.barcodeMap = barcodeMap;
	}
	public void addSingleEndFile(String filename) {
		filesForward.add(filename);
		filesReverse.add(null);
	}
	public void addPairedEndFiles(String filename1, String filename2) {
		filesForward.add(filename1);
		filesReverse.add(filename2);
	}
	
	public void addSingleBarcode(String barcode, String sampleId) throws IOException {
		barcodeMap.addSingleBarcode (barcode, sampleId);
	}
	public void addDualBarcode(String barcode1, String barcode2, String sampleId) throws IOException {
		barcodeMap.addDualBarcode (barcode1, barcode2, sampleId);
	}
	
	public boolean equals(String flowcell, String number) {
		return flowcell.equals(this.flowcell) && number.equals(this.number);
	}
	/**
	 * @return the filesForward
	 */
	public List<String> getFilesForward() {
		return filesForward;
	}
	/**
	 * @return the filesReverse
	 */
	public List<String> getFilesReverse() {
		return filesReverse;
	}
	
	
}
