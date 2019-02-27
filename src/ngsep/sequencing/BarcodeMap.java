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
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import ngsep.sequences.DNASequence;

/**
 * @author Jorge Duitama
 */
public class BarcodeMap {

	private Map<String, String> barcodeSampleMap = new HashMap<>();
	private Set<String> barcodes1 = new HashSet<>();
	private Set<String> barcodes2 = new HashSet<>();
	
	//Sorted list of barcodes
	private String [] sortedBarcodes;
	private int [] [] ntPairFirsts = new int[4][4];
	private int [] [] ntPairLasts = new int[4][4];
	
	public void addSingleBarcode(String barcode, String sampleId) throws IOException {
		String oldSample = barcodeSampleMap.get(barcode);
		if(oldSample!=null && !oldSample.equals(sampleId)) throw new IOException("Barcode "+barcode+" already registered with a sample "+oldSample+" different from: "+sampleId);
		barcodeSampleMap.put(barcode, sampleId);
		barcodes1.add(barcode);
		sortedBarcodes=null;
	}
	public void addDualBarcode(String barcode1, String barcode2, String sampleId) throws IOException {
		String dualKey = buildDualKey(barcode1,barcode2);
		String oldSample = barcodeSampleMap.get(dualKey);
		if(oldSample!=null && !oldSample.equals(sampleId)) throw new IOException("Barcode pair"+barcode1+" - "+barcode2+" already registered with a sample "+oldSample+" different from: "+sampleId);
		barcodeSampleMap.put(dualKey, sampleId);
		barcodes1.add(barcode1);
		barcodes2.add(barcode2);
		sortedBarcodes=null;
	}
	private String buildDualKey(String barcode1, String barcode2) {
		return barcode1+"_"+barcode2;
	}
	public String [] getSampleIdByRead (String read) {
		String [] answer = new String [2];
		answer[1] = findBarcode(read,(byte)1);
		if(answer[1]==null) return null;
		answer[0] = barcodeSampleMap.get(answer[1]);
		return answer;
	}
	/**
	 * 
	 * @param read1
	 * @param read2
	 * @return
	 */
	public String [] getSampleIdByReadPair (String read1, String read2) {
		String [] answer = new String [3];
		answer[1] = findBarcode(read1,(byte)1);
		answer[2] = findBarcode(read2,(byte)2);
		if(answer[1]!=null && answer[2]!=null) {
			answer[0] = barcodeSampleMap.get(buildDualKey(answer[1], answer[2]));
		} else {
			return null;
		}
		return answer;
	}
	
	private void initBarcodeSortedList() {
		Set<String> allBarcodes = new TreeSet<>();
		allBarcodes.addAll(barcodes1);
		allBarcodes.addAll(barcodes2);
		
		sortedBarcodes = allBarcodes.toArray(new String[0]);
		for(int i=0;i<ntPairFirsts.length;i++) {
			Arrays.fill(ntPairFirsts[i], -1);
			Arrays.fill(ntPairLasts[i], -1);
		}
		String alphabet = DNASequence.BASES_STRING;
		for (int k=0;k<sortedBarcodes.length;k++) {
			String barcode = sortedBarcodes[k];
			int i = alphabet.indexOf(barcode.charAt(0));
			int j = alphabet.indexOf(barcode.charAt(1));
			if(ntPairFirsts[i][j]==-1)ntPairFirsts[i][j] = k;
			ntPairLasts[i][j] = k;
		}
	}
	/**
	 * Finds the barcode matching the given sequence
	 * @param readSeq Sequence to barcode
	 * @return String barcode matching the given sequence. Null if the barcode was not found
	 */
	private String findBarcode(String readSeq, byte readIndex) {
		if(sortedBarcodes==null) initBarcodeSortedList();
		String alphabet = DNASequence.BASES_STRING;
		int firstIndex = alphabet.indexOf(readSeq.charAt(0));
		int secondIndex = alphabet.indexOf(readSeq.charAt(1));
		if(firstIndex <0 || secondIndex < 0 ) return null;
		int first = ntPairFirsts[firstIndex][secondIndex];
		int last = ntPairLasts[firstIndex][secondIndex];
		if(first==-1) return null;
		String selectedBarcode = null;
		for(int i=first;i<=last;i++) {
			String barcode = sortedBarcodes[i];
			Set<String> barcodesReadIdx = barcodes1;
			if(readIndex==2) barcodesReadIdx = barcodes2;
			
			boolean barcodeFound = readSeq.startsWith(barcode) && barcodesReadIdx.contains(barcode);
			//The last condition is to taking into account cases where one barcode is substring of another
			//and pick the longest matching barcode
			if(barcodeFound && (selectedBarcode==null || selectedBarcode.length()<barcode.length()) ) {
				selectedBarcode = barcode;
			}
		}
		return selectedBarcode;
	}
}
