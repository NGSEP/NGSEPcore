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
package ngsep.benchmark;
import java.io.BufferedReader;
import java.io.FileReader;

import ngsep.sequences.DNASequence;

/**
 * Simple script to add the AD field to a NGSEP VCF file. Needed for comparyson with SynDip
 * @author Jorge Duitama
 *
 */
public class VCFAddADGTField {

	public static void main(String[] args) throws Exception {
		try (FileReader reader = new FileReader(args[0]); 
			 BufferedReader in = new BufferedReader(reader)) {
			String line = in.readLine();
			while (line != null) {
				if(line.charAt(0)=='#') {
					if(line.charAt(1)!='#') {
						System.out.println("##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">");
					}
					System.out.println(line);
				}
				else {
					String [] items = line.split("\t");
					System.out.print(items[0]);
					for(int i=1;i<8;i++) System.out.print("\t"+items[i]);
					System.out.print("\t"+items[8]+":AD\t"+items[9]);
					String [] items2 = items[9].split(":");
					if(items[8].contains("ADP")) {
						System.out.print(":"+items2[4]);
					} else if(items[8].contains("BSDP")) {
						String [] items3 = items2[4].split(",");
						int iR = DNASequence.BASES_STRING.indexOf(items[3].charAt(0));
						System.out.print(":"+items3[iR]);
						String [] items4 = items[4].split(",");
						for(int j=0;j<items4.length;j++) {
							int iA = DNASequence.BASES_STRING.indexOf(items4[j].charAt(0));
							System.out.print(","+items3[iA]);
						}
					}
					System.out.println();
				}
				line = in.readLine();
			}
		}

	}

}
