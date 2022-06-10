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
package ngsep.transcriptome;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import ngsep.transcriptome.io.GFF3TranscriptomeHandler;

/**
 * 
 * @author Jorge Duitama
 *
 */

public class OrthogroupsFunctionalReport {

	public static void main(String[] args) throws Exception {
		String clustersFile = args[0];
		List<Transcriptome> transcriptomes = new ArrayList<>();
		GFF3TranscriptomeHandler handler = new GFF3TranscriptomeHandler();
		handler.setLoadTextAnnotations(true);
		for(int i=1;i<args.length;i++) {
			Transcriptome t =  handler.loadMap(args[i]);
			transcriptomes.add(t);
		}
		try (FileReader reader = new FileReader(clustersFile);
			 BufferedReader in = new BufferedReader(reader)) {
			
			String line = in.readLine();
			while(line!=null) {
				String [] items = line.split("\t");
				String clusterId = items[0];
				Set<String> annotations = new HashSet<>();
				List<Set<String>> genesPerTranscriptome = new ArrayList<>();
				for(int j=0;j<transcriptomes.size();j++) {
					genesPerTranscriptome.add(new HashSet<>());
				}
				for (int i=1;i<items.length;i++) {
					for(int j=0;j<transcriptomes.size();j++) {
						Transcriptome t = transcriptomes.get(j);
						Gene g = t.getGene(items[i]);
						if(g==null) continue;
						genesPerTranscriptome.get(j).add(items[i]);
						annotations.addAll(g.getTextFunctionalAnnotations());
					}
				}
				System.out.print(""+clusterId);
				for(Set<String> gpt:genesPerTranscriptome) System.out.print("\t"+gpt.size());
				for(Set<String> gpt:genesPerTranscriptome) System.out.print("\t"+gpt);
				System.out.print("\t"+annotations);
				System.out.println();
				line = in.readLine();
			}
		}

	}

}
