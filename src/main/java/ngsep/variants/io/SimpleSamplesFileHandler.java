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
package ngsep.variants.io;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import ngsep.sequences.QualifiedSequenceList;
import ngsep.variants.Sample;

public class SimpleSamplesFileHandler {
	public List<Sample> loadSamples (String filename) throws IOException {
		List<Sample> samples = new ArrayList<Sample>();
		QualifiedSequenceList groupIds = new QualifiedSequenceList();
		FileInputStream fis = null;
		BufferedReader in = null;
		try {
			fis = new FileInputStream(filename);
			in = new BufferedReader(new InputStreamReader(fis));
			String line=in.readLine();
			while (line!=null) {
				String [] items = line.split("\t");
				Sample s = new Sample(items[0]);
				if(items.length>1) {
					String groupId = groupIds.addOrLookupName(items[1]).getName();
					s.setGroup(groupId);
				}
				//if(items.length>2) s.setNormalPloidy(Integer.parseInt(items[2]));
				samples.add(s);
				line=in.readLine();
			}
		} finally {
			if (in!=null) in.close();
			if (fis!=null) fis.close();
		}
		
		return samples;
	}
	public Map<String, Sample> loadSamplesAsMap(String filename)  throws IOException {
		List<Sample> samples = loadSamples(filename);
		Map<String, Sample> answer = new TreeMap<String, Sample>();
		for(Sample s:samples) {
			answer.put(s.getId(), s);
		}
		return answer;
	}
	public void saveSamples(List<Sample> samples, PrintStream out) {
		for(Sample s: samples) {
			out.print(s.getId());
			if(s.getGroup()!=null) {
				out.print("\t"+s.getGroup());
			}
			out.println();
		}
		out.flush();
	}
	public Set<String> fillGroups(String filename, List<Sample> samples) throws IOException {
		Map<String,Sample> samplesWithPop = loadSamplesAsMap(filename);
		Set<String> groups = new TreeSet<String>();
		for(Sample s:samples) {
			Sample sampleWithPop = samplesWithPop.get(s.getId());
			if(sampleWithPop!=null && sampleWithPop.getGroup()!=null) {
				String group = sampleWithPop.getGroup(); 
				s.setGroup(group);
				groups.add(group);
			}
		}
		return groups;
	}
}
