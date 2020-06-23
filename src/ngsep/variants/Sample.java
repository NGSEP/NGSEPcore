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
package ngsep.variants;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

public class Sample {
	private String id;
	private String group;
	private short normalPloidy = GenomicVariant.DEFAULT_PLOIDY;
	private Set<String> readGroups = new HashSet<>();
	
	
	public Sample(String id) {
		super();
		this.id = id;
	}
	public String getId() {
		return id;
	}
	public void setId(String id) {
		this.id = id;
	}
	public String getGroup() {
		return group;
	}
	public void setGroup(String group) {
		this.group = group;
	}
	
	
	public short getNormalPloidy() {
		return normalPloidy;
	}
	public void setNormalPloidy(short normalPloidy) {
		this.normalPloidy = normalPloidy;
	}
	
	/**
	 * @return the readGroups
	 */
	public Set<String> getReadGroups() {
		return readGroups;
	}
	/**
	 * @param readGroups the readGroups to set
	 */
	public void addReadGroup(String readGroup) {
		this.readGroups.add(readGroup);
	}
	/**
	 * Gets a list with the groups represented by the given list of samples
	 * @param samples to look for groups
	 * @return Groups associated with at least one sample
	 */
	public static List<String> getGroupsList(Collection<Sample> samples) {
		Set<String> groups = new TreeSet<String>();
		for(Sample s:samples) groups.add(s.getGroup());
		return new ArrayList<String>(groups);
	}
	/**
	 * Returns the list of sample ids that belong to any of the given groupIds
	 * @param samples List of samples
	 * @param groupIds Ids of the groups to look for
	 * @return String [] List of samples with a group id contained in groupIds
	 */
	public static Set<String> getSampleIds (List<Sample> samples, List<String> groupIds) {
		Set <String> groupsSet = new TreeSet<String>(groupIds);
		Set<String> answer = new TreeSet<String>();
		for(Sample s:samples) {
			if(s.getGroup()!=null && groupsSet.contains(s.getGroup())) {
				answer.add(s.getId());
			}
		}
		return answer;
	}
	/**
	 * Returns the list of sample ids that belong to any of the given groupIds
	 * @param samples List of samples
	 * @param groupIds Ids of the groups to look for
	 * @return String [] List of samples with a group id contained in groupIds
	 */
	public static Set<String> getSampleIds (List<Sample> samples, String [] groupIds) {
		return getSampleIds(samples, Arrays.asList(groupIds));
	}
	/**
	 * Returns a list of sample ids given the list of samples
	 * @param samples List with complete samples information
	 * @return List of sample ids in the same order as that of the input list
	 */
	public static List<String> getSampleIds (List<Sample> samples) {
		List<String> answer = new ArrayList<String>();
		for(Sample s:samples) {
			answer.add(s.getId());
		}
		return answer;
	}
	/**
	 * Returns the ids of the samples that belong to the given group
	 * @param samples List of samples to look for
	 * @param groupId Id of the group
	 * @return Set<String> List of sample ids associated with the group
	 */
	public static Set<String> getSampleIds(List<Sample> samples, String groupId) {
		List<String> groupList = new ArrayList<String>();
		groupList.add(groupId);
		return getSampleIds(samples, groupList);
	}
	
	public static Map<String,List<Sample>> indexSamplesByGroup(List<Sample> samples) {
		Map<String,List<Sample>> samplesDatabase = new TreeMap<String, List<Sample>>();
		for(Sample sample:samples) {
			String groupId = sample.getGroup(); 
			if(groupId==null) continue;
			List<Sample> samplesGroup = samplesDatabase.get(groupId);
			if(samplesGroup==null) {
				samplesGroup = new ArrayList<Sample>();
				samplesDatabase.put(groupId, samplesGroup);
			}
			samplesGroup.add(sample);
		}
		return samplesDatabase;
	}
	
	public static Map<String, List<Integer>> getGroupsWithSampleIdxs(Map<String,Sample> samplesMap, List<String> sampleIdsList) {
		Map<String, List<Integer>> groupSampleIndexes = new TreeMap<String, List<Integer>>();
		if(samplesMap==null) return groupSampleIndexes;
		for(int i=0;i<sampleIdsList.size();i++) {
			Sample sample = samplesMap.get(sampleIdsList.get(i));
			if(sample == null || sample.getGroup()==null) continue;
			List<Integer> idxsGroup = groupSampleIndexes.get(sample.getGroup());
			if(idxsGroup == null) {
				idxsGroup = new ArrayList<Integer>();
				groupSampleIndexes.put(sample.getGroup(), idxsGroup);
			}
			idxsGroup.add(i);
		}
		return groupSampleIndexes;
	}
	
}
