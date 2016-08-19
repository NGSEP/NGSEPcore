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
package ngsep.sequences;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.ListIterator;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;


/**
 * Special collection for qualified sequences which allows to compare sequence names given an order defined
 * by the index where each sequence is located. Useful to sort and compare genomic regions from different
 * sequences taking into account that the comparison between the names usually is not lexicographical
 * @author Jorge Duitama
 *
 */
public class QualifiedSequenceList implements List<QualifiedSequence> {
	private List<QualifiedSequence> sequences = new ArrayList<QualifiedSequence>();
	
	//Cache data to update
	private Map<String, Integer> sequenceIndexes = new HashMap<String, Integer>();
	private List<String> sequenceNames = new ArrayList<String>();
	private long totalLength = 0;
	private boolean sequenceDataUpdated = false;
	
	private boolean allowChanges = true;
	
	
	public QualifiedSequenceList() {
		
	}
	public QualifiedSequenceList (List<QualifiedSequence> sequences) {
		addAll(sequences);
	}
	
	private void updateSequenceData() {
		if(sequenceDataUpdated) return;
		List<String> seqNames = new ArrayList<String>();
		sequenceIndexes.clear();
		totalLength = 0;
		for(int i=0;i<sequences.size();i++) {
			QualifiedSequence s = sequences.get(i); 
			sequenceIndexes.put(s.getName(), i);
			seqNames.add(s.getName());
			totalLength += s.getLength();
		}
		sequenceNames = Collections.unmodifiableList(seqNames);
		sequenceDataUpdated = true;
	}
	
	public boolean isAllowChanges() {
		return allowChanges;
	}
	public void setAllowChanges(boolean allowChanges) {
		this.allowChanges = allowChanges;
	}
	public void add(int index, QualifiedSequence element) {
		validate(element);
		sequences.add(index, element);
		sequenceDataUpdated = false;
	}
	
	public boolean add(QualifiedSequence e) {
		if(e==null) return false;
		String name = e.getName();
		if(indexOf(name)>=0) return false;
		failIfChangesDisabled();
		sequences.add(e);
		sequenceDataUpdated = false;
		
		return true;
	}
	
	public boolean addAll(Collection<? extends QualifiedSequence> c) {
		boolean ret = false;
		for(QualifiedSequence s:c) {
			if(add(s)) ret = true;
		}
		return ret;
	}
	public boolean addAll(int index, Collection<? extends QualifiedSequence> c) {
		failIfChangesDisabled();
		for(QualifiedSequence e:c) {
			validate(e);
			//Added temporarily to fail if two input sequences have the same name
			sequenceIndexes.put(e.getName(), 0);
		}
		boolean ret = sequences.addAll(index,c);
		sequenceDataUpdated = false;
		return ret;
	}
	/**
	 * Creates a default qualified sequence with the given name and adds it to the list
	 * @param sequenceName Sequence name to add
	 * @return QualifiedSequence New sequence created with the given name or existing sequence
	 */
	public QualifiedSequence addOrLookupName (String sequenceName) {
		int index = indexOf(sequenceName); 
		if(index>=0) return sequences.get(index);
		QualifiedSequence seq = new QualifiedSequence(sequenceName);
		add(seq);
		return seq;
	}
	
	@Override
	public boolean remove(Object o) {
		int index = indexOf(o);
		if(index == -1) return false;
		remove(index);
		return true;
	}
	@Override
	public QualifiedSequence remove(int index) {
		failIfChangesDisabled();
		QualifiedSequence ret = sequences.remove(index);
		sequenceDataUpdated = false;
		return ret;
	}
	private void failIfChangesDisabled() {
		if (!allowChanges) throw new RuntimeException("Changes have been disabled for this collection");
	}
	@Override
	public boolean removeAll(Collection<?> c) {
		//Calculate indexes to remove
		Set<Integer> toRemove = new TreeSet<Integer>();
		for(Object o:c) {
			int nextIndex = indexOf(o);
			if(nextIndex>=0) toRemove.add(nextIndex);
		}
		if(toRemove.size()==0) return false;
		failIfChangesDisabled();
		List<Integer> toRemoveList = new ArrayList<Integer>(toRemove);
		//Parallel search to find indexes to retain
		List<Integer> toRetain = new ArrayList<Integer>();
		int j=0;
		int l = sequences.size();
		for(int i=0;i<l;i++) {
			int nextIdxRemove =l;
			if(j<toRemoveList.size()) nextIdxRemove = toRemoveList.get(j);
			if(i<nextIdxRemove) toRetain.add(i);
			else if(i==nextIdxRemove) j++;
			else throw new RuntimeException("Unsorted index to remove "+nextIdxRemove);
		}
		retainIndexes(toRetain);
		return true;
	}
	@Override
	public boolean retainAll(Collection<?> c) {
		Set<Integer> toRetain = new TreeSet<Integer>();
		for(Object o:c) {
			int nextIndex = indexOf(o);
			if(nextIndex>=0) toRetain.add(nextIndex);
		}
		if(sequences.size()==toRetain.size()) return false;
		failIfChangesDisabled();
		retainIndexes(toRetain);
		return true;
	}
	private void retainIndexes(Collection<Integer> toRetain) {
		List<QualifiedSequence> newList = new ArrayList<QualifiedSequence>();
		for(int i:toRetain) {
			newList.add(sequences.get(i));
		}
		sequences = newList;
		sequenceDataUpdated = false;
	}
	@Override
	public QualifiedSequence set(int index, QualifiedSequence element) {
		validate(element);
		failIfChangesDisabled();
		QualifiedSequence old = sequences.set(index, element);
		sequenceDataUpdated = false;
		return old;
	}

	private void validate(QualifiedSequence e) {
		if(e==null || e.getName()==null) throw new IllegalArgumentException("Null values not allowed in QualifiedSequenceList");
		if(sequenceIndexes.containsKey(e.getName())) throw new IllegalArgumentException("Sequence name "+e+ " already included in sequence list");
	}
	public void clear() {
		failIfChangesDisabled();
		sequences.clear();
		sequenceDataUpdated = false;
	}
	@Override
	public boolean contains(Object o) {
		return indexOf(o)>=0;
	}
	@Override
	public boolean containsAll(Collection<?> c) {
		for(Object o:c) {
			if(!contains(o)) return false;
		}
		return true;
	}
	@Override
	public QualifiedSequence get(int index) {
		updateSequenceData();
		return sequences.get(index);
	}
	/**
	 * Returns the sequence with the given name 
	 * @param s Name to look for
	 * @return QualifiedSequence Sequence with the given name. Null if it is not found 
	 */
	public QualifiedSequence get(String s) {
		int index = indexOf(s);
		if(index>=0) return sequences.get(index);
		return null;
	}
	private int lastIndex = -1;
	private Object lastQuery = null;
	@Override
	public synchronized int indexOf(Object o) {
		Integer ret = null;
		if(o==null) return -1;
		if(sequenceDataUpdated && (o==lastQuery || o.equals(lastQuery))) return lastIndex;
		updateSequenceData();
		if(o instanceof String) ret = sequenceIndexes.get(o);
		else if (o instanceof QualifiedSequence) ret = sequenceIndexes.get(((QualifiedSequence)o).getName());
		if(ret == null) return -1;
		lastQuery = o;
		lastIndex = ret;
		return ret;
	}
	@Override
	public boolean isEmpty() {
		return sequences.isEmpty();
	}
	@Override
	public Iterator<QualifiedSequence> iterator() {
		return sequences.iterator();
	}
	@Override
	public int lastIndexOf(Object o) {
		//Since no duplicates are allowed, this does the same as indexOf
		return indexOf(o);
	}
	@Override
	public ListIterator<QualifiedSequence> listIterator() {
		return sequences.listIterator();
	}
	@Override
	public ListIterator<QualifiedSequence> listIterator(int index) {
		return sequences.listIterator(index);
	}
	
	@Override
	public int size() {
		return sequences.size();
	}
	@Override
	public List<QualifiedSequence> subList(int fromIndex, int toIndex) {
		return sequences.subList(fromIndex, toIndex);
	}
	@Override
	public Object[] toArray() {
		return sequences.toArray();
	}
	@Override
	public <T> T[] toArray(T[] a) {
		return sequences.toArray(a);
	}
	public long getTotalLength() {
		updateSequenceData();
		return totalLength;
	}
	public void setTotalLength(long totalLength) {
		this.totalLength = totalLength;
	}
	
	public List<String> getNamesStringList () {
		updateSequenceData();
		return sequenceNames;
	}

}
