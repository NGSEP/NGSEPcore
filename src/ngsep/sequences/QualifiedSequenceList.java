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

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
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
public class QualifiedSequenceList implements List<QualifiedSequence>, Serializable {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private List<QualifiedSequence> sequences = new ArrayList<QualifiedSequence>();
	private Map<String, Integer> sequenceIndexesMap = new HashMap<String, Integer>();
	
	private long totalLength = 0;
	
	private boolean allowChanges = true;
	
	
	public QualifiedSequenceList() {
		
	}
	public QualifiedSequenceList (List<QualifiedSequence> sequences) {
		addAll(sequences);
	}
	
	public boolean isAllowChanges() {
		return allowChanges;
	}
	public void setAllowChanges(boolean allowChanges) {
		this.allowChanges = allowChanges;
	}
	public void add(int index, QualifiedSequence element) {
		Integer indexPresent = validate(element);
		if(indexPresent!=null) throw new IllegalArgumentException("Can not add sequence in index "+index+". A sequence with name: "+element.getName()+" is already present in index "+indexPresent);
		sequences.add(index, element);
		updateSequenceIndexesMap();
	}
	private void updateSequenceIndexesMap() {
		sequenceIndexesMap.clear();
		totalLength = 0;
		for(int i=0;i<sequences.size();i++) {
			QualifiedSequence s = sequences.get(i); 
			sequenceIndexesMap.put(s.getName(), i);
			totalLength += s.getLength();
		}
	}
	
	public boolean add(QualifiedSequence e) {
		Integer index = validate(e);
		if(index!=null) return false;
		failIfChangesDisabled();
		//Attributes updated directly to make this operation in (amortized) constant time
		sequences.add(e);
		sequenceIndexesMap.put(e.getName(), sequences.size()-1);
		totalLength+=e.getLength();
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
			if(validate(e)!=null) throw new IllegalArgumentException("A sequence with name: "+e.getName()+" is already present in the list");
			//Added temporarily to fail if two input sequences have the same name
			sequenceIndexesMap.put(e.getName(), 0);
		}
		boolean ret = sequences.addAll(index,c);
		updateSequenceIndexesMap();
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
		if(ret!=null) updateSequenceIndexesMap();
		return ret;
	}
	private void failIfChangesDisabled() {
		if (!allowChanges) throw new RuntimeException("Changes have been disabled for this collection");
	}
	@Override
	public boolean removeAll(Collection<?> c) {
		//Calculate indexes to remove
		Set<Integer> toRemove = new HashSet<Integer>();
		for(Object o:c) {
			int nextIndex = indexOf(o);
			if(nextIndex>=0) toRemove.add(nextIndex);
		}
		if(toRemove.size()==0) return false;
		failIfChangesDisabled();
		List<Integer> toRetain = new ArrayList<Integer>();
		int l = sequences.size();
		for(int i=0;i<l;i++) {
			if(!toRemove.contains(i)) toRetain.add(i);
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
		updateSequenceIndexesMap();
	}
	@Override
	public QualifiedSequence set(int index, QualifiedSequence element) {
		failIfChangesDisabled();
		Integer indexPresent = validate(element);
		if(indexPresent!=null && indexPresent!=index) {
			throw new IllegalArgumentException("Can not set sequence in index "+index+" A sequence with name: "+element.getName()+" is already present in index "+indexPresent);
		}
		QualifiedSequence old = sequences.set(index, element);
		updateSequenceIndexesMap();
		return old;
	}

	/**
	 * Throws null pointer exceptions if the object has null elements. Returns true if the name of the given sequence is already in the list
	 * @param e Qualified sequence to validate
	 * @return Integer index where a sequence with the given name is present or null if the sequence is not present
	 */
	private Integer validate(QualifiedSequence e) {
		if(e==null) throw new NullPointerException("Null sequences not allowed in QualifiedSequenceList");
		if(e.getName()==null) throw new NullPointerException("Sequences with null names not allowed in QualifiedSequenceList");
		return sequenceIndexesMap.get(e.getName());
	}
	public void clear() {
		failIfChangesDisabled();
		sequences.clear();
		updateSequenceIndexesMap();
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
		if(o==lastQuery || o.equals(lastQuery)) return lastIndex;
		if(o instanceof String) ret = sequenceIndexesMap.get(o);
		else if (o instanceof QualifiedSequence) ret = sequenceIndexesMap.get(((QualifiedSequence)o).getName());
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
		return totalLength;
	}
	public void setTotalLength(long totalLength) {
		this.totalLength = totalLength;
	}
	
	public List<String> getNamesStringList () {
		List<String> answer = new ArrayList<>();
		for(QualifiedSequence seq:sequences) answer.add(seq.getName());
		return answer;
	}
	
	public List<? extends CharSequence> getSequencesDataList () {
		List<CharSequence> answer = new ArrayList<>();
		for(QualifiedSequence seq:sequences) answer.add(seq.getCharacters());
		return answer;
	}

}
