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
package ngsep.genome;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;

public class GenomicRegionSortedCollection<T extends GenomicRegion> implements Collection<T> {
	private QualifiedSequenceList sequences = new QualifiedSequenceList();
	private Map<Integer,List<T>> regionsMap = new HashMap<Integer, List<T>>();
	private Map<Integer,List<T>> longRegionsMap = new HashMap<Integer, List<T>>();
	private int size = 0;
	private int numSpanningLong = 15;
	private boolean sorted = true;
	
	
	public GenomicRegionSortedCollection () {
		
	}
	
	public GenomicRegionSortedCollection (Collection<T> regions) {
		addAll(regions);
	}
	public GenomicRegionSortedCollection (QualifiedSequenceList sequences) {
		this();
		this.sequences.addAll(sequences);
		for(int i=0;i<this.sequences.size();i++) {	 
			regionsMap.put(i, new ArrayList<T>());
			longRegionsMap.put(i, new ArrayList<T>());
		}
	}
	@Override
	public boolean add(T e) {
		if (e==null) return false;
		int index = sequences.indexOf(e.getSequenceName());
		if(index < 0) { 
			QualifiedSequence seq;
			try {
				seq = sequences.addOrLookupName(e.getSequenceName());
			} catch (RuntimeException ex ) {
				throw new IllegalArgumentException("Can not add new genomic region at "+e.getSequenceName()+":"+e.getFirst()+"-"+e.getLast()+" Invalid sequence name",ex);
			}
			index = sequences.indexOf(seq.getName());
			regionsMap.put(index, new ArrayList<T>());
			longRegionsMap.put(index, new ArrayList<T>());
		}
		List<T> regions = regionsMap.get(index);
		regions.add(e);
		size++;
		sorted = false;
		return true;
	}
	@Override
	public boolean addAll(Collection<? extends T> c) {
		boolean changed = false;
		for(T gr:c) {
			if(add(gr)) {
				changed = true;
			}
		}
		sorted =false;
		return changed;
	}
	
	@Override
	public void clear() {
		for(int index:regionsMap.keySet()) {
			regionsMap.get(index).clear();
			longRegionsMap.get(index).clear();
		}
		size = 0;
		sorted = true;
	}
	private static int indexOf(List<? extends GenomicRegion> regions, GenomicRegion gr ) {
		int index = Collections.binarySearch(regions, gr,GenomicRegionPositionComparator.getInstance());
		if(index<0) return -1;
		//Look for the specific object backwards
		for(int i=index;i>=0;i--) {
			GenomicRegion gr2 = regions.get(i);
			if(gr.equals(gr2)) return i;
			else if (gr.getFirst()!=gr2.getFirst()) break;
		}
		//Look for the specific object forward
		for(int i=index+1;i<regions.size();i++) {
			GenomicRegion gr2 = regions.get(i);
			if(gr.equals(gr2)) return i;
			else if (gr.getFirst()!=gr2.getFirst()) break;
		}
		return -1;
	}
	@Override
	public boolean remove(Object o) {
		sort();
		if(o==null || !(o instanceof GenomicRegion)) return false;
		GenomicRegion gr = (GenomicRegion)o;
		int sequenceIndex = sequences.indexOf(gr.getSequenceName());
		if(sequenceIndex<0) return false;
		List<T> regions = regionsMap.get(sequenceIndex);
		int index = GenomicRegionSortedCollection.indexOf (regions,gr);
		if(index<0) return false;
		regions.remove(index);
		size--;
		List<T> longRegions = longRegionsMap.get(sequenceIndex);
		index = GenomicRegionSortedCollection.indexOf (longRegions,gr);
		if(index>=0) longRegions.remove(index);
		return true;
	}
	@Override
	public boolean removeAll(Collection<?> c) {
		boolean changed = false;
		for(Object o:c) {
			if(remove(o)) {
				changed = true;
			}
		}
		return changed;
	}
	/**
	 * Removes the first n elements of the Collection
	 * @param n Number of elements to remove
	 * @return boolean true if the list is changed, false otherwise
	 */
	public boolean removeFirst(int n) {
		sort();
		boolean changed = false;
		int remaining = n;
		for(int i=0;i<sequences.size() && remaining > 0;i++) {
			List<T> regions = regionsMap.get(i);
			List<T> longRegions = longRegionsMap.get(i);
			int nSeq = regions.size();
			if(nSeq<=remaining) {
				regions.clear();
				longRegions.clear();
				remaining -= nSeq;
				size -= nSeq;
			} else {
				List<T> newRegions = new ArrayList<T>();
				for(int j=remaining;j<nSeq;j++) {
					newRegions.add(regions.get(j));
				}
				regions.clear();
				regions.addAll(newRegions);
				//Update long regions in the next sort
				sorted = false;
				size-=remaining;
				remaining = 0;
			}
			changed = true;
		}
		return changed;
	}
	
	@Override
	public boolean retainAll(Collection<?> c) {
		List<T> toRemove = asList();
		toRemove.removeAll(c);
		return removeAll(toRemove);
	}
	/**
	 * Use this method to sort the collection again when the locations of the genomic regions are modified externally
	 */
	public void forceSort() {
		sorted = false;
		sort();
	}
	private void sort() {
		if(!sorted) {
			//System.out.println("Sorting "+regionsForward.size()+" regions");
			//if(regionsForward.size()>0) System.out.println("Type: "+regionsForward.get(0).getClass().getName());
			for(int index:regionsMap.keySet()) {
				List<T> regions = regionsMap.get(index);
				Collections.sort(regions,GenomicRegionPositionComparator.getInstance());
				List<T> longRegions = longRegionsMap.get(index);
				longRegions.clear();
				for(int i=0;i<regions.size();i++) {
					T r = regions.get(i);
		 			//if(i<50) System.out.println("Seqname: "+r.getSequenceName()+". Start: "+r.getStart()+" end:"+r.getEnd());
					int j=i+1;
					for(;j<regions.size() && j<=i+numSpanningLong;j++) {
						GenomicRegion r2 = regions.get(j);
						if(r2.getFirst() > r.getLast()) {
							break;
						}
					}
					if(j>i+numSpanningLong) {
						longRegions.add(r);
					}
				}
				
			}
			//System.out.println("Found "+longRegions.size()+" long regions");
		}
		sorted = true;
	}
	
	public QualifiedSequenceList getSequenceNames() {
		return sequences;
	}
	public GenomicRegionSortedCollection<T> getSequenceRegions(String sequenceName) {
		GenomicRegionSortedCollection<T> answer = new GenomicRegionSortedCollection<T>();
		int sequenceIndex = sequences.indexOf(sequenceName);
		if(sequenceIndex >=0) {
			sort();
			answer.addAll(regionsMap.get(sequenceIndex));
		}
		return answer;
	}
	public GenomicRegionSortedCollection<T> findSpanningRegions(String sequenceName, int position) {
		return findSpanningRegions(sequenceName,position,position);
	}
	public GenomicRegionSortedCollection<T> findSpanningRegions(GenomicRegion region) {
		int index = sequences.indexOf(region.getSequenceName());
		return findSpanningRegions(index, region.getFirst(), region.getLast());
	}
	public GenomicRegionSortedCollection<T> findSpanningRegions(String sequenceName, int first, int last) {
		return findSpanningRegions(sequences.indexOf(sequenceName),first,last);
	}
	public GenomicRegionSortedCollection<T> findSpanningRegions(int sequenceIndex, int first, int last) {
		GenomicRegionSortedCollection<T> answer = new GenomicRegionSortedCollection<T>();
		if(sequenceIndex <0 || sequenceIndex>sequences.size()) return answer;
		sort();
		GenomicRegionSpanComparator spanCmp = GenomicRegionSpanComparator.getInstance();
		GenomicRegion dummyRegion = new GenomicRegionImpl(sequences.get(sequenceIndex).getName(), first, first);
		List<T> regions = regionsMap.get(sequenceIndex);
		int index = Collections.binarySearch(regions, dummyRegion,GenomicRegionPositionComparator.getInstance());
		int indexStart = index;
		if(index<0) {
			index = -index-1;
		}
		indexStart = Math.max(0,index - numSpanningLong);
		
		for(int i=indexStart;i<regions.size();i++) {
			T r = regions.get(i);
			int cmp = spanCmp.compare(r,first, last);
			if (cmp == 0) {
				answer.add(r);
			} else if (cmp > 0) {
				break;
			}
		}
		
		//Look long regions
		List<T> longRegions = longRegionsMap.get(sequenceIndex);
		for(T r: longRegions) {
			int cmp = spanCmp.compare(r,first, last);
			if (cmp == 0 && !answer.contains(r)) {
				answer.add(r);
			} else if (cmp >0) {
				break;
			}
		}
		return answer;
	}
	
	public List<T> asList() {
		sort();
		List<T> answer = new ArrayList<T>();
		for(int i=0;i<sequences.size();i++) {
			answer.addAll(regionsMap.get(i));
		}
		return answer;
	}
	@Override
	public boolean contains(Object o) {
		sort();
		if(o==null || !(o instanceof GenomicRegion)) return false;
		GenomicRegion gr = (GenomicRegion)o;
		int sequenceIndex = sequences.indexOf(gr.getSequenceName());
		if(sequenceIndex<0) return false;
		List<T> regions = regionsMap.get(sequenceIndex);
		int index = GenomicRegionSortedCollection.indexOf (regions,gr);
		return index>=0;
	}
	
	@Override
	public boolean containsAll(Collection<?> c) {
		for(Object o:c) {
			if(!contains(o)) return false;
		}
		return true;
	}
	@Override
	public boolean isEmpty() {
		return size==0;
	}
	@Override
	public Iterator<T> iterator() {
		return asList().iterator();
	}
	
	@Override
	public int size() {
		return size;
	}
	@Override
	public Object[] toArray() {
		return asList().toArray();
	}
	@Override
	public <U> U[] toArray(U[] a) {
		return asList().toArray(a);
	}
}