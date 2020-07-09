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
package ngsep.math;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

/**
 * Helper class to collect counts for a collection of objects and select the best ranked
 * @author Jorge Duitama
 * @param <T> Objects to be counted
 */
public class CountsRankHelper <T> {
	private Map<T,CountObjectPair<T>> countsMap = new TreeMap<T,CountObjectPair<T>>();
	public void add (T o) {
		CountObjectPair<T> count = countsMap.get(o);
		if(count == null) {
			countsMap.put(o, new CountObjectPair<T>(o));
		} else {
			count.addCount();
		}
	}
	public LinkedHashMap<T,Integer> selectBest(int max) {
		
		List<CountObjectPair<T>> allCounts = new ArrayList<CountObjectPair<T>>(countsMap.values());
		
		Collections.sort(allCounts,new Comparator<CountObjectPair<T>>() {

			@Override
			public int compare(CountObjectPair<T> o1, CountObjectPair<T> o2) {
				return o2.getCount() - o1.getCount();
			}
		});
		LinkedHashMap<T,Integer> answer = new LinkedHashMap<>();
		int n = allCounts.size();
		for(int i=0;i<max && i<n;i++) {
			CountObjectPair<T> pair = allCounts.get(i); 
			T object = pair.getObject();
			answer.put(object,pair.getCount());
		}
		return answer;
	}
	public int getNumDifferent() {
		return countsMap.size();
	}
}
class CountObjectPair <T> {
	private T o;
	private int count;
	public CountObjectPair (T o) {
		this.o = o;
		this.count = 1;
	}
	public void addCount() {
		count++;
	}
	public T getObject() {
		return o;
	}
	public int getCount() {
		return count;
	}
}
