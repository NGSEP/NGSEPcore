package ngsep.math;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

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
	public Set<T> selectBest(int max) {
		if(max >= countsMap.size()) return new TreeSet<T>(countsMap.keySet());
		List<CountObjectPair<T>> allCounts = new ArrayList<CountObjectPair<T>>(countsMap.values());
		
		Collections.sort(allCounts,new Comparator<CountObjectPair<T>>() {

			@Override
			public int compare(CountObjectPair<T> o1, CountObjectPair<T> o2) {
				return o2.getCount() - o1.getCount();
			}
		});
		Set<T> answer = new TreeSet<T>();
		for(int i=0;i<max;i++) {
			T object = (T)allCounts.get(i).getObject();
			answer.add(object);
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
