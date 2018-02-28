package ngsep.math;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

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
		for(int i=0;i<max;i++) {
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
