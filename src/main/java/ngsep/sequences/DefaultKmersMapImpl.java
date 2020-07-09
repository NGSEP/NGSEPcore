package ngsep.sequences;

import java.io.PrintStream;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Map.Entry;

import ngsep.math.Distribution;

public class DefaultKmersMapImpl implements KmersMap {

	private Map<CharSequence, Short> kmersMap = new HashMap<>();
	@Override
	public int size() {
		return kmersMap.size();
	}

	@Override
	public int getCount(CharSequence kmer) {
		Short count = kmersMap.get(kmer);
		if(count == null) return 0;
		return count;
	}

	@Override
	public void setCount(CharSequence kmer, int count) {
		if(count>Short.MAX_VALUE) count = Short.MAX_VALUE;
		kmersMap.put(kmer, (short) count);
	}

	@Override
	public void addOcurrance(CharSequence kmer) {
		Short count = kmersMap.get(kmer);
		if(count == null) {
			count = 0;
		}
		if(count < Short.MAX_VALUE) count ++;
		kmersMap.put(kmer, count);
	}

	@Override
	public Distribution calculateAbundancesDistribution() {
		Distribution kmerSpectrum = new Distribution(1, 200, 1);
		Iterator<Entry<CharSequence, Short>> it = kmersMap.entrySet().iterator();
		while (it.hasNext()) {
			Entry<CharSequence, Short> entry = it.next();
		    Short value = entry.getValue();
		    kmerSpectrum.processDatapoint(value);
		}
		return kmerSpectrum;
	}

	@Override
	public void filterKmers(int minAbundance) {
		Iterator<Entry<CharSequence, Short>> it = kmersMap.entrySet().iterator();
		while (it.hasNext()) {
			Entry<CharSequence, Short> entry = it.next();
		    if(entry.getValue() < minAbundance) {
		    	it.remove();
		    }
		}
		
	}

	@Override
	public void save(PrintStream out) {
		Iterator<Entry<CharSequence, Short>> it = kmersMap.entrySet().iterator();
		while (it.hasNext()) {
			Entry<CharSequence, Short> entry = it.next();
			out.println(entry.getKey().toString()+"\t"+entry.getValue());
		}
		
	}

}
