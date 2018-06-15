package ngsep.assembly;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Map;
import java.util.Map.Entry;
import java.util.OptionalInt;
import java.util.stream.Stream;

import htsjdk.samtools.SamReaderFactory.Option;

public class MaxPath extends GraphSimplificator {

    public MaxPath(int N, Map<Integer, Map<Integer, Integer>> Edges) {
	super(N, Edges);
    }

    private ArrayList<Integer> path=new ArrayList<>();
    private int[] fn = new int[N];
    private int[] sig = new int[N];

    public ArrayList<Integer> fn() {
	int[] top = getTopologicalOrder();
	
	for (int i = top.length - 1; i >= 0; i--) {
	    int n = top[i]; 
	    int max = 0;
	    int posMax = -1;
	    for (Entry<Integer, Integer> x : Edges.get(n).entrySet()) {
		int peso = fn[x.getKey()] + x.getValue();
		if (max < peso) {
		    posMax = x.getKey();
		    max = peso;
		}
	    }
	    fn[n] = max;
	    sig[n] = posMax;
	}
	int  pos = getMaxPosition();
	while(pos!=-1) {
	    path.add(pos);
	    pos=sig[pos];  
	}
	System.out.println(path.size());
	return path;
    }
    
    private int getMaxPosition() {
	int temp = 0;
	int pos = 0;
	for(int i = 0; i<fn.length;i++) {
	    if(fn[i]>temp) {
		pos=i;
		temp=fn[i];
	    }
	}
	return pos;
    }
    
    public ArrayList<Integer> getPath(){
	return path;
    }
}
