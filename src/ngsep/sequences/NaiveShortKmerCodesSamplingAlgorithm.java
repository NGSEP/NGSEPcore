package ngsep.sequences;

public class NaiveShortKmerCodesSamplingAlgorithm implements ShortKmerCodesSamplingAlgorithm {

	@Override
	public boolean[] sample(Integer[] hashcodes) {
		int n = hashcodes.length;
		boolean [] answer = new boolean [n];
		for(int i=0;i<n;i++) {
			answer[i] = hashcodes[i]!=null; 
		}
		return answer;
	}

}
