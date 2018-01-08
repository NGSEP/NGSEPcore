package ngsep.sequences;

import java.io.IOException;
import java.io.PrintStream;
import java.util.Iterator;
import java.util.Map;

import ngsep.main.CommandsDescriptor;
import ngsep.sequences.io.FastqFileReader;

public class FastqFileErrorCorrector {
	
	private Map<CharSequence, Short> kmersMap;
	private int kmerSize = KmersCounter.DEFAULT_KMER_SIZE;
	private int minCount = 5;
	private int correctedErrors = 0;
	/**
	 * @return the minCount
	 */
	public int getMinCount() {
		return minCount;
	}

	/**
	 * @param minCount the minCount to set
	 */
	public void setMinCount(int minCount) {
		this.minCount = minCount;
	}
	
	/**
	 * @param minCount the minCount to set
	 */
	public void setMinCount(Integer minCount) {
		this.setMinCount(minCount.intValue());
	}

	public static void main(String[] args) throws Exception {
		FastqFileErrorCorrector instance = new FastqFileErrorCorrector();
		//int i = CommandsDescriptor.getInstance().loadOptions(instance, args);
		int i=0;
		String inFilename = args[i++];
		String outFilename = args[i++];
		instance.process(inFilename,outFilename);

	}

	public void process(String inFilename, String outFilename) throws IOException {
		correctedErrors = 0;
		System.out.println("Calculating k-mers map from: "+inFilename);
		KmersCounter counter = new KmersCounter();
		counter.setKmerSize(kmerSize);
		counter.processFile(inFilename);
		kmersMap = counter.getKmersMap();
		kmerSize = counter.getKmerSize();
		System.out.println("Extracted "+kmersMap.size()+" k-mers from: "+inFilename);
		System.out.println("Processing file: "+inFilename);
		try (FastqFileReader reader = new FastqFileReader(inFilename);
			PrintStream out = new PrintStream(outFilename)) {
			Iterator<RawRead> it = reader.iterator();
			while (it.hasNext()) {
				RawRead read = it.next();
				processRead (read);
				read.save(out);
			}
		}
		System.out.println("Corrected "+correctedErrors+" potential errors. Output written to "+outFilename);
	}

	public void processRead(RawRead read) {
		for(int h=0;h<3;h++) {
			String readStr = read.getCharacters().toString();
			char [] readChars = readStr.toCharArray();
			CharSequence [] readKmers = KmersCounter.extractKmers(readStr, kmerSize , true);
			int [] readKmerCounts = new int [readKmers.length];
			
			for(int i=0;i<readKmers.length;i++) {
				CharSequence kmer = readKmers[i];
				if (kmer==null) readKmerCounts[i] = 0;
				else {
					Short count = kmersMap.get(kmer);
					readKmerCounts[i] = (count!=null)?count:0;
				}
			}
			int lastRepresented= -1;
			boolean corrected = false;
			for(int i=0;i<readKmers.length;i++) {
				if(readKmerCounts[i] >= minCount) {
					if(i-1!=lastRepresented) {
						corrected = corrected || correctErrors (readChars,lastRepresented,i);
					}
					lastRepresented = i;
				}
			}
			if (corrected) {
				read.setCharacters(new String(readChars));
			} else break;
		}
		
	}

	private boolean correctErrors(char [] readChars, int lastRepresented, int nextRepresented) {
		int first = 0;
		if(lastRepresented>=0) first = lastRepresented+kmerSize;
		int last = nextRepresented-1;
		if(last-first < 0) return false;
		double bestScore = getScore(readChars,lastRepresented+1,last);
		int bestI=-1;
		char bestBP = 0;
		for(int i=first;i<=last;i++) {
			char origBP = readChars[i];
			for(int j=0;j<4;j++) {
				char changeBP = DNASequence.BASES_STRING.charAt(j);
				if(changeBP != origBP) {
					readChars[i] = changeBP;
					double score = getScore(readChars,lastRepresented+1,last);
					if(score > bestScore) {
						bestScore = score;
						bestI = i;
						bestBP = changeBP;
					}
				}
			}
			readChars[i] = origBP;
		}
		if(bestI>=0) {
			readChars[bestI] = bestBP;
			correctedErrors++;
			return true;
		}
		return false;
	}

	private double getScore(char[] readChars, int first, int last) {
		CharSequence [] readKmers = KmersCounter.extractKmers(new String (readChars), kmerSize , first, last, true);
		double score = 0;
		for(int i=first;i<=last&& i<readKmers.length;i++) {
			CharSequence kmer = readKmers[i];
			if (kmer!=null) {
				Short count = kmersMap.get(kmer);
				if(count!=null) score+=count;
			}
		}
		return score;
	}

}
