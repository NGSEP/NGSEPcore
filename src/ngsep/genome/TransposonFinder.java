package ngsep.genome;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import ngsep.alignments.ReadAlignment;
import ngsep.genome.io.SimpleGenomicRegionFileHandler;
import ngsep.math.Distribution;
import ngsep.sequences.DefaultKmersMapImpl;
import ngsep.sequences.QualifiedSequence;

public class TransposonFinder {
	
	/**
	 * Genome of interest
	 */
	private ReferenceGenome genome;
	
	/**
	 * FMIndex of the genome of interest
	 */
	private ReferenceGenomeFMIndex fm;
		
	private Map<String, List<GenomicRegion>> STRs;
		
	private Map<String, List<Transposon>> transposons;
	
	private Distribution distrHits = new Distribution(0, 100, 1);
		
	private int lengthKmer;

	private int minHitSize;

	private boolean useSTRs;
	
	public void run() throws IOException {
		
		// Kmers per subsequence
		int numSequences = genome.getNumSequences();
		for (int i = 0; i < numSequences; i++) {
			QualifiedSequence qs = genome.getSequenceByIndex(i);
			CharSequence seq = qs.getCharacters();
			processSequence(seq, qs.getName(), fm);			
		}
				
		// Save the LTR in a text file
		saveLTRs(transposons, "output.txt");
	}
	
	public void saveLTRs(Map<String, List<Transposon>> LTRs, String filename) throws IOException {
		System.out.printf("Writing LTR predictions on %s \n", filename);
		BufferedWriter writer = new BufferedWriter(new FileWriter(new File(filename)));
		for(String actSeq : LTRs.keySet()) {
			for(Transposon LTR: LTRs.get(actSeq)) {
				System.out.printf("%s\t%d\t%d\t%d\n", actSeq, LTR.getFirst(), LTR.getLast(), LTR.getScore());
				writer.write(String.format("%s\t%d\t%d\t%d", actSeq, LTR.getFirst(), LTR.getLast(), LTR.getScore()));
				writer.newLine();
			}
		}
		writer.close();
	}

	/**
	 * Process a sequence in the genome finding the over-represented genomic regions
	 * @param seq Actual sequence in the genome
	 * @param name Name of the sequence
	 * @param fm FM index of the whole genome
	 */
	public void processSequence(CharSequence seq, String name, ReferenceGenomeFMIndex fm) {
		System.out.printf("Processing Sequence %s \n", name);
		List<Transposon> repetitiveRegions = new ArrayList();
		boolean seen = false;
		int count = 0; // Count of intermediate kmers that are not over-represented
		int maxCount = 30; // TODO: tratar de cambiarlo a distancia
		Transposon actTransposon = null;
		// Take into account known STR
		List<GenomicRegion> actSTRs = null;
		int indexSTR = 0;
		if(useSTRs) {
			actSTRs = STRs.get(name);
		}
		//Subsequence 20bp
		for (int i = 0; i + lengthKmer < seq.length(); i+=10) {
			CharSequence kmer = seq.subSequence(i, (i+lengthKmer));
			List<ReadAlignment> hits = fm.search(kmer.toString());
			distrHits.processDatapoint(hits.size());
			// If the kmer is more than the min hit size
			if(hits.size() > minHitSize ) {
				if(!seen) {
					if(useSTRs && indexSTR < actSTRs.size()) {
						GenomicRegion STR = actSTRs.get(indexSTR);
						if(STR.getFirst() >  i || STR.getLast()<i) { //Check overlap 
							// Is not a tandem repeat
							actTransposon = new Transposon(name, i, (i+lengthKmer), "LTR", hits.size());
							seen = true;
						}
					}
					else {
						actTransposon = new Transposon(name, i, (i+lengthKmer), "LTR", hits.size());
						seen = true;						
					}
				}
				else {
					int actHits = actTransposon.getScore();
					actTransposon.setLast((i+lengthKmer));
					actTransposon.setScore(actHits + hits.size());
					count = 0;
				}
			}
			else if(seen) {
				count ++;
			}
			if(count > maxCount) {
				seen = false;
				count = 0;
				if(actTransposon.length() - (lengthKmer + 1) > 1000) { 
					repetitiveRegions.add(actTransposon);
				}
			}
			// Check that the index of the STR is greater than i
			if(useSTRs && indexSTR < actSTRs.size()) {
				GenomicRegion STR = actSTRs.get(indexSTR);
				while(STR.getFirst() < i && (++indexSTR) < actSTRs.size() - 1) {
					STR = actSTRs.get(indexSTR);
				}
			}
		}
		System.out.printf("Found %d repetitive regions \n",repetitiveRegions.size());
		transposons.put(name, repetitiveRegions);
	}

	public static void main(String[] args) throws IOException {
		long start = System.currentTimeMillis();
		TransposonFinder instance = new TransposonFinder();
		// Load the genome from a .fa file
		instance.genome = new ReferenceGenome(args[0]);
		// Load short tandem repeats from the .fa file
		SimpleGenomicRegionFileHandler loader = new SimpleGenomicRegionFileHandler();
		instance.STRs = loader.loadRegionsAsMap(args[1]);
		// Put as "arguments" kmer length and min hit size
		instance.lengthKmer = 20;
		instance.minHitSize = 10;
		instance.transposons = new LinkedHashMap();
		instance.useSTRs = false;
		// FM Index
		instance.fm = new ReferenceGenomeFMIndex(instance.genome);
		// Find transposable elements
		instance.run();
		System.out.println("Total time consumed: " + (System.currentTimeMillis() - start));
	}
}
