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
		
	private Map<String, List<GenomicRegion>> crossRegions;
	
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
		
		// TODO: Check the form of LRT in the candidates
		
		// Save the LTR in a text file
		saveLTRs(crossRegions, "output.txt");
	}
	
	public void saveLTRs(Map<String, List<GenomicRegion>> LTRs, String filename) throws IOException {
		System.out.printf("Writing LTR predictions on %s \n", filename);
		BufferedWriter writer = new BufferedWriter(new FileWriter(new File(filename)));
		writer.write("LTR Prediction");
		writer.newLine();
		writer.write(String.format("%s\t%s\t%s", "Id", "Initial Position", "Final Position"));
		writer.newLine();
		for(String actSeq : LTRs.keySet()) {
			for(GenomicRegion LTR: LTRs.get(actSeq)) {
				System.out.printf("%s\t%d\t%d\n", actSeq, LTR.getFirst(), LTR.getLast());
				writer.write(String.format("%s\t%d\t%d", actSeq, LTR.getFirst(), LTR.getLast()));
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
		List<GenomicRegion> repetitiveRegions = new ArrayList();
		boolean seen = false;
		int count = 0; // Count of intermediate kmers that are not over-represented
		int maxCount = 10; // TODO: tratar de cambiarlo a distancia
		GenomicRegionImpl actGenomicRegion = null;
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
						if(Math.abs(STR.getFirst() - i) > 200) {
							// Is not a tandem repeat
							actGenomicRegion = new GenomicRegionImpl(name, i, (i+lengthKmer));
							seen = true;
						}
					}
					else {
						actGenomicRegion = new GenomicRegionImpl(name, i, (i+lengthKmer));
						seen = true;						
					}
				}
				else if (seen && count <= maxCount) {
					actGenomicRegion.setLast((i+lengthKmer));					
				}
			}
			else if(seen) {
				count ++;
			}
			if(count > maxCount) {
				seen = false;
				count = 0;
				if(actGenomicRegion.length() > (lengthKmer + 1)) { 
					repetitiveRegions.add(actGenomicRegion);
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
		crossRegions.put(name, repetitiveRegions);
	}

	public static void main(String[] args) throws IOException {
		TransposonFinder instance = new TransposonFinder();
		// Load the genome from a .fa file
		instance.genome = new ReferenceGenome(args[0]);
		// Load short tandem repeats from the .fa file
		SimpleGenomicRegionFileHandler loader = new SimpleGenomicRegionFileHandler();
		instance.STRs = loader.loadRegionsAsMap(args[1]);
		// Put as "arguments" kmer length and min hit size
		instance.lengthKmer = 20;
		instance.minHitSize = 10;
		instance.crossRegions = new LinkedHashMap();
		instance.useSTRs = true;
		// FM Index
		instance.fm = new ReferenceGenomeFMIndex(instance.genome);
		// Find transposable elements
		instance.run();
	}
}
