package ngsep.genome;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import ngsep.genome.io.SimpleGenomicRegionFileHandler;
import ngsep.math.Distribution;
import ngsep.sequences.UngappedSearchHit;
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

	private Map<String, List<GenomicRegion>> goldStandard;

	private Map<String, List<Transposon>> transposons;

	private Map<String, List<Annotation>> hitNumbers;

	private Distribution distrHits = new Distribution(0, 100, 1);

	private int lengthKmer;

	private int minHitSize;

	private boolean useSTRs;

	private boolean useGS;

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
		System.out.printf("Writing gt statistics on %s \n", "gt_" + filename);
		BufferedWriter w = new BufferedWriter(new FileWriter(new File("gt_" + filename)));
		for(String actSeq : LTRs.keySet()) {
			for(Annotation ann : hitNumbers.get(actSeq)) {
				System.out.printf("%s\t%d\t%d\t%d\t%f\n", actSeq, ann.first, ann.last, ann.maxHit, ann.medianHit);
				w.write(String.format("%s\t%d\t%d\t%d\t%f", actSeq, ann.first, ann.last, ann.maxHit, ann.medianHit));
				w.newLine();
				if(ann.maxHit <= 10){
					CharSequence act_char_seq = this.genome.getReference(actSeq, ann.first, ann.last);
					System.out.printf("Act seq: %s\n", act_char_seq.toString());
					w.write(String.format("%s", act_char_seq.toString()));
					w.newLine();
				}
			}
		}
		writer.close();
		w.close();
	}

	private class Annotation {
		public int maxHit;
		public String name;
		public int first;
		public int last;
		public double medianHit;

		public Annotation(String name, int first, int last, int maxHit, int sumHits, int numKmers){
			this.name = name;
			this.first = first;
			this.last = last;
			this.maxHit = maxHit;
			if(numKmers == 0){
				this.medianHit = 0;
			}
			else{
				this.medianHit = sumHits / (double) numKmers;
			}
		}

		public String toString(){
			return String.format("%s\t%d\t%d\t%d", this.name, this.first, this.last, this.maxHit);
		}
	}

	/**
	 * Process a sequence in the genome finding the over-represented genomic regions
	 * @param seq Actual sequence in the genome
	 * @param name Name of the sequence
	 * @param fm FM index of the whole genome
	 */
	public void processSequence(CharSequence seq, String name, ReferenceGenomeFMIndex fm) {
		System.out.printf("Processing Sequence %s \n", name);
		List<Transposon> repetitiveRegions = new ArrayList<>();
		boolean seen = false;
		int count = 0; // Count of intermediate kmers that are not over-represented
		int maxCount = 30; // TODO: tratar de cambiarlo a distancia
		Transposon actTransposon = null;
		// Take into account known STR
		List<GenomicRegion> actSTRs = null;
		int indexSTR = 0;
		List<GenomicRegion> actGS = null;
		int indexGS = 0;
		int maxHitGS = 0;
		int countHits = 0;
		int numKmers = 0;
		List<Annotation> hitsGS = new ArrayList<>();
		if(useSTRs) {
			actSTRs = STRs.get(name);
		}
		if(useGS){
			actGS = goldStandard.get(name);
		}
		//Subsequence 20bp
		for (int i = 0; i + lengthKmer < seq.length(); i+=10) {
			String kmer = seq.subSequence(i, (i+lengthKmer)).toString().toUpperCase();
			List<UngappedSearchHit> hits = fm.exactSearch(kmer);
			distrHits.processDatapoint(hits.size());
			// Save the number of hits for the GS
			if(useGS && actGS != null && indexGS < actGS.size()) {
				GenomicRegion gs = actGS.get(indexGS);
				if(gs.getLast() < i){
					Annotation newGS = new Annotation(name, gs.getFirst(), gs.getLast(), maxHitGS, countHits, numKmers);
					hitsGS.add(newGS);
					indexGS ++;
					maxHitGS = 0;
					countHits = 0;
					numKmers = 0;
				}
				else if(gs.getFirst() >= i){
					//Inside a groundtruth transposon
					numKmers ++;
					countHits += hits.size();
					if(hits.size() > maxHitGS){
						maxHitGS = hits.size();
					}
				}
			}
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
		hitNumbers.put(name, hitsGS);
	}

	public static void main(String[] args) throws IOException {
		long start = System.currentTimeMillis();
		TransposonFinder instance = new TransposonFinder();
		// Load the genome from a .fa file
		instance.genome = new ReferenceGenome(args[0]);
		// Load short tandem repeats from the .fa file
		instance.useSTRs = false;
		if(instance.useSTRs){
			SimpleGenomicRegionFileHandler loader = new SimpleGenomicRegionFileHandler();
			instance.STRs = loader.loadRegionsAsMap(args[1]);
		}
		else{
			instance.STRs = null;
		}
		// Put as "arguments" kmer length and min hit size
		instance.lengthKmer = 20;
		instance.minHitSize = 10;
		instance.transposons = new LinkedHashMap<>();
		instance.hitNumbers = new LinkedHashMap<>();
		instance.useGS = false;
		if(instance.useGS){
			SimpleGenomicRegionFileHandler gs_loader = new SimpleGenomicRegionFileHandler();
			instance.goldStandard = gs_loader.loadRegionsAsMap(args[2]);
		}
		else{
			instance.goldStandard = null;
		}
		// FM Index
		instance.fm = new ReferenceGenomeFMIndex(instance.genome);
		// Find transposable elements
		instance.run();
		System.out.println("Total time consumed: " + (System.currentTimeMillis() - start));
	}
}
