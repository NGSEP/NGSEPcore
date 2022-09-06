package ngsep.discovery;

import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import JSci.maths.statistics.GammaDistribution;
import JSci.maths.statistics.NormalDistribution;
import JSci.maths.statistics.SampleStatistics;
import ngsep.alignments.ReadAlignment;
import ngsep.alignments.io.ReadAlignmentFileReader;
import ngsep.genome.GenomicRegion;
import ngsep.genome.GenomicRegionSortedCollection;
import ngsep.genome.GenomicRegionSpanComparator;
import ngsep.genome.ReferenceGenome;
import ngsep.sequences.QualifiedSequenceList;
import ngsep.variants.CalledGenomicVariant;
import ngsep.variants.CalledGenomicVariantImpl;
import ngsep.variants.GenomicVariant;
import ngsep.variants.GenomicVariantAnnotation;
import ngsep.variants.GenomicVariantImpl;
import ngsep.variants.Sample;
import ngsep.vcf.VCFFileHeader;
import ngsep.vcf.VCFFileWriter;
import ngsep.vcf.VCFRecord;

public class LongReadStructuralVariantDetector implements LongReadVariantDetectorAlgorithm {
	
	public static final int DEF_MIN_MQ_UNIQUE_ALIGNMENT = ReadAlignment.DEF_MIN_MQ_UNIQUE_ALIGNMENT;
	public static final String FILTER_SETTING_MISSING = ".";
	public static final String KEY_SEPARATOR = ",";
	public static final String DEFAULT_VARIANT_ID_PREFIX = "NGSEP.";
	
	public static final String[] GENOTYPE_ALLELES = {"REF", "ALT"};
	public static final double DEF_HET_RATE = 0.5;
	public static final double HETEROCIGOZITY_TRESHOLD = 0.8;
	public static final double DEF_SHAPE_PARAMETER_ERROR_DIST = 0.5;
	public static final double DEF_LOG_REF_PROB_SV = Math.log10(0.9999);
	public static final double DEF_LOG_ALT_PROB_SV = Math.log10(0.0001);
	public static final double DEF_LOG_ERROR_PROB_SV = Math.log10(0.00001);

	
	public static final int INV_DETERMINING_MAX_DISTANCE = 800;
	public static final int DEL_INTER_DETERMINING_MAX_DISTANCE = 90000;
	public static final String MAX_CLIQUE_FINDER_ALGORITHM = "MCC";
	public static final String DBSCAN_ALGORITHM = "DBSCAN";
	public static final String STRONGLY_CONNECTED_COMPONENTS_ALGORITHM = "SCC";
	
	private GenomicRegionSortedCollection<GenomicVariant> signatures;
	private GenomicRegionSortedCollection<GenomicVariant> variants;
	private Map<String, List<SimplifiedReadAlignment>> alignments;
	private ReferenceGenome refGenome;
	private String clusteringAlgorithm = DBSCAN_ALGORITHM;
	private int lengthToDefineSVEvent = 50;
	private int minMQ = DEF_MIN_MQ_UNIQUE_ALIGNMENT;
	
	LongReadStructuralVariantDetector(){
	}
	
	
	public void setRefGenome(ReferenceGenome refGenome) {
		this.refGenome = refGenome;
	}


	public void setRefGenome(String referenceFile) throws IOException {
		this.refGenome = new ReferenceGenome(referenceFile);
	}
	
	public void setClusteringAlgorithm(String algorithm) throws Exception {
		if(DBSCAN_ALGORITHM.equals(algorithm) || MAX_CLIQUE_FINDER_ALGORITHM.equals(algorithm) ||
				STRONGLY_CONNECTED_COMPONENTS_ALGORITHM.equals(algorithm)) {
			this.clusteringAlgorithm = algorithm;
		}
		else {
			throw new Exception("Algorithm must be either DBSCAN, SCC or MCC");
		}
	}
	
	private void readAlignments(String alignmentFile,  int SVLength) throws IOException{
		try(ReadAlignmentFileReader alignmentReader = new ReadAlignmentFileReader(alignmentFile, refGenome)){
			setIndelTresholdSize(SVLength);
			int filterFlags = ReadAlignment.FLAG_READ_UNMAPPED;
			alignmentReader.setFilterFlags(filterFlags);
			alignmentReader.setMinMQ(minMQ);
			QualifiedSequenceList qualifiedSequenceNames = refGenome.getSequencesList();
			signatures = new GenomicRegionSortedCollection<>(qualifiedSequenceNames);
			alignments = new HashMap<>();
			Iterator<ReadAlignment> it = alignmentReader.iterator();
			ReadAlignment aln = it.next();
			while(it.hasNext()) {
				aln = it.next();
				computeAlignment(aln);
			}
			for(List <SimplifiedReadAlignment> alnRegions : alignments.values()) {
				if(alnRegions.size() > 1) {
					findInterAlnSignatures(alnRegions);
				}
			}
			signatures.forceSort();
		}
	}
	
	private void computeAlignment(ReadAlignment aln) {
		// TODO Auto-generated method stub
		List<SimplifiedReadAlignment> alnRegions = alignments.computeIfAbsent(aln.getReadName(), v -> new ArrayList<>());
		SimplifiedReadAlignment simpleAln = new SimplifiedReadAlignment(aln);
		if (!aln.isSecondary()) simpleAln.setOriginFromSecondaryAlignment();
		int alnIdx = alnRegions.size();
		Map<Integer,GenomicVariant> filteredCalls = findIntraAlnSignatures(aln, alnIdx);
		alnRegions.add(simpleAln);
		simpleAln.setIndelCalls(filteredCalls);
	}
	
	public GenomicRegionSortedCollection<GenomicVariant> getSignatures() {
		return signatures;
	}
	
	public ReferenceGenome getRefGenome() {
		return refGenome;
	}
	
	private GenomicVariant createIndelIntraAlnSignature(GenomicVariant indel) {
		int first = indel.getFirst();
		int last = indel.getLast();
		int length = indel.length();
		String seqName = indel.getSequenceName(); 
		byte type = first + 1 == last ? GenomicVariant.TYPE_LARGEINS : GenomicVariant.TYPE_LARGEDEL;
		Signature signature = new Signature(seqName, first, last, type);
		signature.setLength(length);
		signature.setSignatureType(Signature.INTRAALIGNMENT);
		return signature;
	}
	
	private Map<Integer, GenomicVariant> findIntraAlnSignatures(ReadAlignment aln, int alnIdx) {
		Map<Integer, GenomicVariant> calls = aln.getIndelCalls();
		Map<Integer, GenomicVariant> filteredCalls = new HashMap<>();
		if(calls != null){
			for (Map.Entry<Integer, GenomicVariant> call : calls.entrySet()) {
				GenomicVariant indel = call.getValue();
				if(indel.length() < lengthToDefineSVEvent) continue;
				Signature sign = (Signature) createIndelIntraAlnSignature(indel);
				if(aln.isSecondary()) sign.setOriginFromSecondaryAlignment();
				filteredCalls.put(call.getKey(), sign);
				sign.setReadName(aln.getReadName());
				sign.setAlnIdx(alnIdx);
				signatures.add(sign);
			}
		}
		return filteredCalls; 
	}
	
	private void findInterAlnSignatures(List<SimplifiedReadAlignment> alnRegions) {
		// TODO Auto-generated method stub
		int nAln = alnRegions.size();
		SimplifiedReadAlignment firstAln = alnRegions.get(0);
		SimplifiedReadAlignment lastAln = alnRegions.get(nAln - 1);
		if(nAln == 2) {
			computeInterAlnIndel(firstAln, lastAln, 0);
				}
		else if(nAln >= 3) {
			for(int i = 0; i < nAln - 1; i++) {
				int j = i + 1;
				SimplifiedReadAlignment currentAln = alnRegions.get(i);
				SimplifiedReadAlignment secondAln = alnRegions.get(j);
				computeInterAlnIndel(currentAln, secondAln, i);				
			}
			for(int i = 0; i < nAln - 2; i++) {
				int j = i + 1;
				int k = i + 2;
				SimplifiedReadAlignment currentAln = alnRegions.get(i);
				SimplifiedReadAlignment secondAln = alnRegions.get(j);
				SimplifiedReadAlignment thirdAln = alnRegions.get(k);
				computeInversions(currentAln, secondAln, thirdAln, j);
			}
		}
	}

	private boolean fulfillsInversionAlnOrientations(SimplifiedReadAlignment firstAln,
			SimplifiedReadAlignment secondAln, SimplifiedReadAlignment thirdAln) {
		// TODO Auto-generated method stub
		if(firstAln.isPositiveStrand() == thirdAln.isPositiveStrand() 
				&& firstAln.isPositiveStrand() != secondAln.isPositiveStrand()) return true;
		return false;
	}

	private int estimateInterAlnLength(SimplifiedReadAlignment firstAln, SimplifiedReadAlignment lastAln) {
		int lengthToSubstract;
		int softClipValue;
		if(firstAln.getSoftClipEnd() > lastAln.getSoftClipStart()) {
			softClipValue = firstAln.getSoftClipEnd();
			lengthToSubstract = lastAln.length();
		}
		else {
			softClipValue = lastAln.getSoftClipStart();
			lengthToSubstract = firstAln.length();
		}
		return softClipValue - lengthToSubstract;
	}
	
	private void computeInterAlnIndel(SimplifiedReadAlignment firstAln, SimplifiedReadAlignment lastAln, int firstAlnIdx) {
		int distance = Math.abs(lastAln.getFirst() - firstAln.getLast());
		String sequenceName = firstAln.getSequenceName();
		int first = firstAln.getLast() + 1;
		int last = lastAln.getFirst();
		int length = last - first + 1;
		int interAlnDistance = estimateInterAlnLength(firstAln, lastAln);
		if(distance >= lengthToDefineSVEvent && interAlnDistance <= 100 && distance < DEL_INTER_DETERMINING_MAX_DISTANCE) { //distance <= DEL_INTER_DETERMINING_MAX_DISTANCE && getInterAlnDistance(firstAln, lastAln) < lengthToDefineSVEvent) {
			byte type = GenomicVariant.TYPE_LARGEDEL;
			if(length >= lengthToDefineSVEvent) {
				Signature signature = (Signature) createIndelInterAlnSignature(sequenceName, first, last, length, type);
				signature.setReadName(firstAln.getReadName());
				signature.setAlnIdx(firstAlnIdx);
				signatures.add(signature);
				Map<Integer, GenomicVariant> indelCalls = firstAln.getIndelCalls();
				indelCalls.put(first, signature);
			}
		}
		else if(interAlnDistance >= lengthToDefineSVEvent && distance <= lengthToDefineSVEvent && 
				firstAln.isPositiveStrand() == lastAln.isPositiveStrand()){
			byte type = GenomicVariant.TYPE_LARGEINS;
			last = first + 1;
			length = interAlnDistance;
			Signature signature = (Signature) createIndelInterAlnSignature(sequenceName, first, last, length, type);
			signatures.add(signature);
			GenomicRegionSpanComparator spanCmp = GenomicRegionSpanComparator.getInstance();
			Map<Integer, GenomicVariant> indelCalls;
			if(spanCmp.compare(firstAln, signature) == 0) {
				indelCalls = firstAln.getIndelCalls();
				indelCalls.put(first, signature);
				signature.setReadName(firstAln.getReadName());
				signature.setAlnIdx(firstAlnIdx);
			}
			else {
				indelCalls = lastAln.getIndelCalls();
				indelCalls.put(first, signature);
				signature.setReadName(lastAln.getReadName());
				signature.setAlnIdx(firstAlnIdx + 1);
			}
		}
	}
	
	private void computeInversions(SimplifiedReadAlignment firstAln, SimplifiedReadAlignment secondAln,
			SimplifiedReadAlignment thirdAln, int middleAlnIdx) {
		// TODO Auto-generated method stub
		int firstDistance = Math.abs(secondAln.getFirst() - firstAln.getLast());
		int secondDistance = Math.abs(thirdAln.getFirst() - secondAln.getLast());
		if(secondDistance <= INV_DETERMINING_MAX_DISTANCE && firstDistance <= INV_DETERMINING_MAX_DISTANCE && 
				fulfillsInversionAlnOrientations(firstAln, secondAln, thirdAln)) {
			String sequenceName = secondAln.getSequenceName();
			int first = secondAln.getFirst();
			int last = secondAln.getLast();
			int length = last - first + 1;
			byte type = GenomicVariant.TYPE_INVERSION;
			if(length >= lengthToDefineSVEvent) {
				Signature signature = (Signature) createIndelInterAlnSignature(sequenceName, first, last, length, type);
				signature.setReadName(secondAln.getReadName());
				signature.setAlnIdx(middleAlnIdx);
				signatures.add(signature);
				Map<Integer, GenomicVariant> indelCalls = secondAln.getIndelCalls();
				indelCalls.put(first, signature);
			}
		}
	}
	
	private GenomicVariant createIndelInterAlnSignature(String seqName, int first, int last, int length,
			byte type) {
		// TODO Auto-generated method stub
		Signature signature = new Signature(seqName, first, last, type);
		signature.setLength(length);
		signature.setSignatureType(Signature.INTERALIGNMENT);
		return signature;
	}
	
	@Override
	public GenomicRegionSortedCollection<GenomicVariant> callVariants(Map<String, List<List<Integer>>> clusters){
		List<GenomicVariant> variants = new ArrayList<>();
		List<String> keys = new ArrayList<>(clusters.keySet());
		QualifiedSequenceList sequences = new QualifiedSequenceList(refGenome.getSequencesList());
		GenomicRegionSortedCollection<GenomicVariant> sortedVariants = new GenomicRegionSortedCollection<>(sequences);
		double secondaryAlnExcludingTreshold = 0.5;
		int[] nId = new int[20];
		for(String k:keys) {
			List<List<Integer>> chrClusters = clusters.get(k);
			List<GenomicVariant> chrSignList = signatures.getSequenceRegions(k).asList();
			for(List<Integer> cluster:chrClusters) {
				List<Signature> clusterSigns = new ArrayList<>();
				boolean comesFromSecondaryAln = false;
				int secondaryAlnSignCount = 0;
				for(int idx:cluster) {
					Signature clusteredSign = (Signature) chrSignList.get(idx);
					if(clusteredSign.comesFromSecondaryAlignment()) secondaryAlnSignCount++; 
					clusterSigns.add(clusteredSign);
				}
				double secondaryAlnSignPercentage = (double) secondaryAlnSignCount/cluster.size();
				if(secondaryAlnSignPercentage >= secondaryAlnExcludingTreshold) comesFromSecondaryAln = true; 
				if(comesFromSecondaryAln) continue;
				Collections.sort(clusterSigns, Comparator.comparingInt(s -> s.getFirst()));
				processClusterToVariant(clusterSigns, k, variants, nId);
			}
		}
		sortedVariants.addAll(variants);
		sortedVariants.forceSort();
		return sortedVariants;
	}
	
	private void processClusterToVariant(List<Signature> clusterSigns, String sequenceName, 
			List<GenomicVariant> variants, int[] varNumber){
		GenomicVariant firstVar = clusterSigns.get(0);
		SampleStatistics calcFirst = new SampleStatistics();
		SampleStatistics calcLast = new SampleStatistics();
		int endOfSpan = -1;
		for(GenomicVariant v:clusterSigns) {
			int currentEndOfSpan = v.getFirst() + v.length() + 1;
			calcFirst.update(v.getFirst());
			calcLast.update(currentEndOfSpan);
		}
		int first = (int) calcFirst.getMean();
		int last = (int) calcLast.getMean();
		endOfSpan = last;
		byte type = firstVar.getType();
		List<String> alleles = new ArrayList<>();
		String ref = "";
		String alt = "";
		if(GenomicVariant.TYPE_LARGEINS == type) {
			last = first+1;
			ref = "" + refGenome.getReferenceBase(sequenceName, first);
			alt = "<" + GenomicVariant.TYPENAME_LARGEINS + ">";
		}
		else if(GenomicVariant.TYPE_LARGEDEL == type) {
			ref = "" + refGenome.getReferenceBase(sequenceName, first);
			alt = "<" + GenomicVariant.TYPENAME_LARGEDEL + ">";
		}
		else if(GenomicVariant.TYPE_INVERSION == type) {
			ref = "" + refGenome.getReferenceBase(sequenceName, first);
			alt = "<" + GenomicVariant.TYPENAME_INVERSION + ">";
		}
		alleles.add(ref);
		alleles.add(alt);
		GenomicVariantImpl variant = new GenomicVariantImpl(sequenceName, first, alleles);
		variant.setLength(endOfSpan - first + 1);
		variant.setLast(last);
		variant.setType(type);
		//short variantScore = calculateVariantScore(clusterSigns, variant.length());
		//variant.setVariantQS(variantScore);
		String id = DEFAULT_VARIANT_ID_PREFIX + GenomicVariantImpl.getVariantTypeName(type) + "." + varNumber[type];
		variant.setId(id);
		if(variant.length() < lengthToDefineSVEvent) return;
		variants.add(variant);
		int nSgins = clusterSigns.size();
		String variantId = variant.getId();
		for(int i = 0; i < nSgins; i++) {
			Signature signature = (Signature) clusterSigns.get(i);
			String readName = signature.getReadName();
			int alnIdx = signature.getAlnIdx();
			SimplifiedReadAlignment aln = alignments.get(readName).get(alnIdx);
			aln.addCallByVariantId(variantId, signature);
		}
		boolean testForDuplications = GenomicVariant.TYPE_LARGEINS == variant.getType();
		if(testForDuplications) determineInsertionAsDuplication(variant, clusterSigns);
		varNumber[variant.getType()]++;
	}
	
	private void determineInsertionAsDuplication(GenomicVariant variant, List<Signature> clusterSigns) {
		int nSigns = clusterSigns.size();
		boolean isDuplication;
		int tandemRepeatSupportingPairings = 0;
		int interspersedTandemRepeatSupportingPairings = 0;
		int possibleDupPairings = (nSigns*(nSigns+1))/2 - nSigns;
		for(int i = 0; i < nSigns; i++) {
			Signature sign = clusterSigns.get(i);
			for(int j = i; j < nSigns; j++) {
				if(i==j) continue;
				Signature nextSign = clusterSigns.get(j);
				int distance = nextSign.getFirst() - sign.getLast();
				if(distance < 0) continue;
				if(distance <= lengthToDefineSVEvent*2) tandemRepeatSupportingPairings++;
				else interspersedTandemRepeatSupportingPairings++;
			}
		}
		isDuplication = (tandemRepeatSupportingPairings + interspersedTandemRepeatSupportingPairings) > 
			possibleDupPairings*0.8;
		if(isDuplication) {
			variant.setType(GenomicVariant.TYPE_DUPLICATION);
			List<String> alleles = Arrays.asList(variant.getAlleles());
			alleles.set(1, "<" + GenomicVariant.TYPENAME_DUPLICATION + ">");
			GenomicVariantImpl impl = (GenomicVariantImpl) variant;
			impl.setAlleles(alleles);
		}
	}
	/**
	private short calculateVariantScore(List<GenomicVariant> candidates, int length) {
		double score = 0;
		SampleStatistics calcPos = new SampleStatistics();
		SampleStatistics calcSpan = new SampleStatistics();
		for(GenomicRegion sign:candidates) {
			calcPos.update(sign.getFirst());
			calcSpan.update(sign.length());
		}
		int n = candidates.size();
		double numSign = Math.min(80, n);
		double stdSpan = Math.sqrt(calcSpan.getVariance());
		double stdPos = Math.sqrt(calcPos.getVariance());
		double spanDeviationScore = 1 - Math.min(1, (double) stdSpan/length);
		double posDeviationScore = 1 - Math.min(1, (double) stdPos/length);
		score = numSign + (spanDeviationScore*((double) numSign/8)) + (posDeviationScore*((double) numSign/8));
		return (short) score;
	}**/

	private GenomicRegionSortedCollection<CalledGenomicVariant> makeBayesianGenotypeCalls(List<GenomicVariant> variantsList) {
		GenomicRegionSortedCollection<CalledGenomicVariant> genotypeCalls = new GenomicRegionSortedCollection<>(refGenome.getSequencesList());
		GenomicRegionSortedCollection<SimplifiedReadAlignment> sortedAlns = new GenomicRegionSortedCollection<>(refGenome.getSequencesList());
		for (List<SimplifiedReadAlignment> alns : alignments.values()) {
			sortedAlns.addAll(alns);
		}
		sortedAlns.forceSort();
		for(GenomicVariant variant : variantsList) {
			int genotype;
			List<SimplifiedReadAlignment> spanningAlns = sortedAlns.findSpanningRegions(variant).asList();
			int varCoveredAlns = spanningAlns.size();
			if(varCoveredAlns == 0) genotype = CalledGenomicVariant.GENOTYPE_UNDECIDED;
			else {
				genotype = assignBayesianGenotype(variant, spanningAlns);
			}
			if(genotype != CalledGenomicVariant.GENOTYPE_UNDECIDED && 
					genotype != CalledGenomicVariant.GENOTYPE_HOMOREF) {
				CalledGenomicVariant calledVariant = new CalledGenomicVariantImpl(variant, genotype);
				calledVariant.setGenotypeQuality(variant.getVariantQS());
				genotypeCalls.add(calledVariant);
			}
		}
		return genotypeCalls;
	}
	
	private int assignBayesianGenotype(GenomicVariant variant, List<SimplifiedReadAlignment> spanningAlns) {
		// TODO Auto-generated method stub
		int genotype = -1;
		List<GenomicVariant> calls = new ArrayList<>();
		Set<String> readNames = new HashSet<>();
		for(SimplifiedReadAlignment alignment : spanningAlns) {
			computeSpanningAlnCall(alignment, variant, calls, readNames);
		}
		double [][] genotypeProbabilities = calculateCalledVariantGenotypePosteriorProbabilities(variant, calls);
		//from the probabilities choose the most probable one
		genotype = decideGenotype(genotypeProbabilities, variant);
		return genotype;
	}
	
	private boolean alignmentCoversVariant(SimplifiedReadAlignment alignment, GenomicVariant variant) {
		// TODO Auto-generated method stub
		int toleranceDiff;
		int firstDiff = Math.abs(alignment.getFirst() - variant.getFirst());
		int lastDiff = Math.abs(alignment.getLast() - variant.getLast());
		if(GenomicVariant.TYPE_LARGEINS == variant.getType()) {
			toleranceDiff = 200;
			if(firstDiff < toleranceDiff || lastDiff < toleranceDiff) return false;
		}
		else if(GenomicVariant.TYPE_LARGEDEL == variant.getType()){
			toleranceDiff = 1000;
			int lastFirstDiff = alignment.getLast() - variant.getFirst();
			int firstLastDiff = variant.getLast() - alignment.getFirst();
			if(lastFirstDiff < toleranceDiff || firstLastDiff < toleranceDiff ) return false;
		}
		return true;
	}

	private void computeSpanningAlnCall(SimplifiedReadAlignment alignment, GenomicVariant variant,
			List<GenomicVariant> calls, Set<String> readNames) {
		// TODO Auto-generated method stub
		Signature candidateCall;
		String variantId = variant.getId();
		boolean addCall = true;
		boolean interAlnReadWasVisited = false;
		if(alignment.hasCallsForVariant(variantId)) {
			candidateCall = (Signature) alignment.getCallByVariantId(variantId);
		}
		else {
			candidateCall = new Signature(variant.getSequenceName(), variant.getFirst(), variant.getLast(),
					GenomicVariant.TYPE_UNDETERMINED);
			candidateCall.setReadName(alignment.getReadName());
		}
		//Add call to calls depending on many factors
		if(candidateCall.getSignatureType() == Signature.INTERALIGNMENT) 
			interAlnReadWasVisited = !readNames.add(candidateCall.getReadName());
		boolean isRefCall = candidateCall.getType() == GenomicVariant.TYPE_UNDETERMINED;
		if(isRefCall) {
			boolean alignmentCoversVariant = alignmentCoversVariant(alignment, variant);
			if(interAlnReadWasVisited || !alignmentCoversVariant) addCall = false;
		}
		if(addCall) calls.add(candidateCall);
	}
	
	private double[][] calculateCalledVariantGenotypePosteriorProbabilities(GenomicVariant variant,
		List<GenomicVariant> calls) {		
		CountsHelper helper = calculateCountsSV(GENOTYPE_ALLELES, variant.length(), calls, 
				CountsHelper.DEF_MAX_BASE_QS, DEF_HET_RATE, false);
		double[][] probabilities = helper.getPosteriorProbabilities(CountsHelper.DEF_HETEROZYGOSITY_RATE_DIPLOID);
		return probabilities;
	}

	public static CountsHelper calculateCountsSV(String[] genotypeAlleles, int avgLength, List<GenomicVariant> calls,
			byte maxBaseQS, double heterozygousProportion, boolean verbose) {
		// TODO Auto-generated method stub
		CountsHelper helper = new CountsHelper(genotypeAlleles);
		helper.setVerbose(verbose);
		if(maxBaseQS>0) helper.setMaxBaseQS(maxBaseQS);
		helper.setHeterozygousProportion(heterozygousProportion);
		SampleStatistics calcLength = new SampleStatistics();
		for(GenomicVariant call: calls) {
			if(call.getType() == GenomicVariant.TYPE_UNDETERMINED) continue;
			calcLength.update(Math.abs(call.length()));
		}
		double lengthVariance = calcLength.getVariance();
		if(lengthVariance<1) lengthVariance = 1;
		NormalDistribution norm = new NormalDistribution(avgLength, lengthVariance);
		GammaDistribution gamma = new GammaDistribution(DEF_SHAPE_PARAMETER_ERROR_DIST);
		for(GenomicVariant call: calls) {
			String callAllele = call.getType() == GenomicVariant.TYPE_UNDETERMINED ? GENOTYPE_ALLELES[0] : GENOTYPE_ALLELES[1];
			helper.updateCountsSV(callAllele, norm, gamma, call.length());
		}
		return helper;
	}
	
	private int decideGenotype(double[][] genotypeProbabilities, GenomicVariant variant) {
		// TODO Auto-generated method stub
		int genotype = -1;
		double bestGenotypeProb=0;
		int bestGenotypeIdx = 0;
		int z = 0;
		for(int i = 0; i < genotypeProbabilities.length; i++) {
			for(int j = 0; j < genotypeProbabilities.length; j++) {
				double currentGenotypeProb = genotypeProbabilities[i][j];
				if(currentGenotypeProb > bestGenotypeProb) {
					bestGenotypeProb = currentGenotypeProb;
					bestGenotypeIdx = z;
				}
				z++;
			}
		}
		if(bestGenotypeIdx == 1 || bestGenotypeIdx == 2) {
			genotype = CalledGenomicVariant.GENOTYPE_HETERO;
		}
		else if(bestGenotypeIdx == 3) {
			genotype = CalledGenomicVariant.GENOTYPE_HOMOALT;
		}
		else if(bestGenotypeIdx == 0) {
		genotype = CalledGenomicVariant.GENOTYPE_HOMOREF ;
		}
		double qualityScore = genotype == CalledGenomicVariant.GENOTYPE_HETERO ? genotypeProbabilities[0][1] +
				genotypeProbabilities[1][0] : bestGenotypeProb;
		int phredScore = calculatePhredScore(1 - qualityScore);
		variant.setVariantQS((short) phredScore);
		return genotype;
	}

	private int calculatePhredScore(double qualityScore) {
		// TODO Auto-generated method stub
		if(qualityScore == 0) qualityScore = Double.MIN_NORMAL;
		double phredSc = Math.log10(qualityScore)*(-10);
		return (int) phredSc;
	}
	

	private void filterIntersectingVariants(GenomicRegionSortedCollection<CalledGenomicVariant> sortedVariants) {
		// TODO Auto-generated method stub
		int n = sortedVariants.size();
		boolean[] visited = new boolean[n];
		List<CalledGenomicVariant> variantsToKeep = new ArrayList<>();
		List<CalledGenomicVariant> variants = sortedVariants.asList();
		for (int i = 0; i < n; i++) {
			CalledGenomicVariant variant = variants.get(i);
			List<CalledGenomicVariant> spanningVariants = sortedVariants.findSpanningRegions(variant).asList();
			if(!visited[i]) {
				if(spanningVariants.size() < 2) {
					variantsToKeep.add(variant);
					continue;
				}
				Collections.sort(spanningVariants, Comparator.comparingInt(v -> (int) v.getVariantQS()));
				variantsToKeep.add(spanningVariants.get(spanningVariants.size() - 1));
				for(GenomicVariant visitedVariant : spanningVariants) {
					int idx = variants.indexOf(visitedVariant);
					visited[idx] = true;
				}
			}
		}
		sortedVariants.retainAll(variantsToKeep);
	}

	private List<GenomicVariantAnnotation> annotateStructuralVariant(GenomicVariant variant){
		List<GenomicVariantAnnotation> annotations = new ArrayList<>();
		int length = variant.getType() == GenomicVariant.TYPE_LARGEDEL ? variant.length()*(-1) : variant.length();
		GenomicVariantAnnotation svLengthAnnot = new GenomicVariantAnnotation(variant,
				GenomicVariantAnnotation.ATTRIBUTE_SVLEN, length);
		GenomicVariantAnnotation svEndAnnot = new GenomicVariantAnnotation(variant,
				GenomicVariantAnnotation.ATTRIBUTE_END, variant.getLast());
		GenomicVariantAnnotation svTypeAnnot = new GenomicVariantAnnotation(variant,
				GenomicVariantAnnotation.ATTRIBUTE_SVTYPE, GenomicVariantImpl.getVariantTypeName(variant.getType()));
		annotations.add(svLengthAnnot);
		annotations.add(svEndAnnot);
		annotations.add(svTypeAnnot);
		return annotations;
	}
	
	private VCFFileHeader createVCFHeader(String sampleId) throws IOException {
		VCFFileHeader header = VCFFileHeader.makeDefaultEmptyHeader();
		//header.addDefaultSample("C24");
		header.addSample(new Sample(sampleId), true);
		String end = "##INFO=<ID=END,Number=1,Type=Integer,Description=" + "\"End position of the structural variant\"" + ">";
		String svtype = "##INFO=<ID=SVTYPE,Number=1,Type=String,Description="
				+ "\"Type of SV:DEL=Deletion, INS=Insertion, DUP=Duplication, INV=Inversion\"" +">";
		String svlen = "##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=" + "\"Difference in length between REF and ALT alleles\"" + ">";
		header.loadHeaderLine(end);
		header.loadHeaderLine(svtype);
		header.loadHeaderLine(svlen);
		return header;
	}
	
	private List<String> getVariantFilters(GenomicVariant variants){
		List<String> filters = new ArrayList<>();
		filters.add(FILTER_SETTING_MISSING);
		return filters;
	}
	
	private List<VCFRecord> buildRecords(List<CalledGenomicVariant> genotypeCalls, VCFFileHeader header){
		List<VCFRecord> records = new ArrayList<>();
		for(int i = 0; i < genotypeCalls.size(); i++) {
			GenomicVariant variant = genotypeCalls.get(i);
			List<CalledGenomicVariant> calls = new ArrayList<>();
			CalledGenomicVariant call = (CalledGenomicVariantImpl) variant;
			calls.add(call);
			List<GenomicVariantAnnotation> infoFields = annotateStructuralVariant(variant);
			List<String> filters = getVariantFilters(variant);
			VCFRecord record = new VCFRecord(variant, filters,
					 infoFields, VCFRecord.DEF_FORMAT_ARRAY_MINIMAL,  calls, header);
			records.add(record);
		}
		return records;
	}
	
	private void printVCFFile(List<VCFRecord> records,VCFFileHeader header, String file) throws IOException {
		try(PrintStream pr = new PrintStream(file);){
			VCFFileWriter writer = new VCFFileWriter();
			writer.printHeader(header, pr);
			writer.printVCFRecords(records, pr);
		}	
	}
	
	private void setIndelTresholdSize(int SVLength) {
		this.lengthToDefineSVEvent = SVLength;
	}
	
	public void saveVCFResultsFile(List<CalledGenomicVariant> genotypeCalls, String sampleID) throws IOException {
		VCFFileHeader header = createVCFHeader(sampleID);
		List<VCFRecord> records = buildRecords(genotypeCalls, header);
		String saveFile = sampleID + ".variants.vcf";
		printVCFFile(records, header, saveFile);
	}
	
	public List<CalledGenomicVariant> run(String alnFile) throws IOException {
		readAlignments(alnFile, lengthToDefineSVEvent);
		Map<String, List<List<Integer>>> clusters = new LinkedHashMap<>();
		if(MAX_CLIQUE_FINDER_ALGORITHM.equals(clusteringAlgorithm)) {
			LongReadVariantDetectorClusteringAlgorithm caller = new MaxCliqueClusteringDetectionAlgorithm(refGenome, signatures, lengthToDefineSVEvent);
			clusters = caller.callVariantClusters();
		}
		else if(DBSCAN_ALGORITHM.equals(clusteringAlgorithm)) {
			LongReadVariantDetectorClusteringAlgorithm caller = new DBSCANClusteringDetectionAlgorithm(refGenome, signatures, lengthToDefineSVEvent);
			clusters = caller.callVariantClusters();
		}
		else if(STRONGLY_CONNECTED_COMPONENTS_ALGORITHM.equals(clusteringAlgorithm)) {
			LongReadVariantDetectorClusteringAlgorithm caller = new SCCClusteringDetectionAlgorithm(refGenome, signatures, lengthToDefineSVEvent);
			clusters = caller.callVariantClusters();
		}
		variants = callVariants(clusters);
		List<GenomicVariant> variantsList = variants.asList();
		GenomicRegionSortedCollection<CalledGenomicVariant> genotypeCallsCollection = makeBayesianGenotypeCalls(variantsList);
		filterIntersectingVariants(genotypeCallsCollection);
		List<CalledGenomicVariant> genotypeCalls = genotypeCallsCollection.asList();
		return genotypeCalls;
	}

	/** @throws Exception, IOException */
	public static void main(String[] args) throws Exception, IOException {
		// TODO Auto-generated method stub
		LongReadStructuralVariantDetector caller = new LongReadStructuralVariantDetector();
		String algorithm = args[0];
		String refFile = args[1];
		String alnFile = args[2];
		String prefix = args[3];
		int svLength = Integer.parseInt(args[4]);
		caller.setRefGenome(refFile);
		caller.setClusteringAlgorithm(algorithm);
		caller.setIndelTresholdSize(svLength);
		List<CalledGenomicVariant> calledVariants = caller.run(alnFile);
		caller.saveVCFResultsFile(calledVariants, prefix);
	}
	
	class Signature implements GenomicVariant {
		
		public static final byte INTRAALIGNMENT = 0;
		public static final byte INTERALIGNMENT = 1;

		private final String[] SIGNATURE_TYPES = {"INTRAALIGNMENT", "INTERALIGNMENT"};
		private List<String> DEFAULT_ALLELES_BIALLELIC = Arrays.asList(GENOTYPE_ALLELES);
		
		private String id;
		private String sequenceName;
		private String readName;
		private String variantId;
		private int alnIdx;

		private int first;
		private int last;
		private int length;
		private byte signatureType;
		private boolean comesFromSecondaryAlignment;
		
		private List<String> alleles =new ArrayList<String>();
		private boolean negativeStrand = false;
		private short variantQS=0;
		private byte type = GenomicVariant.TYPE_UNDETERMINED;
		
		/**
		 * Creates a biallelic genomic variant with default allele names
		 * @param sequenceName
		 * @param first
		 * @param last
		 * @param type Type of the variant
		 */
		public Signature(String sequenceName, int first, int last, byte type) {
			this.setSequenceName(sequenceName);
			this.setFirst(first);
			this.setLast(last);
			this.setType(type);
			this.setLength(last-first+1);
			alleles = DEFAULT_ALLELES_BIALLELIC;
		}
		

		@Override
		public String getSequenceName() {
			// TODO Auto-generated method stub
			return sequenceName;
		}

		@Override
		public int getFirst() {
			// TODO Auto-generated method stub
			return first;
		}

		@Override
		public int getLast() {
			// TODO Auto-generated method stub
			return last;
		}

		@Override
		public int length() {
			// TODO Auto-generated method stub
			return length;
		}
		
		public byte getSignatureType() {
			return signatureType;
		}
		
		@Override
		public boolean isPositiveStrand() {
			// TODO Auto-generated method stub
			return !negativeStrand;
		}

		@Override
		public boolean isNegativeStrand() {
			// TODO Auto-generated method stub
			return negativeStrand;
		}

		@Override
		public String [] getAlleles() {
			// TODO Auto-generated method stub
			return alleles.toArray(new String[0]);
		}

		@Override
		public String getReference() {
			// TODO Auto-generated method stub
			return null;
		}

		@Override
		public String getId() {
			// TODO Auto-generated method stub
			return this.id;
		}


		@Override
		public short getVariantQS() {
			// TODO Auto-generated method stub
			return this.variantQS;
		}


		@Override
		public boolean isCompatible(GenomicVariant variant) {
			// TODO Auto-generated method stub
			return false;
		}

		@Override
		public boolean isBiallelic() {
			// TODO Auto-generated method stub
			return false;
		}

		@Override
		public boolean isSNV() {
			// TODO Auto-generated method stub
			return false;
		}

		@Override
		public byte getType() {
			// TODO Auto-generated method stub
			return this.type;
		}
		
		public String getSignatureTypeAsString() {
			return SIGNATURE_TYPES[this.signatureType];
		}
		
		
		private void setLength(int length) {
			// TODO Auto-generated method stub
			this.length = length;
		}

		private void setLast(int last) {
			// TODO Auto-generated method stub
			this.last = last;
		}

		private void setFirst(int first) {
			// TODO Auto-generated method stub
			this.first = first;
		}

		private void setSequenceName(String sequenceName) {
			// TODO Auto-generated method stub
			this.sequenceName = sequenceName;
		}
		
		@Override
		public void setId(String id) {
			// TODO Auto-generated method stub
			this.id = id;
		}
		
		@Override
		public void setVariantQS(short qualityScore) {
			// TODO Auto-generated method stub
			this.variantQS = qualityScore;
		}
		
		@Override
		public void setType(byte type) {
			// TODO Auto-generated method stub
			this.type = type;
		}

		public void setSignatureType(byte signatureType) {
			this.signatureType = signatureType;
		}


		public boolean comesFromSecondaryAlignment() {
			return comesFromSecondaryAlignment;
		}


		public void setOriginFromSecondaryAlignment() {
			this.comesFromSecondaryAlignment = true;
		}


		public String getReadName() {
			return readName;
		}


		public void setReadName(String readName) {
			this.readName = readName;
		}


		public String getVariantId() {
			return variantId;
		}


		public void setVariantId(String variantId) {
			this.variantId = variantId;
		}
		
		public int getAlnIdx() {
			return alnIdx;
		}


		public void setAlnIdx(int alnIdx) {
			this.alnIdx = alnIdx;
		}
	}
	
	class SimplifiedReadAlignment implements GenomicRegion{
		
		private int sequenceIndex;
		private String sequenceName;
		private int first;
		private int last;
		private byte alignmentQuality;
		private int readLength;
		private int softClipStart;
		private int softClipEnd;
		private boolean isNegativeStrand;
		private Map<Integer,GenomicVariant> signatureCalls; 
		private Map<String, GenomicVariant> callsPerVariantID;
		private String readName;
		private boolean comesFromSecondaryAlignment;
		
		public SimplifiedReadAlignment(ReadAlignment aln) {
			super();
			this.sequenceIndex = aln.getSequenceIndex();
			this.sequenceName = aln.getSequenceName();
			this.first = aln.getFirst();
			this.last = aln.getLast();
			this.alignmentQuality = aln.getAlignmentQuality();
			this.readLength = aln.getReadLength();
			this.softClipStart = aln.getSoftClipStart();
			this.softClipEnd = aln.getSoftClipEnd();
			//this.indelCalls = aln.getIndelCalls();
			this.readName = aln.getReadName();
			this.isNegativeStrand = aln.isNegativeStrand(); 
			this.setCallsPerVariantID(new HashMap<>());
		}
		
		public boolean hasCallsForVariant(String variantId) {
			return this.callsPerVariantID.containsKey(variantId);
		}
		
		public GenomicVariant getCallByVariantId(String variantId) {
			return this.callsPerVariantID.get(variantId);
		}
		
		public GenomicVariant addCallByVariantId(String variantId, GenomicVariant call) {
			return this.callsPerVariantID.put(variantId, call);
		}
		
		public GenomicVariant getIndelCall(int referencePos) {
			if(signatureCalls == null) return null;
			return signatureCalls.get(referencePos);
		}
		
		public int getSequenceIndex() {
			return sequenceIndex;
		}

		public byte getAlignmentQuality() {
			return alignmentQuality;
		}

		public int getSoftClipStart() {
			return softClipStart;
		}

		public int getSoftClipEnd() {
			return softClipEnd;
		}

		@Override
		public String getSequenceName() {
			// TODO Auto-generated method stub
			return sequenceName;
		}

		@Override
		public int getFirst() {
			// TODO Auto-generated method stub
			return first;
		}

		@Override
		public int getLast() {
			// TODO Auto-generated method stub
			return last;
		}

		@Override
		public int length() {
			// TODO Auto-generated method stub
			return last - first + 1;
		}
		public String getReadName() {
			return readName;
		}
		
		public int getReadLength() {
			return readLength;
		}
		
		@Override
		public boolean isPositiveStrand() {
			// TODO Auto-generated method stub
			return !isNegativeStrand;
		}

		@Override
		public boolean isNegativeStrand() {
			// TODO Auto-generated method stub
			return isNegativeStrand;
		}
		
		public void setIndelCalls(Map<Integer,GenomicVariant> indelCalls) {
			this.signatureCalls = indelCalls;
		}
		
		public Map<Integer,GenomicVariant> getIndelCalls() {
			return signatureCalls;
		}

		public Map<String, GenomicVariant> getCallsPerVariantID() {
			return callsPerVariantID;
		}

		public void setCallsPerVariantID(Map<String, GenomicVariant> callsPerVariantID) {
			this.callsPerVariantID = callsPerVariantID;
		}
		
		public boolean comesFromSecondaryAlignment() {
			return comesFromSecondaryAlignment;
		}

		public void setOriginFromSecondaryAlignment() {
			this.comesFromSecondaryAlignment = true;
		}
	}
}
	
