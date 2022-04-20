package ngsep.discovery;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.concurrent.ThreadLocalRandom;

import JSci.maths.statistics.GammaDistribution;
import JSci.maths.statistics.NormalDistribution;
import JSci.maths.statistics.SampleStatistics;
import ngsep.alignments.ReadAlignment;
import ngsep.alignments.io.ReadAlignmentFileReader;
import ngsep.genome.GenomicRegion;
import ngsep.genome.GenomicRegionComparator;
import ngsep.genome.GenomicRegionImpl;
import ngsep.genome.GenomicRegionSortedCollection;
import ngsep.genome.GenomicRegionSpanComparator;
import ngsep.genome.ReferenceGenome;
import ngsep.graphs.CliquesFinder;
import ngsep.math.Distribution;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;
import ngsep.variants.CalledGenomicVariant;
import ngsep.variants.CalledGenomicVariantImpl;
import ngsep.variants.CalledInversion;
import ngsep.variants.CalledLargeIndel;
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
	
	public static final String[] GENOTYPE_ALLELES = {"REF", "ALT"};
	public static final double DEF_HET_RATE = 0.5;
	public static final double HETEROCIGOZITY_TRESHOLD = 0.8;
	public static final double DEF_SHAPE_PARAMETER_ERROR_DIST = 0.35;
	public static final double DEF_LOG_REF_PROB_SV = Math.log10(0.9);
	public static final double DEF_LOG_ALT_PROB_SV = Math.log10(0.1);
	public static final double DEF_LOG_ERROR_PROB_SV = Math.log10(0.00001);

	
	public static final int INV_DETERMINING_MAX_DISTANCE = 1000;
	public static final int DEL_INTER_DETERMINING_MAX_DISTANCE = 90000;
	public static final String MAX_CLIQUE_FINDER_ALGORITHM = "MCC";
	public static final String DBSCAN_ALGORITHM = "DBSCAN";
	
	private GenomicRegionSortedCollection<GenomicVariant> signatures;
	private GenomicRegionSortedCollection<GenomicVariant> variants;
	//private Map<String, List<GenomicRegion>> alignmentRegions;
	private Map<String, List<SimplifiedReadAlignment>> alignments;
	private ReferenceGenome refGenome;
	private int lengthToDefineSVEvent;
	private int minMQ = DEF_MIN_MQ_UNIQUE_ALIGNMENT;
	
	LongReadStructuralVariantDetector(){
	}
	
	public void setRefGenome(String referenceFile) throws IOException {
		this.refGenome = new ReferenceGenome(referenceFile);
	}
	
	private void readAlignments(String alignmentFile,  int SVLength) throws IOException{
		try(ReadAlignmentFileReader alignmentReader = new ReadAlignmentFileReader(alignmentFile, refGenome)){
			setIndelTresholdSize(SVLength);
			//alignmentReader.setLoadMode(ReadAlignmentFileReader.LOAD_MODE_ALIGNMENT_NAME);
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
					//System.out.println("#Read name=" + alnRegions.get(0).getReadName());
					findInterAlnSignatures(alnRegions);
				}
			}
			signatures.forceSort();
		}
	}
	
	private void computeAlignment(ReadAlignment aln) {
		// TODO Auto-generated method stub
		/**
		GenomicRegion alnRegion = new GenomicRegionImpl(aln.getSequenceName(), aln.getFirst(), aln.getLast());
		List<GenomicRegion> alnRegions = alignmentRegions.computeIfAbsent(aln.getReadName(), v -> new ArrayList<>());
		alnRegions.add(alnRegion);
		**/
		Map<Integer,GenomicVariant> filteredCalls = findIntraAlnSignatures(aln);
		List<SimplifiedReadAlignment> alnRegions = alignments.computeIfAbsent(aln.getReadName(), v -> new ArrayList<>());
		if (!aln.isSecondary()) {
			SimplifiedReadAlignment simpleAln = new SimplifiedReadAlignment(aln);
			simpleAln.setIndelCalls(filteredCalls);
			alnRegions.add(simpleAln);
		}
	}
	
	public GenomicRegionSortedCollection<GenomicVariant> getSignatures() {
		return signatures;
	}
	
	public ReferenceGenome getRefGenome() {
		return refGenome;
	}
	
	private GenomicVariant createIndelInterAlnSignature(String seqName, int first, int last, int length,
			byte type) {
		// TODO Auto-generated method stub
		Signature signature = new Signature(seqName, first, last, type);
		signature.setLength(length);
		signature.setSignatureType(Signature.INTERALIGNMENT);
		return signature;
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
	
	private Map<Integer, GenomicVariant> findIntraAlnSignatures(ReadAlignment aln) {
		Map<Integer, GenomicVariant> calls = aln.getIndelCalls();
		Map<Integer, GenomicVariant> filteredCalls = new HashMap<>();
		/**
		if(aln.getFirst() == 14383838) {
			System.out.println("~secondary=" + aln.isSecondary() + " supplementary=" + aln.isSupplementary());
			System.out.println("~Read alignment=" +  aln.toString());
		}**/
		//System.out.println("~Read alignment=" + aln.getFirst() + " " + aln.getLast() + " "+ aln.length() + " SCStart " + aln.getSoftClipStart() + 
			//	" SCEnd " + aln.getSoftClipEnd());
		if(calls != null){// && !aln.isSecondary()) {
			for (Map.Entry<Integer, GenomicVariant> call : calls.entrySet()) {
				GenomicVariant indel = call.getValue();
				if(indel.length() < lengthToDefineSVEvent) continue;
				Signature sign = (Signature) createIndelIntraAlnSignature(indel);
				if(aln.isSecondary()) sign.setOriginFromSecondaryAlignment();
				filteredCalls.put(call.getKey(), sign);
				signatures.add(sign);
				//System.out.println("secondary=" + aln.isSecondary());
				/**System.out.println(" Sign first: " + sign.getFirst() + " last: " + sign.getLast() + " chr: "
					+ sign.getSequenceName() + " with length: " + sign.length());
				System.out.println(GenomicVariantImpl.getVariantTypeName(sign.getType()));**/
			}
			/**
			if(aln.isUnique()) {
				int softClipStart = aln.getSoftClipStart();
				int softClipEnd = aln.getSoftClipEnd();
				if(softClipStart >= lengthToDefineSVEvent) {
					int first = aln.getFirst();
					int last = first + 1;
					int length = softClipStart;
					String seqName = aln.getSequenceName(); 
					byte type = GenomicVariant.TYPE_LARGEINS;
					GenomicVariantImpl SCFirstSignature = new GenomicVariantImpl(seqName, first, last, type);
					SCFirstSignature.setLength(length);
					filteredCalls.put(first, SCFirstSignature);
					signatures.add(SCFirstSignature);
				}
				if(softClipEnd >= lengthToDefineSVEvent) {
					int first = aln.getLast();
					int last = first + 1;
					int length = softClipEnd;
					String seqName = aln.getSequenceName(); 
					byte type = GenomicVariant.TYPE_LARGEINS;
					GenomicVariantImpl SCLastSignature = new GenomicVariantImpl(seqName, first, last, type);
					SCLastSignature.setLength(length);
					filteredCalls.put(first, SCLastSignature);
					signatures.add(SCLastSignature);
				}
			}**/
		}
		return filteredCalls; 
	}
	
	private void findInterAlnSignatures(List<SimplifiedReadAlignment> alnRegions) {
		// TODO Auto-generated method stub
		//System.out.println("Read aln list: ");
		int nAln = alnRegions.size();
		SimplifiedReadAlignment firstAln = alnRegions.get(0);
		SimplifiedReadAlignment lastAln = alnRegions.get(nAln - 1);
		String sequenceName = firstAln.getSequenceName();
		if(nAln == 2) {
			// && Math.abs(firstAln.getLast() - lastAln.getFirst())
			//Two aln that are very far away could be a DUP
			computeInterAlnIndel(firstAln, lastAln);
				}
		else if(nAln >= 3) {
			SimplifiedReadAlignment midAln = alnRegions.get(1);
			for(int i = 0; i < 2; i++) {
				int j = i + 1;
				SimplifiedReadAlignment currentAln = alnRegions.get(i);
				SimplifiedReadAlignment secondAln = alnRegions.get(j);
				computeInterAlnIndel(currentAln, secondAln);				
			}
			if(sequenceName.equals(midAln.getSequenceName()) && sequenceName.equals(lastAln.getSequenceName())) {
				//A DEL in two aln could be an INS far away
				int firstDistance = Math.abs(midAln.getFirst() - firstAln.getLast());
				int secondDistance = Math.abs(midAln.getLast() - lastAln.getFirst());
				if(secondDistance > INV_DETERMINING_MAX_DISTANCE && firstDistance > INV_DETERMINING_MAX_DISTANCE) {
					int first = midAln.getFirst();
					int last = midAln.getLast();
					int length = last - first + 1;
					byte type = GenomicVariant.TYPE_INVERSION;
					if(length >= lengthToDefineSVEvent) {
						GenomicVariant signature = createIndelInterAlnSignature(sequenceName, first, last, length, type);
						signatures.add(signature);
						//System.out.println("INV=" + signature.getSequenceName() + " " + signature.getFirst() + " " + signature.getLast() + " " +
							//+ signature.length() + " " + GenomicVariantImpl.getVariantTypeName(signature.getType()));
					}
				}
			}
		}
		for(int i = 0; i < alnRegions.size(); i++) {
			SimplifiedReadAlignment current = alnRegions.get(i);
			//System.out.println("#" + current.getSequenceName() + "\t" + current.getFirst() + "\t" + current.getLast() +
				//	"\t" + current.length() + "\tScStart" + current.getSoftClipStart() + "\tScEnd" + current.getSoftClipEnd());
		}
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
	
	private void computeInterAlnIndel(SimplifiedReadAlignment firstAln, SimplifiedReadAlignment lastAln) {
		int distance = Math.abs(lastAln.getFirst() - firstAln.getLast());
		String sequenceName = firstAln.getSequenceName();
		int first = firstAln.getLast() + 1;
		int last = lastAln.getFirst();
		int length = last - first + 1;
		int interAlnDistance = estimateInterAlnLength(firstAln, lastAln);
		if(distance >= lengthToDefineSVEvent && interAlnDistance <= 100) { //distance <= DEL_INTER_DETERMINING_MAX_DISTANCE && getInterAlnDistance(firstAln, lastAln) < lengthToDefineSVEvent) {
			//System.out.println("distance=" + distance);
			byte type = GenomicVariant.TYPE_LARGEDEL;
			if(length >= lengthToDefineSVEvent) {
				GenomicVariant signature = createIndelInterAlnSignature(sequenceName, first, last, length, type);
				signatures.add(signature);
				Map<Integer, GenomicVariant> indelCalls = firstAln.getIndelCalls();
				indelCalls.put(first, signature);
				//System.out.println("DEL=" + signature.getSequenceName() + " " + signature.getFirst() + " " + signature.getLast() + " " +
					//	+ signature.length() + " " + GenomicVariantImpl.getVariantTypeName(signature.getType()));
			}
		}
		else if(interAlnDistance >= lengthToDefineSVEvent && distance <= lengthToDefineSVEvent){
			byte type = GenomicVariant.TYPE_LARGEINS;
			last = first + 1;
			length = interAlnDistance;
			GenomicVariant signature = createIndelInterAlnSignature(sequenceName, first, last, length, type);
			signatures.add(signature);
			GenomicRegionSpanComparator spanCmp = GenomicRegionSpanComparator.getInstance();
			Map<Integer, GenomicVariant> indelCalls;
			if(spanCmp.compare(firstAln, signature) == 0) {
				indelCalls = firstAln.getIndelCalls();
				indelCalls.put(first, signature);
			}
			else {
				indelCalls = lastAln.getIndelCalls();
				indelCalls.put(first, signature);
			}
			//System.out.println("INS=" + signature.getSequenceName() + " " + signature.getFirst() + " " + signature.getLast() + " " +
				//	+ signature.length() + " " + GenomicVariantImpl.getVariantTypeName(signature.getType()));
			//GenomicVariant call = indelCalls.get(first);
			//System.out.println(call);
			//System.out.println("Call=" + call.getSequenceName() + " " + call.getFirst() + " " + call.getLast() + " " +
				//	+ call.length() + " " + GenomicVariantImpl.getVariantTypeName(call.getType()));
		}
	}
	
	/**
	private int getInterAlnDistance(SimplifiedReadAlignment aln1, SimplifiedReadAlignment aln2) {
		System.out.println("aln2start=" + aln2.getSoftClipStart() + " aln1=" + (aln1.getReadLength() - aln1.getSoftClipEnd()) +
				" return="+(aln2.getSoftClipStart() - (aln1.getReadLength() - aln1.getSoftClipEnd())));
		return aln2.getSoftClipStart() - (aln1.getReadLength() - aln1.getSoftClipEnd());
	}**/
	
	@Override
	public GenomicRegionSortedCollection<GenomicVariant> callVariants(Map<String, List<List<Integer>>> clusters){
		List<GenomicVariant> variants = new ArrayList<>();
		List<String> keys = new ArrayList<>(clusters.keySet());
		QualifiedSequenceList sequences = new QualifiedSequenceList(refGenome.getSequencesList());
		GenomicRegionSortedCollection<GenomicVariant> sortedVariants = new GenomicRegionSortedCollection<>(sequences);
		double secondaryAlnExcludingTreshold = 0.5;
		for(String k:keys) {
			List<List<Integer>> chrClusters = clusters.get(k);
			List<GenomicVariant> chrSignList = signatures.getSequenceRegions(k).asList();
			for(List<Integer> cluster:chrClusters) {
				List<GenomicVariant> clusterSigns = new ArrayList<>();
				boolean comesFromSecondaryAln = false;
				int secondaryAlnSignCount = 0;
				for(int idx:cluster) {
					Signature clusteredSign = (Signature) chrSignList.get(idx);
					if(clusteredSign.comesFromSecondaryAlignment()) secondaryAlnSignCount++; 
					clusterSigns.add(clusteredSign);
				}
				/**boolean trace = false;
				System.out.println("$Cluster with size: " + clusterSigns.size());
				for (GenomicVariant sign:clusterSigns) {
					System.out.println("$Signature in cluster: " + sign.getSequenceName() + " begin: " + sign.getFirst() + " end: " +
							sign.getLast() + " length: " + sign.length() + " type: " + 
							GenomicVariantImpl.getVariantTypeName(sign.getType()));
					if(sign.getFirst() == 1716430) trace = true;
				}**/
				double secondaryAlnSignPercentage = (double) secondaryAlnSignCount/cluster.size();
				if(secondaryAlnSignPercentage >= secondaryAlnExcludingTreshold) comesFromSecondaryAln = true; 
				//if(trace) System.out.println("comes from secondary " + comesFromSecondaryAln);
				if(comesFromSecondaryAln) continue;
				Collections.sort(clusterSigns, Comparator.comparingInt(s -> s.getFirst()));
				GenomicVariantImpl variant = processClusterToVariant(clusterSigns, k);
				//if(trace) System.out.println("Variant created: " + variant.getSequenceName() + variant.getFirst() + variant.getLast() +
					//	variant.length());
				//if(variant.getFirst() == 14363746 && variant.getLast() == 14384751 && Math.abs(variant.length()) == 21006
					//	&& variant.getType() == GenomicVariant.TYPE_LARGEDEL) {
				//}
				if(variant.length() >= lengthToDefineSVEvent) variants.add(variant);
			}
		}
		sortedVariants.addAll(variants);
		sortedVariants.forceSort();
		/**List<GenomicVariant> test = sortedVariants.findSpanningRegions("chr1" , 1716430, 1716440).asList();
		for(GenomicVariant gv : test) {
			System.out.println("$variant in region before filter: " + gv.getSequenceName() + " begin: " + gv.getFirst() + " end: " +
					gv.getLast() + " length: " + gv.length() + " type: " + 
					GenomicVariantImpl.getVariantTypeName(gv.getType()));
		}
		List<GenomicVariant> test2 = sortedVariants.findSpanningRegions("chr1" , 17536260, 17561570).asList();
		for(GenomicVariant gv : test2) {
			System.out.println("$variant in region after filter: " + gv.getSequenceName() + " begin: " + gv.getFirst() + " end: " +
					gv.getLast() + " length: " + gv.length() + " type: " + 
					GenomicVariantImpl.getVariantTypeName(gv.getType()));
		}**/
		//System.out.println("Variants filtered");
		return sortedVariants;
	}
	
	private void filterIntersectingVariants(GenomicRegionSortedCollection<GenomicVariant> sortedVariants) {
		// TODO Auto-generated method stub
		int n = sortedVariants.size();
		boolean[] visited = new boolean[n];
		List<GenomicVariant> variantsToKeep = new ArrayList<>();
		List<GenomicVariant> variants = sortedVariants.asList();
		for (int i = 0; i < n; i++) {
			GenomicVariant variant = variants.get(i);
			List<GenomicVariant> spanningVariants = sortedVariants.findSpanningRegions(variant).asList();
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

	private GenomicVariantImpl processClusterToVariant(List<GenomicVariant> clusterSigns, String sequenceName){
		GenomicVariant firstVar = clusterSigns.get(0);
		SampleStatistics calcFirst = new SampleStatistics();
		SampleStatistics calcLast = new SampleStatistics();
		//detectAndRemoveOutliers(clusterSigns);
		int endOfSpan = -1;
		for(GenomicVariant v:clusterSigns) {
			int currentEndOfSpan = v.getFirst() + v.length() - 1;
			calcFirst.update(v.getFirst());
			calcLast.update(currentEndOfSpan);
			//if(currentEndOfSpan > endOfSpan) endOfSpan = currentEndOfSpan;
		}
		int first = (int) calcFirst.getMean();
		int last = (int) calcLast.getMean();
		endOfSpan = last;
		//System.out.println("first: " + first);
		//System.out.println("last: " + last);
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
		short variantScore = calculateVariantScore(clusterSigns, variant.length());
		variant.setVariantQS(variantScore);
		/**if(variant.getFirst() == 1716434 && variant.getLast() == 1716435 && variant.length() == 303) {
			System.out.println("$WEIRD variant processed from cluster: " + variant.getSequenceName() + " begin: " + variant.getFirst() + " end: " +
					variant.getLast() + " length: " + variant.length() + " type: " + 
					GenomicVariantImpl.getVariantTypeName(variant.getType()));
			for (GenomicVariant sign:clusterSigns) {
				System.out.println("$WEIRD Signature in cluster: " + sign.getSequenceName() + " begin: " + sign.getFirst() + " end: " +
						sign.getLast() + " length: " + sign.length() + " type: " + 
						GenomicVariantImpl.getVariantTypeName(sign.getType()));
			}
		}**/
		return variant;
	}
	/**
	public void filterVariantsBySimilarity(GenomicRegionSortedCollection<GenomicVariant> sortedVariants) {
		Set<GenomicVariant> variantsToKeep = new HashSet<>();
		List<GenomicVariant> variantsList = sortedVariants.asList();
		List<List<GenomicVariant>> variantsPerType = new ArrayList<>();
		GenomicRegionComparator cmpClassInstance = new GenomicRegionComparator(sortedVariants.getSequenceNames());
		List<GenomicVariant> indel = new ArrayList<>();
		List<GenomicVariant> del = new ArrayList<>();
		List<GenomicVariant> ins = new ArrayList<>();
		List<GenomicVariant> inv = new ArrayList<>();
		List<GenomicVariant> dup = new ArrayList<>();
		List<GenomicVariant> tan = new ArrayList<>();
		List<GenomicVariant> und = new ArrayList<>();
		for(int i = 0; i < sortedVariants.size(); i++) {
			GenomicVariant current = variantsList.get(i);
			byte type = current.getType();
			if(GenomicVariant.TYPE_INDEL == type) indel.add(current);
			else if(GenomicVariant.TYPE_LARGEDEL == type) del.add(current);
			else if(GenomicVariant.TYPE_LARGEINS == type) ins.add(current);
			else if(GenomicVariant.TYPE_DUPLICATION  == type) dup.add(current);
			else if(GenomicVariant.TYPE_INVERSION  == type) inv.add(current);
			else if(GenomicVariant.TYPE_CNV == type) tan.add(current);
			else if(GenomicVariant.TYPE_UNDETERMINED == type) und.add(current);
		}
		variantsPerType.add(ins);
		variantsPerType.add(del);
		variantsPerType.add(dup);
		variantsPerType.add(inv);
		variantsPerType.add(tan);
		variantsPerType.add(indel);
		variantsPerType.add(und);
		for(List<GenomicVariant> typedVariantList : variantsPerType) {
			if(typedVariantList.isEmpty()) continue;
			for(int i = 0; i < typedVariantList.size() - 1; i++) {
				GenomicVariant current = typedVariantList.get(i);
				GenomicVariant next = typedVariantList.get(i+1);
				int cmp;
				if(current.getType() == GenomicVariant.TYPE_LARGEINS) {
					cmp = GenomicVariantImpl.testCompatibility(current, next, true, false) ? 0 : 10;
				}
				else cmp = cmpClassInstance.compare(current, next);
				if(cmp == -1 || cmp == 0 || cmp == 1) {
					variantsToKeep.add(current.getVariantQS() > next.getVariantQS() ? current : next);
				}
				else variantsToKeep.add(current);
			}
		}
		sortedVariants.retainAll(variantsToKeep);
		sortedVariants.forceSort();
	}**/
	
	/**
	 * Calculates the likelihood of signatures given genomic variant hypotheses
	 * @param candidates Signature candidates
	 * @param variantHypothesis minimum hypotheses are no-variant or 1 variant.L(D | NoVar) = 1 - return
	 * @return Likelihood of data given the variant hypothesis. 
	 */
	/**
	public double calculateSignatureLikelihood(List<GenomicVariant> candidates, GenomicVariant variantHypothesis) {
		double[] likelihoods;
		int n = candidates.size();
		Collections.sort(candidates, Comparator.comparingInt(s -> s.length()));
		Distribution lengthDistribution = new Distribution(candidates.get(0).getFirst(), candidates.get(n-1).getFirst(),
				1);
		Collections.sort(candidates, Comparator.comparingInt(s -> s.getFirst()));
		Distribution positionDistribution = new Distribution(candidates.get(0).getFirst(), candidates.get(n-1).getFirst(),
				1);
		for(GenomicVariant sign:candidates) {
			calcPos.update(sign.getFirst());
			calcSpan.update(sign.length());
		}
		double stdSpan = Math.sqrt(calcSpan.getVariance());
		double stdPos = Math.sqrt(calcPos.getVariance());
		
		return null;
	}
	**/
	
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
	}
/**
	private void detectAndRemoveOutliers(List<GenomicVariant> clusterSigns) {
		// TODO Auto-generated method stub
		List<GenomicVariant> outliersToRemove = new ArrayList<>();
		int n = clusterSigns.size();
		double [] pdValues = new double[n];
		SampleStatistics clusterCalc = new SampleStatistics();
		double spd;
		SampleStatistics signCalc = new SampleStatistics();
		for(int i = 0; i < n; i++) {
			GenomicVariant si = clusterSigns.get(i);
			signCalc = new SampleStatistics();
			for(int j = 0; j < n; j++) {
				GenomicVariant sj = clusterSigns.get(j);
				if(si.equals(sj)) continue;
				else spd = MaxCliqueClusteringDetectionAlgorithm.calculateSD(si.length(), sj.length());
				signCalc.update(spd);
			}
			double rowAvg = signCalc.getMean();
			clusterCalc.update(rowAvg);
			pdValues[i] = rowAvg;
		}
		double averageSPD = clusterCalc.getMean();
		double stdSPD = Math.sqrt(clusterCalc.getVariance());
		double limit = 2*stdSPD;
		for(int i = 0; i < n; i++) {
			double distance = Math.abs(pdValues[i] - averageSPD);
			if(distance > limit) {
				outliersToRemove.add(clusterSigns.get(i));
			}
		}
		clusterSigns.removeAll(outliersToRemove);
	}
**/
	private GenomicRegionSortedCollection<GenomicVariant> makeGenotypeCalls(List<GenomicVariant> variantsList) {
		GenomicRegionSortedCollection<GenomicVariant> genotypeCalls = new GenomicRegionSortedCollection<>(refGenome.getSequencesList());
		//GenomicRegionComparator cmpClassInstance = new GenomicRegionComparator(variants.getSequenceNames());
		GenomicRegionSortedCollection<SimplifiedReadAlignment> sortedAlns = new GenomicRegionSortedCollection<>(refGenome.getSequencesList());
		for (List<SimplifiedReadAlignment> alns : alignments.values()) {
			sortedAlns.addAll(alns);
		}
		sortedAlns.forceSort();
		for(GenomicVariant variant : variantsList) {
			int genotype;
			int varCalls = 0;
			//GenomicVariantImpl dummyVariant = (GenomicVariantImpl) variant;
			//dummyVariant.setFirst(variant.getFirst() - 1000);
			//System.out.println("*Variant: " + variant.getSequenceName() + " begin: " + variant.getFirst() + " end: " +
				//	variant.getLast() + " length: " + variant.length() + " type: " + 
					//GenomicVariantImpl.getVariantTypeName(variant.getType()));
			List<SimplifiedReadAlignment> spanningAlns = sortedAlns.findSpanningRegions(variant).asList();
			int varCoveredAlns = spanningAlns.size();
			if(varCoveredAlns == 0) genotype = -1;
			else {
				for(SimplifiedReadAlignment alignment : spanningAlns) {
					//System.out.println("*SpanningSRA: " + alignment.getSequenceName() + " begin: " + alignment.getFirst() + " end: " +
						//	alignment.getLast() + " length: " + alignment.length());
					if(hasSVAlleleAlnCall(alignment, variant)) varCalls++;
				}
				//Check calls for intraalignment evidence
				genotype = decideGenotype(varCalls, varCoveredAlns);
				//System.out.println("genotype=" + genotype);
			}
			CalledGenomicVariant calledVariant = new CalledGenomicVariantImpl(variant, genotype);
			genotypeCalls.add(calledVariant);
		}
		return genotypeCalls;
	}
	
	private GenomicRegionSortedCollection<GenomicVariant> makeBayesianGenotypeCalls(List<GenomicVariant> variantsList) {
		GenomicRegionSortedCollection<GenomicVariant> genotypeCalls = new GenomicRegionSortedCollection<>(refGenome.getSequencesList());
		//GenomicRegionComparator cmpClassInstance = new GenomicRegionComparator(variants.getSequenceNames());
		GenomicRegionSortedCollection<SimplifiedReadAlignment> sortedAlns = new GenomicRegionSortedCollection<>(refGenome.getSequencesList());
		for (List<SimplifiedReadAlignment> alns : alignments.values()) {
			sortedAlns.addAll(alns);
		}
		sortedAlns.forceSort();
		for(GenomicVariant variant : variantsList) {
			int genotype;
			//int varCalls = 0;
			//GenomicVariantImpl dummyVariant = (GenomicVariantImpl) variant;
			//dummyVariant.setFirst(variant.getFirst() - 1000);
			//System.out.println("*Variant: " + variant.getSequenceName() + " begin: " + variant.getFirst() + " end: " +
				//	variant.getLast() + " length: " + variant.length() + " type: " + 
					//GenomicVariantImpl.getVariantTypeName(variant.getType()));
			List<SimplifiedReadAlignment> spanningAlns = sortedAlns.findSpanningRegions(variant).asList();
			int varCoveredAlns = spanningAlns.size();
			if(varCoveredAlns == 0) genotype = CalledGenomicVariant.GENOTYPE_UNDECIDED;
			else {
				genotype = assignBayesianGenotype(variant, spanningAlns);
				//System.out.println("genotype=" + genotype);
			}
			CalledGenomicVariant calledVariant = new CalledGenomicVariantImpl(variant, genotype);
			genotypeCalls.add(calledVariant);
		}
		//return genotypeCalls;
		return genotypeCalls;
	}
	
	private int assignBayesianGenotype(GenomicVariant variant, List<SimplifiedReadAlignment> spanningAlns) {
		// TODO Auto-generated method stub
		int genotype = -1;
		List<GenomicVariant> calls = new ArrayList<>();
		for(SimplifiedReadAlignment alignment : spanningAlns) {
			//System.out.println("*SpanningSRA: " + alignment.getSequenceName() + " begin: " + alignment.getFirst() + " end: " +
				//	alignment.getLast() + " length: " + alignment.length());
			GenomicVariant call = computeSpanningAlnCall(alignment, variant);
			calls.add(call);
		}
		double [][] genotypeProbabilities = calculateCalledVariantGenotypePosteriorProbabilities(variant, calls);
		//from the probabilities choose the most probable one
		genotype = decideGenotype(genotypeProbabilities, variant);
		return genotype;
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
		qualityScore = qualityScore*100;
		variant.setVariantQS((short) qualityScore);
		return genotype;
	}

	private GenomicVariant computeSpanningAlnCall(SimplifiedReadAlignment alignment, GenomicVariant variant) {
		// TODO Auto-generated method stub
		GenomicVariant candidateIndelCall;
		int varFirst = variant.getFirst();
		int varLength = variant.length();
		int baseQueueMidWindowSize = (int) (varLength * 0.20);
		int baseQueueBegin = varFirst - baseQueueMidWindowSize;
		int baseQueueLimit = varFirst + baseQueueMidWindowSize;
		for(int i = baseQueueBegin; i <= baseQueueLimit; i++) {
			candidateIndelCall = alignment.getIndelCall(i);
			byte firstIndelCallType;
			if(candidateIndelCall != null) {
				firstIndelCallType = candidateIndelCall.getFirst() + 1 == candidateIndelCall.getLast() ? 
						GenomicVariant.TYPE_LARGEINS : GenomicVariant.TYPE_LARGEDEL;
				if(firstIndelCallType == variant.getType()) {
					candidateIndelCall.setType(firstIndelCallType);
					return candidateIndelCall;
				}
			}
		}
		candidateIndelCall = new GenomicVariantImpl(variant.getSequenceName(), varFirst, variant.getLast(),
				GenomicVariant.TYPE_UNDETERMINED);
		return candidateIndelCall;
	}

	private double[][] calculateCalledVariantGenotypePosteriorProbabilities(GenomicVariant variant,
		List<GenomicVariant> calls) {		
		CountsHelper helper = calculateCountsSV(GENOTYPE_ALLELES, variant.length(), calls, 
				CountsHelper.DEF_MAX_BASE_QS, DEF_HET_RATE, true);
		double[][] probabilities = helper.getPosteriorProbabilities(CountsHelper.DEF_HETEROZYGOSITY_RATE_DIPLOID);
		System.out.println("Posterior probabilities for var: first=" + variant.getFirst() + " last=" + variant.getLast() + 
				" length" + variant.length() + " type" + GenomicVariantImpl.getVariantTypeName(variant.getType()));
		helper.printProbs(probabilities, false);
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

	/**Method implemented by reiterating the alignments in the bam file
	private List<CalledGenomicVariant> makeGenotypeCalls(String alignmentFile, List<GenomicVariant> variantsList) throws IOException{
		List<CalledGenomicVariant> genotypeCalls = new ArrayList<>();
		try(ReadAlignmentFileReader alignmentReader = new ReadAlignmentFileReader(alignmentFile, refGenome)){
			int filterFlags = ReadAlignment.FLAG_READ_UNMAPPED;
			alignmentReader.setFilterFlags(filterFlags);
			alignmentReader.setMinMQ(minMQ);
			GenomicRegionComparator cmpClassInstance = new GenomicRegionComparator(variants.getSequenceNames());
			Iterator<ReadAlignment> it = alignmentReader.iterator();
			ReadAlignment aln = it.next();	
			int varIdx = 0;
			int genotype;
			int currIntersectingVars = 0;
			int [] alnCalls = new int [variantsList.size()];
			int [] alnCoveredPerVar = new int [variantsList.size()];
			boolean lastVar = false;
			while(it.hasNext()) {
				GenomicVariant var = variantsList.get(varIdx);
				//System.out.println("var=" + var.getSequenceName() + " " + var.getFirst() 
				//+ " " + var.getLast() + " "  + var.length());
				//System.out.println("aln=" + aln.getSequenceName() + " " + aln.getFirst() 
				//+ " " + aln.getLast() + " "  + aln.length());
				int cmp;
				//if(aln.getLast() == aln.getFirst()) cmp = modifiedAlnComparison(var, aln, cmpClassInstance);
				//else cmp = cmpClassInstance.compare(var, aln);
				cmp = cmpClassInstance.compare(var, aln);
				//System.out.println("cmp=" + cmp);
				if(cmp > 1 || lastVar) {
					varIdx -= currIntersectingVars;
					currIntersectingVars = 0;
					if(it.hasNext()) aln = it.next();
					lastVar = false;
				}
				else if(cmp < -1) {
					genotype = decideGenotype(alnCalls[varIdx], alnCoveredPerVar[varIdx]);
					//System.out.println("calls% " + alnCalls[varIdx]);
					//System.out.println("all% " + alnCoveredPerVar[varIdx]);
					//System.out.println("het% " + (double) alnCalls[varIdx] / alnCoveredPerVar[varIdx]);
					CalledGenomicVariant calledVariant = new CalledGenomicVariantImpl(var, genotype);
					genotypeCalls.add(calledVariant);
					//System.out.println("genotype " + genotype);
					if(varIdx + 1 != variantsList.size()) varIdx++;
					else lastVar = true;
				}
				else {
					alnCoveredPerVar[varIdx]++;
					//Implement something for interaln calls
					if(hasSVAlleleAlnCall(aln, var)) alnCalls[varIdx]++;
					currIntersectingVars++;
					if(varIdx + 1 != variantsList.size()) varIdx++;
					else lastVar = true;
				}
			}
			//add remaining calls
			for(;varIdx < variantsList.size(); varIdx++) {
				GenomicVariant var = variantsList.get(varIdx);	
				genotype = decideGenotype(alnCalls[varIdx], alnCoveredPerVar[varIdx]);
				CalledGenomicVariant calledVariant = new CalledGenomicVariantImpl(var, genotype);
				genotypeCalls.add(calledVariant);
				varIdx++;
			}
		}
		return genotypeCalls;
	}
	**/
	private boolean hasSVAlleleAlnCall(SimplifiedReadAlignment aln, GenomicVariant var) {
		int varFirst = var.getFirst();
		int varLength = var.length();
		int baseQueueMidWindowSize = (int) (varLength * 0.20);
		int baseQueueBegin = varFirst - baseQueueMidWindowSize;
		int baseQueueLimit = varFirst + baseQueueMidWindowSize;
		int lengthDiffLimit = (int) (varLength * 0.40);
		for(int i = baseQueueBegin; i <= baseQueueLimit; i++) {
			GenomicVariant candidateIndelCall = aln.getIndelCall(i);
			byte firstIndelCallType;
			if(candidateIndelCall != null) {
				firstIndelCallType = candidateIndelCall.getFirst() + 1 == candidateIndelCall.getLast() ? 
						GenomicVariant.TYPE_LARGEINS : GenomicVariant.TYPE_LARGEDEL;
				//System.out.println("call=" + candidateIndelCall.getSequenceName() + " " + candidateIndelCall.getFirst() 
					//	+ " " + candidateIndelCall.getLast() + " "  + candidateIndelCall.length());
				if(firstIndelCallType == var.getType()) {
					int callLength = candidateIndelCall.length();
					int lengthDiff = Math.abs(callLength - varLength);
					if(lengthDiff <= lengthDiffLimit) return true;
				}
			}
		}
		/**
		if(aln.getIndelCalls() != null) {
			for(Map.Entry<Integer, GenomicVariant> entry : aln.getIndelCalls().entrySet()) {
				GenomicVariant tv = entry.getValue();
				System.out.println("&sign=" + tv.getSequenceName() + " " + tv.getFirst() 
				+ " " + tv.getLast() + " "  + tv.length());
			}
		}**/
		return false;
	}
	
	private int decideGenotype(int calls, int coverage) {
		double filterCallTreshold = 0.0;
		double percentage = (double) calls / coverage;
		return percentage > HETEROCIGOZITY_TRESHOLD ? 2 : (percentage > filterCallTreshold ? 1 : -1);
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
	
	private List<VCFRecord> buildRecords(List<GenomicVariant> genotypeCalls, VCFFileHeader header){
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
	
	public void run(String algorithm, String refFile, String alnFile, String sampleId, int SVLength) throws IOException {
		setRefGenome(refFile);
		readAlignments(alnFile, SVLength);
		Map<String, List<List<Integer>>> clusters = new LinkedHashMap<>();
		if(MAX_CLIQUE_FINDER_ALGORITHM.equals(algorithm)) {
			LongReadVariantDetectorClusteringAlgorithm caller = new MaxCliqueClusteringDetectionAlgorithm(refGenome, signatures, SVLength);
			clusters = caller.callVariantClusters();
		}
		else if(DBSCAN_ALGORITHM.equals(algorithm)) {
			LongReadVariantDetectorClusteringAlgorithm caller = new DBSCANClusteringDetectionAlgorithm(refGenome, signatures, SVLength);
			clusters = caller.callVariantClusters();
		}
		variants = callVariants(clusters);
		List<GenomicVariant> variantsList = variants.asList();
		//List<CalledGenomicVariant> genotypeCalls = makeGenotypeCalls(variantsList);
		//GenomicRegionSortedCollection<GenomicVariant> genotypeCallsCollection = makeGenotypeCalls(variantsList);
		GenomicRegionSortedCollection<GenomicVariant> genotypeCallsCollection = makeBayesianGenotypeCalls(variantsList);
		filterIntersectingVariants(genotypeCallsCollection);
		List<GenomicVariant> genotypeCalls = genotypeCallsCollection.asList();
		VCFFileHeader header = createVCFHeader(sampleId);
		//filterVariantsWithoutGenotype(variantsList, genotypeCalls);
		List<VCFRecord> records = buildRecords(genotypeCalls, header);
		String saveFile = sampleId + ".variants.vcf";
		printVCFFile(records, header, saveFile);
	}
	
	private void filterVariantsWithoutGenotype(List<GenomicVariant> variantsList,
			List<CalledGenomicVariant> genotypeCalls) {
		// TODO Auto-generated method stub
		List<GenomicVariant> variantsToRemove = new ArrayList<>();
		List<CalledGenomicVariant> callsToRemove = new ArrayList<>();
		for (int i = 0; i < variantsList.size(); i++) {
			CalledGenomicVariant call = genotypeCalls.get(i);
			GenomicVariant variant = variantsList.get(i);
			if(call.isUndecided()) {
				callsToRemove.add(call);
				variantsToRemove.add(variant);
			}
		}
		genotypeCalls.removeAll(callsToRemove);
		variantsList.removeAll(variantsToRemove);
	}

	/**
	public int modifiedAlnComparison(GenomicVariant var, ReadAlignment aln, GenomicRegionComparator cmpClassInstance) {
		String seqName = aln.getSequenceName();
		int first =  aln.getFirst();
		int length = aln.length();
		int last =  aln.getFirst() + length;
		ReadAlignment modAln = new ReadAlignment(seqName, first, last, length, aln.getFlags());
		return cmpClassInstance.compare(var, modAln);
	}*/
	
	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub
		LongReadStructuralVariantDetector caller = new LongReadStructuralVariantDetector();
		String algorithm = args[0];
		String refFile = args[1];
		String alnFile = args[2];
		String prefix = args[3];
		int svLength = Integer.parseInt(args[4]);
		caller.run(algorithm, refFile, alnFile, prefix, svLength);
	}
	
	class Signature implements GenomicVariant {
		
		public static final byte INTRAALIGNMENT = 0;
		public static final byte INTERALIGNMENT = 1;

		private final String[] SIGNATURE_TYPES = {"INTRAALIGNMENT", "INTERALIGNMENT"};
		private List<String> DEFAULT_ALLELES_BIALLELIC = Arrays.asList(GENOTYPE_ALLELES);
		
		private String id;
		private String sequenceName;
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
			return SIGNATURE_TYPES[this.type];
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
		private Map<Integer,GenomicVariant> indelCalls; 
		private String readName;
		
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
		}
		
		public GenomicVariant getIndelCall(int referencePos) {
			if(indelCalls == null) return null;
			return indelCalls.get(referencePos);
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
			return false;
		}

		@Override
		public boolean isNegativeStrand() {
			// TODO Auto-generated method stub
			return false;
		}
		
		public void setIndelCalls(Map<Integer,GenomicVariant> indelCalls) {
			this.indelCalls = indelCalls;
		}
		
		public Map<Integer,GenomicVariant> getIndelCalls() {
			return indelCalls;
		}
	}
}
	
