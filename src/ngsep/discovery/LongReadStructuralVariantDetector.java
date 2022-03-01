package ngsep.discovery;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
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
	public static final double HETEROCIGOZITY_TRESHOLD = 0.8;
	public static final int INV_DETERMINING_MAX_DISTANCE = 1000;
	public static final int DEL_INTER_DETERMINING_MAX_DISTANCE = 90000;
	public static final String MAX_CLIQUE_FINDER_ALGORITHM = "MCC";
	public static final String DBSCAN_ALGORITHM = "DBSCAN";
	public static final String FILTER_SETTING_MISSING = ".";
	public static final String KEY_SEPARATOR = ",";
	
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
		GenomicVariantImpl signature = new GenomicVariantImpl(seqName, first, last, type);
		signature.setLength(length);
		return signature;
	}
	
	private GenomicVariantImpl createIndelIntraAlnSignature(GenomicVariant indel) {
		int first = indel.getFirst();
		int last = indel.getLast();
		int length = indel.length();
		String seqName = indel.getSequenceName(); 
		byte type = first + 1 == last ? GenomicVariant.TYPE_LARGEINS : GenomicVariant.TYPE_LARGEDEL;
		GenomicVariantImpl signature = new GenomicVariantImpl(seqName, first, last, type);
		signature.setLength(length);
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
				GenomicVariantImpl sign = createIndelIntraAlnSignature(indel);
				if(aln.isSecondary()) sign.setNegativeStrand(true);
				filteredCalls.put(call.getKey(), sign);
				signatures.add(sign);
				//System.out.println("secondary=" + aln.isSecondary());
				System.out.println(" Sign first: " + sign.getFirst() + " last: " + sign.getLast() + " chr: "
					+ sign.getSequenceName() + " with length: " + sign.length());
				System.out.println(GenomicVariantImpl.getVariantTypeName(sign.getType()));
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
						System.out.println("INV=" + signature.getSequenceName() + " " + signature.getFirst() + " " + signature.getLast() + " " +
							+ signature.length() + " " + GenomicVariantImpl.getVariantTypeName(signature.getType()));
					}
				}
			}
		}
		
		for(int i = 0; i < alnRegions.size(); i++) {
			SimplifiedReadAlignment current = alnRegions.get(i);
			System.out.println("#" + current.getSequenceName() + "\t" + current.getFirst() + "\t" + current.getLast() +
					"\t" + current.length() + "\tScStart" + current.getSoftClipStart() + "\tScEnd" + current.getSoftClipEnd());
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
			System.out.println("distance=" + distance);
			byte type = GenomicVariant.TYPE_LARGEDEL;
			if(length >= lengthToDefineSVEvent) {
				GenomicVariant signature = createIndelInterAlnSignature(sequenceName, first, last, length, type);
				signatures.add(signature);
				Map<Integer, GenomicVariant> indelCalls = firstAln.getIndelCalls();
				indelCalls.put(first, signature);
				System.out.println("DEL=" + signature.getSequenceName() + " " + signature.getFirst() + " " + signature.getLast() + " " +
						+ signature.length() + " " + GenomicVariantImpl.getVariantTypeName(signature.getType()));
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
			System.out.println("INS=" + signature.getSequenceName() + " " + signature.getFirst() + " " + signature.getLast() + " " +
					+ signature.length() + " " + GenomicVariantImpl.getVariantTypeName(signature.getType()));
			GenomicVariant call = indelCalls.get(first);
			System.out.println(call);
			System.out.println("Call=" + call.getSequenceName() + " " + call.getFirst() + " " + call.getLast() + " " +
					+ call.length() + " " + GenomicVariantImpl.getVariantTypeName(call.getType()));
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
		double secondaryAlnExcludingTreshold = 0.2;
		for(String k:keys) {
			List<List<Integer>> chrClusters = clusters.get(k);
			List<GenomicVariant> chrSignList = signatures.getSequenceRegions(k).asList();
			for(List<Integer> cluster:chrClusters) {
				List<GenomicVariant> clusterSigns = new ArrayList<>();
				boolean comesFromSecondaryAln = false;
				int secondaryAlnSignCount = 0;
				for(int idx:cluster) {
					GenomicVariant clusteredSign = chrSignList.get(idx);
					if(clusteredSign.isNegativeStrand()) secondaryAlnSignCount++; 
					clusterSigns.add(clusteredSign);
				}
				boolean trace = false;
				System.out.println("$Cluster with size: " + clusterSigns.size());
				for (GenomicVariant sign:clusterSigns) {
					System.out.println("$Signature in cluster: " + sign.getSequenceName() + " begin: " + sign.getFirst() + " end: " +
							sign.getLast() + " length: " + sign.length() + " type: " + 
							GenomicVariantImpl.getVariantTypeName(sign.getType()));
					if(sign.getFirst() == 1716430) trace = true;
				}
				double secondaryAlnSignPercentage = (double) secondaryAlnSignCount/cluster.size();
				if(secondaryAlnSignPercentage >= secondaryAlnExcludingTreshold) comesFromSecondaryAln = true; 
				if(trace) System.out.println("comes from secondary " + comesFromSecondaryAln);
				if(comesFromSecondaryAln) continue;
				Collections.sort(clusterSigns, Comparator.comparingInt(s -> s.getFirst()));
				GenomicVariantImpl variant = processClusterToVariant(clusterSigns, k);
				if(trace) System.out.println("Variant created: " + variant.getSequenceName() + variant.getFirst() + variant.getLast() +
						variant.length());
				//if(variant.getFirst() == 14363746 && variant.getLast() == 14384751 && Math.abs(variant.length()) == 21006
					//	&& variant.getType() == GenomicVariant.TYPE_LARGEDEL) {
				//}
				if(variant.length() >= lengthToDefineSVEvent) variants.add(variant);
			}
		}
		sortedVariants.addAll(variants);
		sortedVariants.forceSort();
		List<GenomicVariant> test = sortedVariants.findSpanningRegions("chr1" , 1716430, 1716440).asList();
		for(GenomicVariant gv : test) {
			System.out.println("$variant in region before filter: " + gv.getSequenceName() + " begin: " + gv.getFirst() + " end: " +
					gv.getLast() + " length: " + gv.length() + " type: " + 
					GenomicVariantImpl.getVariantTypeName(gv.getType()));
		}
		filterIntersectingVariants(sortedVariants);
		List<GenomicVariant> test2 = sortedVariants.findSpanningRegions("chr1" , 17536260, 17561570).asList();
		for(GenomicVariant gv : test2) {
			System.out.println("$variant in region after filter: " + gv.getSequenceName() + " begin: " + gv.getFirst() + " end: " +
					gv.getLast() + " length: " + gv.length() + " type: " + 
					GenomicVariantImpl.getVariantTypeName(gv.getType()));
		}
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

	public GenomicVariantImpl processClusterToVariant(List<GenomicVariant> clusterSigns, String sequenceName){
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
		short variantScore = calculateVariantScore(clusterSigns);
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
		variant.setVariantQS(variantScore);
		if(variant.getFirst() == 1716434 && variant.getLast() == 1716435 && variant.length() == 303) {
			System.out.println("$WEIRD variant processed from cluster: " + variant.getSequenceName() + " begin: " + variant.getFirst() + " end: " +
					variant.getLast() + " length: " + variant.length() + " type: " + 
					GenomicVariantImpl.getVariantTypeName(variant.getType()));
			for (GenomicVariant sign:clusterSigns) {
				System.out.println("$WEIRD Signature in cluster: " + sign.getSequenceName() + " begin: " + sign.getFirst() + " end: " +
						sign.getLast() + " length: " + sign.length() + " type: " + 
						GenomicVariantImpl.getVariantTypeName(sign.getType()));
			}
		}
		return variant;
	}
	
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
	}
	
	public short calculateSvimVariantScore(List<GenomicVariant> candidates) {
		double score = 0;
		SampleStatistics calcPos = new SampleStatistics();
		SampleStatistics calcSpan = new SampleStatistics();
		for(GenomicRegion sign:candidates) {
			calcPos.update(sign.getFirst());
			calcSpan.update(sign.length());
		}
		int n = candidates.size();
		double numSign = Math.min(80, n);
		double spanMean = calcSpan.getMean();
		double normPos = (numSign/8)*(1-Math.min(1, 
						(Math.sqrt(calcPos.getVariance()))/spanMean));
		double normSpan = (numSign/8)*(1-Math.min(1, 
						(Math.sqrt(calcSpan.getVariance()))/spanMean));
		score = normSpan + normPos + numSign;
		return (short) score;
	}
	

	public short calculateVariantScore(List<GenomicVariant> candidates) {
		double score = 0;
		SampleStatistics calcPos = new SampleStatistics();
		SampleStatistics calcSpan = new SampleStatistics();
		for(GenomicRegion sign:candidates) {
			calcPos.update(sign.getFirst());
			calcSpan.update(sign.length());
		}
		int n = candidates.size();
		double numSign = Math.min(80, n);
		double spanMean = calcSpan.getMean();
		double normPos = (numSign/8)*(1-Math.min(1, 
						(Math.sqrt(calcPos.getVariance()))/spanMean));
		double normSpan = (numSign/8)*(1-Math.min(1, 
						(Math.sqrt(calcSpan.getVariance()))/spanMean));
		score = normSpan + normPos + numSign;
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
	private List<CalledGenomicVariant> makeGenotypeCalls(List<GenomicVariant> variantsList) {
		List<CalledGenomicVariant> genotypeCalls = new ArrayList<>();
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
				System.out.println("call=" + candidateIndelCall.getSequenceName() + " " + candidateIndelCall.getFirst() 
						+ " " + candidateIndelCall.getLast() + " "  + candidateIndelCall.length());
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
	
	private List<VCFRecord> buildRecords(List<GenomicVariant> variants, List<CalledGenomicVariant> genotypeCalls, VCFFileHeader header){
		List<VCFRecord> records = new ArrayList<>();
		for(int i = 0; i < variants.size(); i++) {
			GenomicVariant variant = variants.get(i);
			CalledGenomicVariant calledVariant = genotypeCalls.get(i);
			List<CalledGenomicVariant> calls = new ArrayList<>();
			calls.add(calledVariant);
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
	
	public static Distribution calculateDistribution(List<Integer> values) {
		Collections.sort(values);
		Distribution dist = new Distribution((double) values.get(0), (double) values.get(values.size()-1), 1);
		for(double value : values) {
			dist.processDatapoint(value);
		}
		return dist;
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
		List<CalledGenomicVariant> genotypeCalls = makeGenotypeCalls(variantsList);
		VCFFileHeader header = createVCFHeader(sampleId);
		//filterVariantsWithoutGenotype(variantsList, genotypeCalls);
		List<VCFRecord> records = buildRecords(variantsList, genotypeCalls, header);
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
	
