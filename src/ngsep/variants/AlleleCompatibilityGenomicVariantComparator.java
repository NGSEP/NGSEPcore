package ngsep.variants;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import ngsep.genome.GenomicRegionComparator;
import ngsep.genome.ReferenceGenome;
import ngsep.sequences.QualifiedSequenceList;

public class AlleleCompatibilityGenomicVariantComparator implements Comparator<GenomicVariant> {

	private ReferenceGenome genome;
	private QualifiedSequenceList seqNames;
	private GenomicRegionComparator internalComparator;
	
	public AlleleCompatibilityGenomicVariantComparator(ReferenceGenome genome) {
		this.genome = genome;
		seqNames = genome.getSequencesMetadata();
		internalComparator = new GenomicRegionComparator(seqNames);
	}

	
	@Override
	public int compare(GenomicVariant v1, GenomicVariant v2) {
		int cmp1 = internalComparator.compare(v1, v2);
		if(cmp1>2 || cmp1<-2) return cmp1;
		if(v1 instanceof SNV && v2 instanceof SNV) return cmp1;
		//if(v2.getFirst()==753845) System.out.println("Variant 1: "+v1.getSequenceName()+":"+v1.getFirst()+"-"+v1.getLast()+"Comparator result: "+cmp1);
		if(v2.getFirst() - v1.getLast() > 2 || v1.getFirst() - v2.getLast() > 2) return cmp1;
		int firstRegion = Math.min(v1.getFirst(), v2.getFirst());
		int lastRegion = Math.max(v1.getLast(), v2.getLast());
		//if(v2.getFirst()==753845) System.out.println("Getting alleles between: "+firstRegion+" and "+lastRegion);
		List<String> alleleStrs1 = new ArrayList<>(buildAlleleStrings(v1,firstRegion,lastRegion));
		List<String> alleleStrs2 = new ArrayList<>(buildAlleleStrings(v2,firstRegion,lastRegion));
		//if(v2.getFirst()==753845) System.out.println("Alleles 1: "+alleleStrs1.toString()+" Alleles 2: "+alleleStrs2.toString());
		if(alleleStrs1.size()!=alleleStrs2.size()) return cmp1!=0?cmp1: alleleStrs1.size()-alleleStrs2.size();
		for(int i=0;i<alleleStrs1.size();i++) {
			String a1 = alleleStrs1.get(i);
			String a2 = alleleStrs2.get(i);
			if(!a1.equals(a2))  return cmp1!=0?cmp1:a1.compareTo(a2);
		}
		return 0;
	}
	private Set<String> buildAlleleStrings(GenomicVariant v, int firstRegion, int lastRegion) {
		String seqName = v.getSequenceName();
		CharSequence leftC = genome.getReference(seqName, Math.max(1, firstRegion-3), Math.max(1, v.getFirst()-1));
		int length = seqNames.get(seqName).getLength();
		CharSequence rightC = genome.getReference(seqName, Math.min(length, v.getLast()+1), Math.min(length, lastRegion+3));
		String left = leftC!=null?leftC.toString().toUpperCase():"";
		String right = rightC!=null?rightC.toString().toUpperCase():"";
		Set<String> allelesSet = new TreeSet<>();
		String [] alleles = v.getAlleles();
		for(int i=0;i<alleles.length;i++) {
			allelesSet.add((left+alleles[i]+right).toUpperCase());
		}
		return allelesSet;
	}
		

}
