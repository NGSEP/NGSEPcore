package ngsep.haplotyping;

import java.util.ArrayList;
import java.util.List;

import ngsep.alignments.ReadAlignment;
import ngsep.genome.GenomicRegionPositionComparator;
import ngsep.math.NumberArrays;
import ngsep.variants.CalledGenomicVariant;
import ngsep.variants.GenomicVariant;

public class HaplotypeBlockBuilder {
	/**
	 * Builds a haplotyppe block for molecular phasing with the given data
	 * @param hetVars List of heterozygous variants within a chromosome
	 * @param alignments List of alignments to the same chromosome where the variants are located 
	 * @return HaplotypeBlock Matrix for phasing
	 */
	public static HaplotypeBlock buildBlock (List<? extends GenomicVariant> hetVars, List<ReadAlignment> alignments) {
		HaplotypeBlock block = new HaplotypeBlock(hetVars);
		int i=0;
		for(ReadAlignment aln:alignments) {
			//Advance i
			GenomicVariant firstHetVar = null;
			while(i<hetVars.size()) {
				firstHetVar = hetVars.get(i);
				if(GenomicRegionPositionComparator.getInstance().compare(firstHetVar, aln)>=0) {
					break;
				}
				i++;
			}
			if(i==hetVars.size()) {
				break;
			}
			//Extract relevant calls from alignment
			int lastAln = aln.getLast();
			List<Byte> calls = new ArrayList<>(50);
			for(int j=i;j<hetVars.size();j++) {
				GenomicVariant var = hetVars.get(j);
				if(var.getFirst()>lastAln) {
					break;
				}
				String [] alleles = var.getAlleles();
				String call = aln.getAlleleCall(var.getFirst(), var.getLast()).toString();
				if(alleles[0].equals(call)) {
					calls.add(CalledGenomicVariant.ALLELE_REFERENCE);
				} else if(alleles[1].equals(call)) {
					calls.add(CalledGenomicVariant.ALLELE_ALTERNATIVE);
				} else if (calls.size()==0) {
					i=j+1;
				} else {
					calls.add(CalledGenomicVariant.ALLELE_UNDECIDED);
				}
			}
			//Trim last undecided calls
			for(int j=calls.size()-1;j>=0;j--) {
				Byte call = calls.get(j);
				if(call!=CalledGenomicVariant.ALLELE_UNDECIDED) {
					break;
				}
				calls.remove(call);
			}
			if(calls.size()>1) {
				block.addFragment (i,NumberArrays.toByteArray(calls));
			}
		}
		return block;
	}
}
