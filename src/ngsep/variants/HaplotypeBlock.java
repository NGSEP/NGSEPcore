package ngsep.variants;

import java.util.ArrayList;
import java.util.List;

public class HaplotypeBlock {

	/**
	 * Represents the matrix of fragments and variants,
	 */
	private List<HaplotypeFragment> matrix;

	/**
	 * Represents the list of variants
	 */
	private List<GenomicVariant> variants;
	
	/**
	 * Represents a haplotype.
	 */
	private byte haplotype[];
	
	/**
	 * 
	 * @param fragments
	 * @param variants
	 * @throws Exception
	 */

	public HaplotypeBlock(List<GenomicVariant> variants) throws Exception 
	{
		this.variants = variants;
		matrix = new ArrayList <HaplotypeFragment>();
		haplotype = null;

	}
	
	public byte getAllele(int i, int j)
	{
		HaplotypeFragment row = matrix.get(i);
		byte allele = row.getCall(j);
		return allele;
	}
	
	public byte [] getHaplotype()
	{
		return haplotype;
	}
	
	public GenomicVariant getVariant(int i)
	{
		GenomicVariant temp = variants.get(i);
		return  temp;
	}
	
	/**
	 * 	
	 * @param row1. Row1 is bigger than row2
	 * @param row2
	 * @return
	 */
	public int getHammingDistance(int row1, int row2)
	{
		return 0;
	}
	/**
	 * 
	 * @param row1
	 * @param row2
	 * @return
	 */
	public boolean overlap(int row1, int row2)
	{
		return false;
	}
	


}
