package ngsep.variants;

import java.util.List;

public class HaplotypeBlock {

	/**
	 * Represents the matrix of fragments and variants,
	 */
	private byte fragments [][];

	/**
	 * Represents the list of variants
	 */
	private List<GenomicVariant> variants;
	
	/**
	 * Represents a haplotype.
	 */
	private String haplotype;
	
	/**
	 * 
	 * @param fragments
	 * @param variants
	 * @throws Exception
	 */

	public HaplotypeBlock(byte fragments [][],  List<GenomicVariant> variants) throws Exception 
	{
		if(fragments[0].length==variants.size())
		{
			this.fragments = fragments;
			this.variants = variants;
		} else {
			throw new Exception ("The number of columns in the fragments matrix must be the same of rows in the list of variants");
		}

	}
	
	public String getHaplotype()
	{
		return haplotype;
	}

}
