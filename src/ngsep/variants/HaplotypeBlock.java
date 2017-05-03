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
	 * Returns Hamming distance between two fragments
	 * <b> pre: </b> The matrix of fragments has been initialized.
	 * @param row1. Row1 is bigger than row2
	 * @param row2.
	 * @exception row1 must be bigger than row2.
	 * @return hamming distance between two fragments.
	 */
	public int getHammingDistance(int row1, int row2) 
	{
		int distance = 0;
		HaplotypeFragment fragment1 = matrix.get(row1);
		HaplotypeFragment fragment2 = matrix.get(row2);
		byte[] callsF1 = fragment1.getCalls();
		byte[] callsF2 = fragment2.getCalls();
		for(int i = fragment1.getFirstColum(); i < callsF1.length ; i++ )
		{
			byte callF1 = fragment1.getCall(i);
			for(int j = fragment2.getFirstColum(); j < callsF2.length ; j++)
			{
				byte callF2 = fragment2.getCall(j);
				if(callF1 != callF2)
				{
					distance ++;
				}
			}
		}
		return distance;
	}
	/**
	 * 
	 * @param row1
	 * @param row2
	 * @return
	 */
	public boolean overlap(int row1, int row2)
	{
		boolean overlap = false;		
		HaplotypeFragment fragment2 = matrix.get(row2);
		HaplotypeFragment fragment1 = matrix.get(row1);
		int initPosF1 = fragment1.getFirstColum();
		int initPosF2 = fragment2.getFirstColum();
		int lengthF1 = fragment1.getCalls().length;
		int lengthF2 = fragment2.getCalls().length;
		int lastPosF2 = initPosF2 + lengthF2 -1;
		for (int i = initPosF1; i < lengthF1 ; i++ )
		{
			byte callF1 = fragment1.getCall(i);
			for(int j = lastPosF2 ; j > initPosF2; j--)
			{
				byte callF2 = fragment2.getCall(j);
			}
		}
		
		return overlap;
	}
	


}
