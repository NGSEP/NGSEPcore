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

	public HaplotypeBlock(List<GenomicVariant> variants) 
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
	 * minimo de los comienzos con maximo del final.
	 * @return hamming distance between two fragments.
	 */
	public int getHammingDistance(int row1, int row2) 
	{
		int distance = 0;
		HaplotypeFragment fragment1 = matrix.get(row1);
		HaplotypeFragment fragment2 = matrix.get(row2);
		byte[] callsF1 = fragment1.getCalls();
		byte[] callsF2 = fragment2.getCalls();
		for(int i = fragment1.getFirstColumn(); i < callsF1.length ; i++ )
		{
			byte callF1 = fragment1.getCall(i);
			for(int j = fragment2.getFirstColumn(); j < callsF2.length ; j++)
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
		int initPosF1 = fragment1.getFirstColumn();
		int initPosF2 = fragment2.getFirstColumn();
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
	/**
	 * Dada una fila de la matriz cual es la ultima posicion que no tiene dato perdido.
	 * reciben la fila que van a consultar
	 */
	public int getFirstColumn(int row)
	{
		return 0;
	}
	
	/**
	 * Dada la fila de la matriz decir cual es la primera posicion que no tiene dato perdido.
	 * reciben la fila que van a consultar
	 */
	public int getLastColumn(int row)
	{
		return 0;
	}
	/**
	 * Get number of fragments in the block. getNumFragments()
	 */
	public int getNumFragments()
	{
		return 0;
	}
	/**
 	* GetHamming distance2 (hamming2). if it is the same then add -1, if they are different +1, if it is compared with a null allele then adds nothing
 	* calcular el score de los fragmentos comparados. cambiar por este en getScore
 	*/
	public int getHamming2(int row1 , int row2)
	{
		return 0;
	}
	/**
	 * getNumVariants
	 */
	public int getNumVariants()
	{
		return 0;
	}
	public void setHaplotype(byte [] haplotype)
	{
		
	}
}
