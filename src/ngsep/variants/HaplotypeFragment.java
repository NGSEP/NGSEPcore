package ngsep.variants;

public class HaplotypeFragment {
	
	private int firstColumn;
	private byte calls [];
	
	
	public HaplotypeFragment(int firstColum, byte[] calls) {
		
		this.firstColumn = firstColum;
		this.calls = calls;
	}

	
	
	public int getFirstColumn() {
		return firstColumn;
	}

	public byte[] getCalls() {
		return calls;
	}
	
	public byte getCall(int column)
	{
		byte call = CalledGenomicVariant.ALLELE_UNDECIDED;
		int posIni = getFirstColumn();
		int length = calls.length;
		int posLas = posIni + length - 1;
		if (column >= posIni && column <= posLas)
		{
			int relativePosition = column - posIni;
			call = calls[relativePosition];
		}
		return call;	
		
	}
	/**
	 * Get last column
	 */
	public int getLastColumn()
	{
		return 0;
	}
	
	

}
