package ngsep.variants;

public class HaplotypeFragment {
	
	private int firstColum;
	private byte calls [];
	
	
	public HaplotypeFragment(int firstColum, byte[] calls) {
		
		this.firstColum = firstColum;
		this.calls = calls;
	}

	
	
	public int getFirstColum() {
		return firstColum;
	}

	public byte[] getCalls() {
		return calls;
	}
	
	public byte getCall(int column)
	{
		byte call = CalledGenomicVariant.ALLELE_UNDECIDED;
		int posIni = getFirstColum();
		int length = calls.length;
		int posLas = posIni + length - 1;
		if (column >= posIni && column <= posLas)
		{
			int relativePosition = column - posIni;
			call = calls[relativePosition];
		}
		return call;	
		
	}
	
	

}
