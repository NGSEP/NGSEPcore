package ngsep.sequences;

public class DirectShortKmerCodesHashFunction implements ShortKmerCodesHashFunction {

	@Override
	public int getHash(long code) {
		if(code <Integer.MAX_VALUE) {
			return (int) code;
		}
		throw new IllegalArgumentException("Code "+code+" can not be used for direct hashing");
	}

}
