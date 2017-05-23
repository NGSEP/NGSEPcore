package ngsep.haplotyping;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Iterator;

import ngsep.vcf.VCFFileHeader;
import ngsep.vcf.VCFFileReader;
import ngsep.vcf.VCFRecord;

public class SingleIndividualHaplotyper {

	private String algorithmClassName = "ngsep.haplotyping.RefhapSIHAlgorithm";
	private SIHAlgorithm algClass;
	public static void main(String[] args) throws Exception {
		SingleIndividualHaplotyper instance = new SingleIndividualHaplotyper();
		String vcfFilename = args[0];
		String bamFilename = args[1];
		PrintStream out = new PrintStream(args[2]);
		instance.process (vcfFilename,bamFilename,out);

	}
	/**
	 * 
	 * @param vcfFilename Input VCF
	 * @param bamFilename Input BAM file
	 * @param out Output VCF
	 */
	public void process(String vcfFilename, String bamFilename, PrintStream out) throws Exception{
		List<VCFRecord> sequenceRecords = new ArrayList<>();
		//TODO: Hace un ciclo que recorre el VCF y va agregando los VCFRecord a la lista sequenceRecords
		//Cada vez que cambias de cromosoma, procesas la lista y la desocupas
		//Te puedes guiar por GenotypeImputer
		VCFFileReader inputVCF = null;
		FileInputStream inputBAM = null;
		PrintStream outVCF = null;
		try
		{
			inputVCF = new VCFFileReader(vcfFilename);
			inputBAM = new FileInputStream(bamFilename);
			out = new PrintStream(outVCF);
			Iterator<VCFRecord> iter = inputVCF.iterator();
			VCFFileHeader header = inputVCF.getHeader();
			while(iter.hasNext())
			{
				VCFRecord record = iter.next();
			
			}

		}
		finally{
			if(inputVCF!=null) inputVCF.close();
			if(inputBAM!=null) inputBAM.close();
			if(out!=null)out.close();
		}


	}
	private void processChormosomeVariants (List<VCFRecord> chromosomeVariants, PrintStream out) {
		//TODO: Por ahora solamente imprime el total de variantes y cuantas son heterocigotas
		PrintStream outVariants = null;
		int variants = 0;
		try
		{
			out = new PrintStream(outVariants);
			variants = chromosomeVariants.size();
			for(int i = 0; i < variants; i++)
			{
				VCFRecord current = chromosomeVariants.get(i);
						
			}
		}
		finally
		{
			if(outVariants!=null) outVariants.close();
		}

	}

}
