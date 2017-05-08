package ngsep.haplotyping;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

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
	public void process(String vcfFilename, String bamFilename, PrintStream out) {
		List<VCFRecord> sequenceRecords = new ArrayList<>();
		//TODO: Hace un ciclo que recorre el VCF y va agregando los VCFRecord a la lista sequenceRecords
		//Cada vez que cambias de cromosoma, procesas la lista y la desocupas
		//Te puedes guiar por GenotypeImputer
		
	}
	private void processChormosomeVariants (List<VCFRecord> chromosomeVariants, PrintStream out) {
		//TODO: Por ahora solamente imprime el total de variantes y cuantas son heterocigotas
	}

}
