package ngsep.haplotyping;

import java.io.FileInputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Iterator;

import ngsep.variants.CalledGenomicVariant;
import ngsep.vcf.VCFFileHeader;
import ngsep.vcf.VCFFileReader;
import ngsep.vcf.VCFFileWriter;
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
		List<VCFRecord> chromosomeVariants = new ArrayList<>();
		//TODO: Hace un ciclo que recorre el VCF y va agregando los VCFRecord a la lista sequenceRecords
		//Cada vez que cambias de cromosoma, procesas la lista y la desocupas
		//Te puedes guiar por GenotypeImputer
		VCFFileReader inputVCF = null;
		VCFFileWriter vcfWriter = new VCFFileWriter();
		try
		{
			inputVCF = new VCFFileReader(vcfFilename);
			VCFFileHeader header = inputVCF.getHeader();
			vcfWriter = new VCFFileWriter();
			vcfWriter.printHeader(header, out);
			String lastSeqName = null;
			Iterator<VCFRecord> iter = inputVCF.iterator();
			
			while(iter.hasNext())
			{
				VCFRecord record = iter.next();
				if(record.getSequenceName().equals(lastSeqName)) {
					if(chromosomeVariants.size()>0) {
						processChromosomeVariants(lastSeqName, chromosomeVariants, bamFilename, vcfWriter, out);
					}
					chromosomeVariants.clear();
					lastSeqName = record.getSequenceName();
				}
				chromosomeVariants.add(record);
			}
			if(chromosomeVariants.size()>0) {
				processChromosomeVariants(lastSeqName, chromosomeVariants, bamFilename, vcfWriter, out);
			}
		}
		finally{
			if(inputVCF!=null) inputVCF.close();
			if(out!=null)out.close();
		}
	}
	private void processChromosomeVariants (String seqName, List<VCFRecord> chromosomeVariants, String bamFilename, VCFFileWriter vcfWriter, PrintStream out) {
		List<CalledGenomicVariant> hetVars = extractHeterozygousVariants(chromosomeVariants);
		System.out.println("Heterozygous variants: "+hetVars.size());
		//Read bam file. Ignore alignments to chromosomes different than seqName
		//Crear lista de ReadAlignment con los alineamientos del cromosoma y cuando acabes crear un solo 
		//HaplotypeBlock llamando HaplotypeBlockBuilder
		//Y luego pides que 
		
		

	}
	private List<CalledGenomicVariant> extractHeterozygousVariants(List<VCFRecord> chromosomeVariants) {
		List<CalledGenomicVariant> hetVars = new ArrayList<>();
		for(VCFRecord record:chromosomeVariants)
		{
			CalledGenomicVariant call = record.getCalls().get(0);
			if(call.isHeterozygous());
			hetVars.add(call);
		}
		return hetVars;
	}

}
