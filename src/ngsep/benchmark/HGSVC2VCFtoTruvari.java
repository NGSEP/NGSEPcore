package ngsep.benchmark;

import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import ngsep.discovery.LongReadStructuralVariantDetector;
import ngsep.genome.ReferenceGenome;
import ngsep.sequences.QualifiedSequence;
import ngsep.variants.CalledGenomicVariant;
import ngsep.variants.CalledGenomicVariantImpl;
import ngsep.variants.GenomicVariant;
import ngsep.variants.GenomicVariantAnnotation;
import ngsep.variants.GenomicVariantImpl;
import ngsep.variants.Sample;
import ngsep.vcf.VCFFileHeader;
import ngsep.vcf.VCFFileReader;
import ngsep.vcf.VCFFileWriter;
import ngsep.vcf.VCFHeaderLine;
import ngsep.vcf.VCFRecord;

public class HGSVC2VCFtoTruvari {
	public static final String CONTIG_HEADER_TYPE = "contig";
    public static final String CONTIG_LENGTH_HEADER_ATTRIBUTE = "length";
    public static final String CONTIG_HEADER_PREFIX = "##contig=<ID=";
    public static final String CONTIG_HEADER_POST_ID = ",length=";
    public static final String CONTIG_HEADER_SUFFIX = ">";
    public static final String FILTER_SETTING_MISSING = LongReadStructuralVariantDetector.FILTER_SETTING_MISSING;

    public static void main(String[] args) throws IOException {
        processVCF(args[0], args[1], args[2]);
    }
    
    public static void processVCF(String file, String refGenFile, String outFile) throws IOException {
        try(VCFFileReader reader = new VCFFileReader(file);){
        	System.out.println("Loading reference genome");
            ReferenceGenome refGenome = new ReferenceGenome(refGenFile);
            System.out.println("Loaded reference genome");
            List<CalledGenomicVariant> variants = new ArrayList<>();
            VCFFileHeader oldHeader = reader.getHeader();
            System.out.println("Loaded old header");
            
	        Iterator<VCFRecord> it = reader.iterator();
            
            System.out.println("Began processing of vcf records");
            while(it.hasNext()){
            	VCFRecord oldRecord = it.next();
                //System.out.println(old_record.getVariant().getLast());
            	GenomicVariantImpl impl = (GenomicVariantImpl) oldRecord.getVariant();
            	//System.out.println(impl.getSequenceName() + " " + impl.getFirst() + " " + impl.getLast() + " " + impl.length() + " " + 
            	//	GenomicVariantImpl.getVariantTypeName(impl.getType()));
                if(impl.getType() == GenomicVariant.TYPE_LARGEINS) {
                    int last = impl.getFirst() + 1;
                    impl.setLast(last);
                }
                int genotype;
                //System.out.println(oldRecord.getCalls().size()!=1);
                CalledGenomicVariantImpl oldCall = (CalledGenomicVariantImpl) oldRecord.getCalls().get(0);
                if(oldCall.isHeterozygous()) genotype = 1;
                else if (oldCall.isHomozygous() && !oldCall.isHomozygousReference()) genotype = 2;
                else genotype = 0;
                CalledGenomicVariant newCall = new CalledGenomicVariantImpl(impl, genotype);
                newCall.setSampleId(oldHeader.getSampleIds().get(0));
                variants.add(newCall);
            }
            System.out.println("Finished processing of vcf records");
            VCFFileHeader newHeader = VCFFileHeader.makeDefaultEmptyHeader();
            for(QualifiedSequence qSeq : refGenome.getSequencesList()){
                String id = qSeq.getName();
                int len = qSeq.getLength();
                VCFHeaderLine newContigLine = makeHeaderContigLine(id, len);
                newHeader.addHeaderLine(newContigLine);
            }
            String sample = variants.get(0).getSampleId();
            Sample s = new Sample(sample);
            newHeader.addSample(s, true);
            newHeader.addMissingEntries();
            System.out.println("Finished production of new header");
            saveVCFResultsFile(variants, newHeader, sample, outFile);
        }
    }

    private static VCFHeaderLine makeHeaderContigLine(String contigId, int contigLength){
        //String line = CONTIG_HEADER_PREFIX + contigId + CONTIG_HEADER_POST_ID + contigLength + CONTIG_HEADER_SUFFIX;
        VCFHeaderLine answer = new VCFHeaderLine(CONTIG_HEADER_TYPE, contigId);
        answer.setAttribute(VCFHeaderLine.ATTRIBUTE_ID, contigId);
        answer.setAttribute(CONTIG_LENGTH_HEADER_ATTRIBUTE, contigLength + "");
        return answer;
    }

    private static List<VCFRecord> buildRecords(List<CalledGenomicVariant> genotypeCalls, VCFFileHeader header){
        List<VCFRecord> records = new ArrayList<>();
        for(int i = 0; i < genotypeCalls.size(); i++) {
            GenomicVariant variant = genotypeCalls.get(i);
            List<CalledGenomicVariant> calls = new ArrayList<>();
            CalledGenomicVariant call = (CalledGenomicVariantImpl) variant;
            calls.add(call);
            List<GenomicVariantAnnotation> infoFields = GenomicVariantAnnotation.annotateStructuralVariant(variant);
            List<String> filters = getVariantFilters(variant);
            VCFRecord record = new VCFRecord(variant, filters, infoFields, VCFRecord.DEF_FORMAT_ARRAY_MINIMAL,  calls, header);
            records.add(record);
        }
        return records;
    }

    private static List<String> getVariantFilters(GenomicVariant variants){
        List<String> filters = new ArrayList<>();
        filters.add(FILTER_SETTING_MISSING);
        return filters;
    }

    private static void printVCFFile(List<VCFRecord> records,VCFFileHeader header, String file) throws IOException {
        try(PrintStream pr = new PrintStream(file);){
            VCFFileWriter writer = new VCFFileWriter();
            writer.printHeader(header, pr);
            writer.printVCFRecords(records, pr);
        }
    }

    private static void saveVCFResultsFile(List<CalledGenomicVariant> genotypeCalls, VCFFileHeader header, String sampleID, String outFile) throws IOException {
        List<VCFRecord> records = buildRecords(genotypeCalls, header);
        printVCFFile(records, header, outFile);
    }

    
}
