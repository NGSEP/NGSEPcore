/*******************************************************************************
 * NGSEP - Next Generation Sequencing Experience Platform
 * Copyright 2016 Jorge Duitama
 *
 * This file is part of NGSEP.
 *
 *     NGSEP is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     NGSEP is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with NGSEP.  If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/
package ngsep.vcf;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.logging.Logger;

import ngsep.genome.ReferenceGenome;
import ngsep.main.CommandsDescriptor;
import ngsep.main.OptionValuesDecoder;
import ngsep.main.ProgressNotifier;
import ngsep.transcriptome.Transcriptome;
import ngsep.transcriptome.VariantAnnotationParameters;
import ngsep.transcriptome.VariantFunctionalAnnotation;
import ngsep.transcriptome.VariantFunctionalAnnotationType;
import ngsep.transcriptome.io.GFF3TranscriptomeHandler;
import ngsep.variants.GenomicVariant;
import ngsep.variants.GenomicVariantAnnotation;

/**
 * 
 * @author Jorge Duitama
 *
 */
public class VCFFunctionalAnnotator {
	
	// Constants for default values
	public static final int DEF_UPSTREAM=VariantAnnotationParameters.DEF_UPSTREAM;
	public static final int DEF_DOWNSTREAM=VariantAnnotationParameters.DEF_DOWNSTREAM;
	public static final int DEF_SPLICE_DONOR=VariantAnnotationParameters.DEF_SPLICE_DONOR;
	public static final int DEF_SPLICE_ACCEPTOR=VariantAnnotationParameters.DEF_SPLICE_ACCEPTOR;
	public static final int DEF_SPLICE_REGION_INTRON=VariantAnnotationParameters.DEF_SPLICE_REGION_INTRON;
	public static final int DEF_SPLICE_REGION_EXON=VariantAnnotationParameters.DEF_SPLICE_REGION_EXON;
	
	// Logging and progress
	private Logger log = Logger.getLogger(VCFFunctionalAnnotator.class.getName());
	private ProgressNotifier progressNotifier=null;
	
	// Parameters
	private String inputFile = null;
	private ReferenceGenome genome;
	private String transcriptomeFile = null;
	private Transcriptome transcriptome;
	private String outputFile = null;
	private VariantAnnotationParameters parameters = new VariantAnnotationParameters();
	
	// Get and set methods
	public Logger getLog() {
		return log;
	}
	public void setLog(Logger log) {
		if (log == null) throw new NullPointerException("Log can not be null");
		this.log = log;
	}

	public ProgressNotifier getProgressNotifier() {
		return progressNotifier;
	}
	public void setProgressNotifier(ProgressNotifier progressNotifier) {
		this.progressNotifier = progressNotifier;
	}
		
	public String getInputFile() {
		return inputFile;
	}
	public void setInputFile(String inputFile) {
		this.inputFile = inputFile;
	}
	
	public ReferenceGenome getGenome() {
		return genome;
	}
	public void setGenome(ReferenceGenome genome) {
		this.genome = genome;
	}
	public void setGenome(String genomeFile) throws IOException {
		setGenome(OptionValuesDecoder.loadGenome(genomeFile,log));
	}
	
	public String getTranscriptomeFile() {
		return transcriptomeFile;
	}
	public void setTranscriptomeFile(String transcriptomeFile) {
		this.transcriptomeFile = transcriptomeFile;
	}
	
	public String getOutputFile() {
		return outputFile;
	}
	public void setOutputFile(String outputFile) {
		this.outputFile = outputFile;
	}
	
	public int getOffsetUpstream() {
		return parameters.getOffsetUpstream();
	}
	public void setOffsetUpstream(int offsetUpstream) {
		parameters.setOffsetUpstream(offsetUpstream);
	}
	public void setOffsetUpstream(String value) {
		setOffsetUpstream((int) OptionValuesDecoder.decode(value, Integer.class));
	}
	
	public int getOffsetDownstream() {
		return parameters.getOffsetDownstream();
	}
	public void setOffsetDownstream(int offsetDownstream) {
		parameters.setOffsetDownstream(offsetDownstream);
	}
	public void setOffsetDownstream(String value) {
		setOffsetDownstream((int) OptionValuesDecoder.decode(value, Integer.class));
	}
	
	/**
	 * @return
	 * @see ngsep.transcriptome.VariantAnnotationParameters#getSpliceDonorOffset()
	 */
	public int getSpliceDonorOffset() {
		return parameters.getSpliceDonorOffset();
	}
	/**
	 * @param spliceDonorOffset
	 * @see ngsep.transcriptome.VariantAnnotationParameters#setSpliceDonorOffset(int)
	 */
	public void setSpliceDonorOffset(int spliceDonorOffset) {
		parameters.setSpliceDonorOffset(spliceDonorOffset);
	}
	public void setSpliceDonorOffset(String value) {
		setSpliceDonorOffset((int) OptionValuesDecoder.decode(value, Integer.class));
	}
	/**
	 * @return
	 * @see ngsep.transcriptome.VariantAnnotationParameters#getSpliceAcceptorOffset()
	 */
	public int getSpliceAcceptorOffset() {
		return parameters.getSpliceAcceptorOffset();
	}
	/**
	 * @param spliceAcceptorOffset
	 * @see ngsep.transcriptome.VariantAnnotationParameters#setSpliceAcceptorOffset(int)
	 */
	public void setSpliceAcceptorOffset(int spliceAcceptorOffset) {
		parameters.setSpliceAcceptorOffset(spliceAcceptorOffset);
	}
	public void setSpliceAcceptorOffset(String value) {
		this.setSpliceAcceptorOffset((int) OptionValuesDecoder.decode(value, Integer.class));
	}
	
	/**
	 * @return
	 * @see ngsep.transcriptome.VariantAnnotationParameters#getSpliceRegionIntronOffset()
	 */
	public int getSpliceRegionIntronOffset() {
		return parameters.getSpliceRegionIntronOffset();
	}
	/**
	 * @param spliceRegionIntronOffset
	 * @see ngsep.transcriptome.VariantAnnotationParameters#setSpliceRegionIntronOffset(int)
	 */
	public void setSpliceRegionIntronOffset(int spliceRegionIntronOffset) {
		parameters.setSpliceRegionIntronOffset(spliceRegionIntronOffset);
	}
	public void setSpliceRegionIntronOffset(String value) {
		this.setSpliceRegionIntronOffset((int) OptionValuesDecoder.decode(value, Integer.class));
	}
	/**
	 * @return
	 * @see ngsep.transcriptome.VariantAnnotationParameters#getSpliceRegionExonOffset()
	 */
	public int getSpliceRegionExonOffset() {
		return parameters.getSpliceRegionExonOffset();
	}
	/**
	 * @param spliceRegionExonOffset
	 * @see ngsep.transcriptome.VariantAnnotationParameters#setSpliceRegionExonOffset(int)
	 */
	public void setSpliceRegionExonOffset(int spliceRegionExonOffset) {
		parameters.setSpliceRegionExonOffset(spliceRegionExonOffset);
	}
	public void setSpliceRegionExonOffset(String value) {
		this.setSpliceRegionExonOffset((int) OptionValuesDecoder.decode(value, Integer.class));
	}
	
	public static void main(String[] args) throws Exception {
		VCFFunctionalAnnotator instance = new VCFFunctionalAnnotator();
		CommandsDescriptor.getInstance().loadOptions(instance, args);
		instance.run();
	}
	
	public void run() throws IOException {
		logParameters();
		if (genome == null) throw new IOException("The file with the reference genome is a required parameter");
		loadMap(transcriptomeFile, genome);
		if(inputFile==null) {
			
			if(outputFile == null) annotate(System.in, System.out);
			else annotate(System.in, outputFile);
		} else {
			if(outputFile == null) annotate(inputFile,System.out);
			else {
				annotate(inputFile,outputFile);
			}
		}
		log.info("Process finished");
	}

	private void logParameters() {
		ByteArrayOutputStream os = new ByteArrayOutputStream();
		PrintStream out = new PrintStream(os);
		if(inputFile != null) out.println("Input file: "+inputFile);
		else out.println("Read variants from standard input");
		out.println("GFF transcriptome file: "+transcriptomeFile);
		if (genome!=null) out.println("Loaded reference genome from: "+genome.getFilename());
		if(outputFile != null) out.println("Output file: "+outputFile);
		else out.println("Annotated VCF will be written to standard output");
		
		out.println("Upstream offset: "+getOffsetUpstream());
		out.println("Downstream offset: "+getOffsetDownstream());
		out.println("Splice donor offset: "+getSpliceDonorOffset());
		out.println("Splice acceptor offset: "+getSpliceAcceptorOffset());
		out.println("Splice region intron offset: "+getSpliceRegionIntronOffset());
		out.println("Splice region exon offset: "+getSpliceRegionExonOffset());
		
		log.info(""+os.toString());
	}
	public void loadMap(String transcriptomeMap,String transcriptomeCDNA) throws IOException {
		GFF3TranscriptomeHandler handler = new GFF3TranscriptomeHandler();
		handler.setLog(log);
		transcriptome = handler.loadMap(transcriptomeMap);
		handler.loadSequences(transcriptome, transcriptomeCDNA);
	}
	public void loadMap(String transcriptomeMap, ReferenceGenome genome) throws IOException {
		GFF3TranscriptomeHandler handler = new GFF3TranscriptomeHandler(genome.getSequencesMetadata());
		handler.setLog(log);
		transcriptome = handler.loadMap(transcriptomeMap);
		transcriptome.fillSequenceTranscripts(genome, log);
	}
	public void annotate(String variantsFile,String outputFile) throws IOException {
		try (VCFFileReader in = new VCFFileReader(variantsFile);
			 PrintStream out = new PrintStream(outputFile)){
			annotate(in,out);
		}
	}
	public void annotate(InputStream input,String outputFile) throws IOException {
		try (VCFFileReader in = new VCFFileReader(input);
			 PrintStream out = new PrintStream(outputFile)){
			annotate(in,out);
		}
	}
	public void annotate(InputStream input, PrintStream out) throws IOException {
		try (VCFFileReader in = new VCFFileReader(input)) {
			annotate(in,out);
		}
	}
	public void annotate(String variantsFile,PrintStream out) throws IOException {
		try (VCFFileReader in = new VCFFileReader(variantsFile)){
			annotate(in,out);
		}
	}
	public void annotate(VCFFileReader in,PrintStream out) throws IOException {
		VCFFileWriter writer = new VCFFileWriter(); 
		in.setLog(log);
		writer.printHeader(in.getHeader(),out);
		Iterator<VCFRecord> it = in.iterator();
		int n=0;
		while (it.hasNext()) {
			VCFRecord record = it.next();
			if(record.getVariant().getAlleles().length>=2) annotate(record);
			writer.printVCFRecord(record, out);
			n++;
			if (progressNotifier!=null && n%1000==0) {
				int progress = n/1000;
				if (!progressNotifier.keepRunning(progress)) {
					out.flush();
					return;
				}
			}
		}
		out.flush();
	}

	public void annotate(VCFRecord record) {
		GenomicVariant v = record.getVariant();
		List<VariantFunctionalAnnotation> annotations = transcriptome.calculateAnnotations(v, parameters);
		if(annotations.size()>0) {
			Collections.sort(annotations,VariantFunctionalAnnotationType.getPriorityComparator());
			VariantFunctionalAnnotation ann = annotations.get(0);
			record.addAnnotation(new GenomicVariantAnnotation(v, GenomicVariantAnnotation.ATTRIBUTE_TRANSCRIPT_ANNOTATION, ann.getTypeName()));
			if(ann.getTranscript()!=null) {
				record.addAnnotation(new GenomicVariantAnnotation(v, GenomicVariantAnnotation.ATTRIBUTE_TRANSCRIPT_ID, ann.getTranscript().getId()));
				record.addAnnotation(new GenomicVariantAnnotation(v, GenomicVariantAnnotation.ATTRIBUTE_GENE_NAME, ann.getTranscript().getGeneName()));
				if(ann.getCodonNumber()>0) {
					record.addAnnotation(new GenomicVariantAnnotation(v, GenomicVariantAnnotation.ATTRIBUTE_TRANSCRIPT_CODON, ""+ann.getCodonNumber()+"."+ann.getCodonPosition()));
				}
				if(ann.getAminoacidChange()!=null) {
					record.addAnnotation(new GenomicVariantAnnotation(v, GenomicVariantAnnotation.ATTRIBUTE_TRANSCRIPT_AMINOACID_CHANGE, ann.getAminoacidChange()));
				}
			}
		}
	}
	
}
