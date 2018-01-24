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

import java.io.IOException;
import java.io.PrintStream;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.logging.Logger;

import ngsep.genome.ReferenceGenome;
import ngsep.main.CommandsDescriptor;
import ngsep.main.ProgressNotifier;
import ngsep.transcriptome.Transcript;
import ngsep.transcriptome.Transcriptome;
import ngsep.transcriptome.VariantAnnotationParameters;
import ngsep.transcriptome.VariantFunctionalAnnotation;
import ngsep.transcriptome.io.GFF3TranscriptomeHandler;
import ngsep.variants.GenomicVariant;
import ngsep.variants.GenomicVariantAnnotation;

public class VCFFunctionalAnnotator {
	private Logger log = Logger.getLogger(VCFFunctionalAnnotator.class.getName());
	private Transcriptome transcriptome;
	private VariantAnnotationParameters parameters = new VariantAnnotationParameters();
	
	private ProgressNotifier progressNotifier=null;
	
	public static void main(String[] args) throws Exception {
		VCFFunctionalAnnotator annotator = new VCFFunctionalAnnotator();
		int i= CommandsDescriptor.getInstance().loadOptions(annotator, args);
		String variantsFile = args[i++];
		String transcriptomeMap = args[i++];
		String sequenceFasta = args[i++]; 
		
		annotator.setLog(Logger.getLogger(VCFFunctionalAnnotator.class.getName()));
		annotator.loadMap(transcriptomeMap, new ReferenceGenome(sequenceFasta));
		annotator.annotate(variantsFile, System.out);
	}	
	public int getOffsetUpstream() {
		return parameters.getOffsetUpstream();
	}

	public void setOffsetUpstream(int offsetUpstream) {
		parameters.setOffsetUpstream(offsetUpstream);
	}
	
	public void setOffsetUpstream(Integer offsetUpstream) {
		setOffsetUpstream(offsetUpstream.intValue());
	}
	
	public int getOffsetDownstream() {
		return parameters.getOffsetDownstream();
	}

	public void setOffsetDownstream(int offsetDownstream) {
		parameters.setOffsetDownstream(offsetDownstream);
	}
	
	public void setOffsetDownstream(Integer offsetDownstream) {
		setOffsetDownstream(offsetDownstream.intValue());
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
	public void setSpliceDonorOffset(Integer spliceDonorOffset) {
		setSpliceDonorOffset(spliceDonorOffset.intValue());
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
	public void setSpliceAcceptorOffset(Integer spliceAcceptorOffset) {
		this.setSpliceAcceptorOffset(spliceAcceptorOffset.intValue());
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
	public void setSpliceRegionIntronOffset(Integer spliceRegionIntronOffset) {
		this.setSpliceRegionIntronOffset(spliceRegionIntronOffset.intValue());
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
	
	public void setSpliceRegionExonOffset(Integer spliceRegionExonOffset) {
		this.setSpliceRegionExonOffset(spliceRegionExonOffset.intValue());
	}
	
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

	public void loadMap(String transcriptomeMap,String transcriptomeCDNA) throws IOException {
		GFF3TranscriptomeHandler handler = new GFF3TranscriptomeHandler();
		handler.setLog(log);
		transcriptome = handler.loadMap(transcriptomeMap);
		handler.loadSequences(transcriptome, transcriptomeCDNA);
	}
	public void loadMap(String transcriptomeMap, ReferenceGenome genome) throws IOException {
		GFF3TranscriptomeHandler handler = new GFF3TranscriptomeHandler();
		handler.setLog(log);
		transcriptome = handler.loadMap(transcriptomeMap);
		transcriptome.fillSequenceTranscripts(genome);
	}
	public void annotate(String variantsFile,PrintStream out) throws IOException {
		VCFFileReader in = null;
		try {
			VCFFileWriter writer = new VCFFileWriter();
			in = new VCFFileReader(variantsFile);
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
		} finally {
			if(in!=null) in.close();
		}
		out.flush();
	}

	public void annotate(VCFRecord record) {
		GenomicVariant v = record.getVariant();
		List<VariantFunctionalAnnotation> annotations = transcriptome.calculateAnnotations(v, parameters);
		if(annotations.size()>0) {
			Collections.sort(annotations,new ImpactVariantFunctionalAnnotatorComparator());
			VariantFunctionalAnnotation ann = annotations.get(0);
			record.addAnnotation(new GenomicVariantAnnotation(v, GenomicVariantAnnotation.ATTRIBUTE_TRANSCRIPT_ANNOTATION, ann.getAnnotation()));
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
class ImpactVariantFunctionalAnnotatorComparator implements Comparator<VariantFunctionalAnnotation> {

	private Map<String, Integer> prioritiesMap = new HashMap<>();
	public ImpactVariantFunctionalAnnotatorComparator () {
		prioritiesMap = new HashMap<>();
		
		prioritiesMap.put(VariantFunctionalAnnotation.ANNOTATION_SPLICE_DONOR, 1);
		prioritiesMap.put(VariantFunctionalAnnotation.ANNOTATION_SPLICE_ACCEPTOR, 2);
		prioritiesMap.put(VariantFunctionalAnnotation.ANNOTATION_SPLICE_REGION, 3);
		prioritiesMap.put(VariantFunctionalAnnotation.ANNOTATION_FRAMESHIFT, 4);
		prioritiesMap.put(VariantFunctionalAnnotation.ANNOTATION_NONSENSE, 5);
		prioritiesMap.put(VariantFunctionalAnnotation.ANNOTATION_START_LOSS, 6);
		
		prioritiesMap.put(VariantFunctionalAnnotation.ANNOTATION_MISSENSE, 7);
		prioritiesMap.put(VariantFunctionalAnnotation.ANNOTATION_INFRAME_DEL, 8);
		prioritiesMap.put(VariantFunctionalAnnotation.ANNOTATION_INFRAME_INS, 9);
		prioritiesMap.put(VariantFunctionalAnnotation.ANNOTATION_STOP_LOSS, 10);
		prioritiesMap.put(VariantFunctionalAnnotation.ANNOTATION_SYNONYMOUS, 11);
		
		prioritiesMap.put(VariantFunctionalAnnotation.ANNOTATION_CODING, 12);
		prioritiesMap.put(VariantFunctionalAnnotation.ANNOTATION_5P_UTR, 13);
		prioritiesMap.put(VariantFunctionalAnnotation.ANNOTATION_3P_UTR, 14);
		prioritiesMap.put(VariantFunctionalAnnotation.ANNOTATION_NONCODINGRNA, 15);
		
		prioritiesMap.put(VariantFunctionalAnnotation.ANNOTATION_UPSTREAM, 16);
		prioritiesMap.put(VariantFunctionalAnnotation.ANNOTATION_DOWNSTREAM, 17);
		prioritiesMap.put(VariantFunctionalAnnotation.ANNOTATION_INTRON, 18);
		prioritiesMap.put(VariantFunctionalAnnotation.ANNOTATION_INTERGENIC, 19);
	}
	@Override
	public int compare(VariantFunctionalAnnotation a1, VariantFunctionalAnnotation a2) {
		int p1 = prioritiesMap.get(a1.getAnnotation());
		int p2 = prioritiesMap.get(a2.getAnnotation());
		if(p1!=p2) return p1 - p2;
		Transcript t1 = a1.getTranscript();
		Transcript t2 = a2.getTranscript();
		if(t1!=null && t2!=null) return t1.getId().compareTo(t2.getId());
		return 0;
	}
	
}
