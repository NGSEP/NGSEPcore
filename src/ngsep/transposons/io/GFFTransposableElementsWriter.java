package ngsep.transposons.io;

import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

import ngsep.genome.io.GFF3GenomicFeatureLine;
import ngsep.genome.io.GFF3Writer;
import ngsep.transposons.TransposableElementAnnotation;
import ngsep.transposons.TransposableElementFamily;

public class GFFTransposableElementsWriter {
	public void printGFF(List<TransposableElementAnnotation> transposonAnnotations, String outFilename) throws IOException {
		try(PrintStream out = new PrintStream(outFilename)) {
			printGFF(transposonAnnotations, out);
		}
		
	}
	public void printGFF(List<TransposableElementAnnotation> transposonAnnotations, PrintStream out) {
		GFF3Writer writer = new GFF3Writer();
		writer.printHeader(out);
		int i = 1;
		for(TransposableElementAnnotation t:transposonAnnotations) {
			List<GFF3GenomicFeatureLine> lines = new ArrayList<GFF3GenomicFeatureLine>();
			GFF3GenomicFeatureLine mainLine = new GFF3GenomicFeatureLine(t.getSequenceName(), t.getFirst(), t.getLast(), GFF3GenomicFeatureLine.FEATURE_TYPE_REPEAT_REGION);
			String id =  "TE_"+i;
			mainLine.setNegativeStrand(t.isNegativeStrand());
			mainLine.setSource("NGSEP");
			mainLine.setId(id);
			TransposableElementFamily family = t.getInferredFamily();
			if(family == null) family = t.getSourceFamily();
			if(family!=null) mainLine.addAnnotation("classification", family.toString());
			String ontologyId = writer.getOntologyId(GFF3GenomicFeatureLine.FEATURE_TYPE_REPEAT_REGION);
			if(ontologyId!=null) mainLine.addAnnotation(GFF3GenomicFeatureLine.ATTRIBUTE_ONTOLOGY, ontologyId);
			lines.add(mainLine);
			String method = t.getQueryName()!=null?"homology":"structural";
			mainLine.addAnnotation("method", method);
			int first2 = t.getFirst();
			String tsd = t.getTsd();
			if(tsd!=null) {
				mainLine.addAnnotation("tsd", tsd);
				first2 = t.getFirst()+tsd.length();
				GFF3GenomicFeatureLine tsdLine1 = new GFF3GenomicFeatureLine(t.getSequenceName(), t.getFirst(), first2-1, GFF3GenomicFeatureLine.FEATURE_TYPE_TSD);
				tsdLine1.setNegativeStrand(t.isNegativeStrand());
				tsdLine1.setSource("NGSEP");
				tsdLine1.setId(id+"_tsd1");
				tsdLine1.addAnnotation(GFF3GenomicFeatureLine.ATTRIBUTE_PARENT, id);
				ontologyId = writer.getOntologyId(GFF3GenomicFeatureLine.FEATURE_TYPE_TSD);
				if(ontologyId!=null) tsdLine1.addAnnotation(GFF3GenomicFeatureLine.ATTRIBUTE_ONTOLOGY, ontologyId);
				lines.add(tsdLine1);
			}
			if(t.getLeftEndRepeat()>first2) {
				String type = t.isPositiveStrand()?GFF3GenomicFeatureLine.FEATURE_TYPE_5PLTR:GFF3GenomicFeatureLine.FEATURE_TYPE_3PLTR;
				GFF3GenomicFeatureLine ltrLine1 = new GFF3GenomicFeatureLine(t.getSequenceName(), first2, t.getLeftEndRepeat(), type);
				first2=t.getLeftEndRepeat()+1;
				ltrLine1.setNegativeStrand(t.isNegativeStrand());
				ltrLine1.setSource("NGSEP");
				ltrLine1.setId(id+"_"+type);
				ltrLine1.addAnnotation(GFF3GenomicFeatureLine.ATTRIBUTE_PARENT, id);
				ontologyId = writer.getOntologyId(type);
				if(ontologyId!=null) ltrLine1.addAnnotation(GFF3GenomicFeatureLine.ATTRIBUTE_ONTOLOGY, ontologyId);
				lines.add(ltrLine1);
			}
			int lastCenter = t.getLast();
			if(tsd!=null) lastCenter-=tsd.length();
			if(t.getRightStartRepeat()>t.getLeftEndRepeat()) lastCenter=t.getRightStartRepeat()-1;
			GFF3GenomicFeatureLine centerLine = new GFF3GenomicFeatureLine(t.getSequenceName(), first2, lastCenter, GFF3GenomicFeatureLine.FEATURE_TYPE_TRANSPOSON);
			centerLine.setNegativeStrand(t.isNegativeStrand());
			centerLine.setSource("NGSEP");
			centerLine.setId(id+"_TR");
			centerLine.addAnnotation(GFF3GenomicFeatureLine.ATTRIBUTE_PARENT, id);
			ontologyId = writer.getOntologyId(GFF3GenomicFeatureLine.FEATURE_TYPE_TRANSPOSON);
			if(family!=null) centerLine.addAnnotation("classification", family.toString());
			if(ontologyId!=null) centerLine.addAnnotation(GFF3GenomicFeatureLine.ATTRIBUTE_ONTOLOGY, ontologyId);
			lines.add(centerLine);
			int lastLTR = t.getLast();
			if(tsd!=null) lastLTR-=tsd.length();
			if(t.getRightStartRepeat()>=0 && t.getRightStartRepeat()<t.getLast()) {
				String type = t.isPositiveStrand()?GFF3GenomicFeatureLine.FEATURE_TYPE_3PLTR:GFF3GenomicFeatureLine.FEATURE_TYPE_5PLTR;
				GFF3GenomicFeatureLine ltrLine2 = new GFF3GenomicFeatureLine(t.getSequenceName(), t.getRightStartRepeat(), lastLTR, type);
				first2=t.getLeftEndRepeat()+1;
				ltrLine2.setNegativeStrand(t.isNegativeStrand());
				ltrLine2.setSource("NGSEP");
				ltrLine2.setId(id+"_"+type);
				ltrLine2.addAnnotation(GFF3GenomicFeatureLine.ATTRIBUTE_PARENT, id);
				ontologyId = writer.getOntologyId(type);
				if(ontologyId!=null) ltrLine2.addAnnotation(GFF3GenomicFeatureLine.ATTRIBUTE_ONTOLOGY, ontologyId);
				lines.add(ltrLine2);
			}
			if(tsd!=null) {
				GFF3GenomicFeatureLine tsdLine1 = new GFF3GenomicFeatureLine(t.getSequenceName(), lastLTR+1, t.getLast(), GFF3GenomicFeatureLine.FEATURE_TYPE_TSD);
				tsdLine1.setNegativeStrand(t.isNegativeStrand());
				tsdLine1.setSource("NGSEP");
				tsdLine1.setId(id+"_tsd2");
				tsdLine1.addAnnotation(GFF3GenomicFeatureLine.ATTRIBUTE_PARENT, id);
				ontologyId = writer.getOntologyId(GFF3GenomicFeatureLine.FEATURE_TYPE_TSD);
				if(ontologyId!=null) tsdLine1.addAnnotation(GFF3GenomicFeatureLine.ATTRIBUTE_ONTOLOGY, ontologyId);
				lines.add(tsdLine1);
			}
			writer.printGffFeatureLines(lines, out);
			i++;
		}
	}
}
