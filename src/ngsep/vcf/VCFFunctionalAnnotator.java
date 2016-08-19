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
import java.util.Iterator;
import java.util.List;
import java.util.logging.Logger;

import ngsep.genome.ReferenceGenome;
import ngsep.main.CommandsDescriptor;
import ngsep.main.ProgressNotifier;
import ngsep.transcriptome.Transcriptome;
import ngsep.transcriptome.io.GFF3TranscriptomeHandler;
import ngsep.variants.GenomicVariant;
import ngsep.variants.GenomicVariantAnnotation;

public class VCFFunctionalAnnotator {
	private Logger log = Logger.getLogger(VCFFunctionalAnnotator.class.getName());
	public static final int DEF_UPSTREAM=1000;
	public static final int DEF_DOWNSTREAM=300;
	private Transcriptome transcriptome;
	private int offsetUpstream = DEF_UPSTREAM;
	private int offsetDownstream = DEF_DOWNSTREAM;
	
	private ProgressNotifier progressNotifier=null;
	
	public static void main(String[] args) throws Exception {
		if (args.length == 0 || args[0].equals("-h") || args[0].equals("--help")){
			CommandsDescriptor.getInstance().printHelp(VCFFunctionalAnnotator.class);
			return;
		}
		VCFFunctionalAnnotator annotator = new VCFFunctionalAnnotator();
		int i=0;
		while(i<args.length && args[i].charAt(0)=='-') {
			if("-u".equals(args[i])) {
				i++;
				annotator.setOffsetUpstream(Integer.parseInt(args[i]));
			} else if("-d".equals(args[i])) {
				i++;
				annotator.setOffsetDownstream(Integer.parseInt(args[i]));
			} else {
				System.err.println("Unrecognized option :" + args[i]);
				CommandsDescriptor.getInstance().printHelp(VCFFunctionalAnnotator.class);
				return;
			}
			i++;
		}
		String variantsFile = args[i++];
		String transcriptomeMap = args[i++];
		String sequenceFasta = args[i++]; 
		
		annotator.setLog(Logger.getLogger(VCFFunctionalAnnotator.class.getName()));
		annotator.loadMap(transcriptomeMap, new ReferenceGenome(sequenceFasta));
		annotator.annotate(variantsFile, System.out);
	}

	
	
	public int getOffsetUpstream() {
		return offsetUpstream;
	}



	public void setOffsetUpstream(int offsetUpstream) {
		this.offsetUpstream = offsetUpstream;
	}



	public int getOffsetDownstream() {
		return offsetDownstream;
	}



	public void setOffsetDownstream(int offsetDownstream) {
		this.offsetDownstream = offsetDownstream;
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
		List<GenomicVariantAnnotation> annotations = transcriptome.calculateAnnotation(v, offsetUpstream, offsetDownstream);
		for(GenomicVariantAnnotation ann:annotations) record.addAnnotation(ann);
	}
	
}
