package ngsep.transcriptome;

import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;
import java.util.logging.Logger;

import ngsep.genome.ReferenceGenome;
import ngsep.main.ProgressNotifier;
import ngsep.sequences.DNASequence;
import ngsep.transcriptome.io.GFF3TranscriptomeHandler;
import ngsep.variants.CalledGenomicVariant;
import ngsep.variants.CalledSNV;
import ngsep.variants.GenomicVariant;
import ngsep.variants.SNV;
import ngsep.vcf.VCFFileReader;
import ngsep.vcf.VCFRecord;

public class MutatedPeptidesExtractor {

	private Logger log = Logger.getLogger(MutatedPeptidesExtractor.class.getName());
	private ProgressNotifier progressNotifier=null;
	
	public static final int DEF_MAX_LENGTH = 10;
	
	private int maxLength=DEF_MAX_LENGTH;
	private String mutatedSampleId;
	private String controlSampleId;
	private Transcriptome transcriptome;
	private ProteinTranslator translator = new ProteinTranslator();
	
	public static void main(String[] args) throws Exception {
		MutatedPeptidesExtractor instance = new MutatedPeptidesExtractor();
		String vcfFile = args[0];
		String transcriptomeMap = args[1];
		String referenceFile = args[2];
		instance.loadTranscriptome(transcriptomeMap, new ReferenceGenome(referenceFile));
		instance.findMutatedPeptides(vcfFile, System.out);
	}
	
	/**
	 * @return the log
	 */
	public Logger getLog() {
		return log;
	}



	/**
	 * @param log the log to set
	 */
	public void setLog(Logger log) {
		this.log = log;
	}



	/**
	 * @return the progressNotifier
	 */
	public ProgressNotifier getProgressNotifier() {
		return progressNotifier;
	}

	/**
	 * @param progressNotifier the progressNotifier to set
	 */
	public void setProgressNotifier(ProgressNotifier progressNotifier) {
		this.progressNotifier = progressNotifier;
	}

	/**
	 * @return the maxLength
	 */
	public int getMaxLength() {
		return maxLength;
	}

	/**
	 * @param maxLength the maxLength to set
	 */
	public void setMaxLength(int maxLength) {
		this.maxLength = maxLength;
	}
	
	public void setMaxLength(Integer maxLength) {
		this.setMaxLength(maxLength.intValue());
	}
	
	/**
	 * @return the mutatedSampleId
	 */
	public String getMutatedSampleId() {
		return mutatedSampleId;
	}

	/**
	 * @param mutatedSampleId the mutatedSampleId to set
	 */
	public void setMutatedSampleId(String mutatedSampleId) {
		this.mutatedSampleId = mutatedSampleId;
	}

	/**
	 * @return the controlSampleId
	 */
	public String getControlSampleId() {
		return controlSampleId;
	}

	/**
	 * @param controlSampleId the controlSampleId to set
	 */
	public void setControlSampleId(String controlSampleId) {
		this.controlSampleId = controlSampleId;
	}

	public void loadTranscriptome(String transcriptomeMap, ReferenceGenome genome) throws IOException {
		GFF3TranscriptomeHandler handler = new GFF3TranscriptomeHandler(genome.getSequencesMetadata());
		handler.setLog(log);
		transcriptome = handler.loadMap(transcriptomeMap);
		transcriptome.fillSequenceTranscripts(genome, log);
	}
	
	public void findMutatedPeptides (String vcfFile, PrintStream out) throws IOException  {
		
		try (VCFFileReader reader = new VCFFileReader(vcfFile)){
			//Ids of the samples
			List<String> sampleIds = reader.getSampleIds();
			if(sampleIds.size()==0) return;
			int mutSampleIdx = -1;
			int controlSampleIdx = -1;
			
			if(mutatedSampleId!=null) mutSampleIdx = sampleIds.indexOf(mutatedSampleId);
			if(controlSampleId!=null) controlSampleIdx = sampleIds.indexOf(controlSampleId);
			
			
			//If you want the full header
			List<SNV> selectedSNVs = new ArrayList<>();
			String lastSeqId = null;
			
			Iterator<VCFRecord> it = reader.iterator();
			while(it.hasNext()) {
				VCFRecord record = it.next();
				//Basic information from the variant
				GenomicVariant variant = record.getVariant();
				if(!(variant  instanceof SNV))continue;
				
				if(!variant.getSequenceName().equals(lastSeqId)) {
					if(lastSeqId!=null) findMutatedPeptides(selectedSNVs,out);
					selectedSNVs.clear();
					lastSeqId = variant.getSequenceName();
				}
				List<CalledGenomicVariant> genotypeCalls = record.getCalls();
				CalledGenomicVariant call = genotypeCalls.get(0);
				if(mutSampleIdx>=0) call = genotypeCalls.get(mutSampleIdx);
				if((call instanceof CalledSNV) && call.isHeterozygous()) {
					CalledGenomicVariant controlCall = null;
					if(controlSampleIdx>=0) controlCall = genotypeCalls.get(controlSampleIdx);
					if(controlCall==null || controlCall.isHomozygousReference()) selectedSNVs.add((SNV)variant);
				}
			}
			if(lastSeqId!=null) findMutatedPeptides(selectedSNVs,out);
		}
	}
	
	public void findMutatedPeptides(List<SNV> selectedSNVs, PrintStream out) {
		//Find candidates for report
		List<ReportRecord> candidates = new ArrayList<ReportRecord>();
		for(int i=0;i<selectedSNVs.size();i++) {
			SNV snv = selectedSNVs.get(i);
			List<ReportRecord> snvCandidates = getCandidates(snv, i);
			candidates.addAll(snvCandidates);
		}
		//Perform epitope predictions
		for(ReportRecord record:candidates) {
			SNV snv = record.getSnv();
			Transcript t = record.getTranscript();
			out.println(">"+snv.getSequenceName()+"_"+snv.getPosition()+"_TID_"+t.getId()+"_REF");
			out.println(record.getReferenceProtein());
			out.println();
			out.println(">"+snv.getSequenceName()+"_"+snv.getPosition()+"_TID_"+t.getId()+"_MUT");
			out.println(record.getMutatedProtein());
			out.println();
		}
	}
	private List<ReportRecord> getCandidates(SNV snv, int posList) {
		
		List <SNV> snvsZeroAllele = new ArrayList<SNV>();
		List <SNV> snvsOneAllele = new ArrayList<SNV>();
		snvsOneAllele.add(snv);
		//TODO: Add homozygous not reference and close heterozygous in phase  
		/*List <CalledSNV> nearbySNVS = getNearbySnvs(block,posList);
		for(CalledSNV nearSNV:nearbySNVS) {
			if(nearSNV.isHomozygous()) {
				snvsZeroAllele.add(nearSNV);
				snvsOneAllele.add(nearSNV);
			}
		}*/
		List<ReportRecord> candidates = new ArrayList<ReportRecord>();
		Set<String> cdnaSeqs = new TreeSet<String>();
		List<Transcript> transcripts = transcriptome.getTranscripts(snv.getSequenceName(), snv.getPosition()).asList(); 
		if(transcripts.size()==0) {
			return candidates;
		}
		for(Transcript t:transcripts) {
			int relPos = t.getRelativeTranscriptPosition(snv.getPosition());
			if(relPos<0) {
				continue;
			}
			String cDNASeqRef = t.getCDNASequence().toString();
			int beginPos = relPos-3*(maxLength-1);
			int endPos = relPos+3*maxLength;
			if(beginPos<0) {
				beginPos = 0;
			}
			beginPos-=beginPos%3;
			if(endPos>cDNASeqRef.length()) {
				endPos = cDNASeqRef.length();
			}
			endPos-=endPos%3;
			cDNASeqRef = cDNASeqRef.substring(beginPos,endPos);
			if(cdnaSeqs.contains(cDNASeqRef)) continue;
			cdnaSeqs.add(cDNASeqRef);
			//System.out.println("Sequence  ref: "+cDNASeqRef);
			char [] cdnaSeqZeroAllele = cDNASeqRef.toCharArray();
			applySNVs(snvsZeroAllele,t,beginPos,false,cdnaSeqZeroAllele);
			//applySNVs(triallelicZero,t,beginPos,false,cdnaSeqZeroAllele);
			//applySNVs(triallelicOne,t,beginPos,true,cdnaSeqZeroAllele);
			String seqZeroAllele = new String(cdnaSeqZeroAllele);
			//System.out.println("Sequence zero: "+seqZeroAllele);
			char [] cdnaSeqOneAllele = cDNASeqRef.toCharArray();
			applySNVs(snvsOneAllele,t,beginPos,false,cdnaSeqOneAllele);
			//applySNVs(triallelicZero,t,beginPos,true,cdnaSeqOneAllele);
			//applySNVs(triallelicOne,t,beginPos,false,cdnaSeqOneAllele);
			String seqOneAllele=new String(cdnaSeqOneAllele);
			//System.out.println("Sequence  one: "+seqOneAllele);
			String referenceProtein = getProtein(seqZeroAllele, 0);
			String modifiedProtein = getProtein(seqOneAllele, 0);
			if(!referenceProtein.equals(modifiedProtein)) {
				ReportRecord record = new ReportRecord(snv,t);
				record.setReferenceProtein(referenceProtein);
				record.setMutatedProtein(modifiedProtein);
				record.setTranscriptPosition(relPos/3+1);
				candidates.add(record);
			}
		}
		return candidates;
	}
	private void applySNVs(List<SNV> snvs, Transcript transcript, int cDNARefSubseqStart, boolean triallelicFirst, char [] cdnaSeqMod) {
		for(SNV snv: snvs) { 
			applySNV(snv, transcript, cDNARefSubseqStart,triallelicFirst, cdnaSeqMod);
		}
	}
	private void applySNV(SNV snv, Transcript transcript, int cDNARefSubseqStart, boolean triallelicFirst, char [] cdnaSeqMod) {
		int relPos = transcript.getRelativeTranscriptPosition(snv.getPosition())-cDNARefSubseqStart;
		if(relPos>=0 && relPos<cdnaSeqMod.length) {
			char ref = snv.getReferenceBase();
			char alt = snv.getAlternativeBase();
			//TODO: Handle triallelic
			/*if(snv.isTriallelic()) {
				alt = snv.getGenotype().charAt(1);
				if(triallelicFirst) {
					alt = snv.getGenotype().charAt(0);
				}
			}*/
			if(transcript.isNegativeStrand()) {
				ref = DNASequence.getComplement(ref);
				alt = DNASequence.getComplement(alt);
			}
			if(cdnaSeqMod[relPos]!=ref) {
				System.err.println("Unexpected reference base "+cdnaSeqMod[relPos]+" found. Expected: "+ref+" for SNV in chromosome: "+snv.getSequenceName()+" and position:"+snv.getPosition());
			}
			cdnaSeqMod[relPos]=alt;
		}
	}
	
	public String getProtein(String dnaSequence, int initialPosition) {
		return String.valueOf(translator.getProteinSequence(dnaSequence.toCharArray(),initialPosition, false));
	}
	public String getPeptide(String protein,int start, int length) {
		if(protein.length() > start) {
			return protein.substring(start,Math.min(protein.length(), start+length));
		}
		return "";
	}

}

class ReportRecord {
	private SNV snv;
	private Transcript transcript;
	private String referenceProtein="";
	private String mutatedProtein="";
	private int transcriptPosition;
	public ReportRecord (SNV snv, Transcript transcript) {
		this.snv = snv;
		this.transcript = transcript;
	}
	
	/**
	 * @return the snv
	 */
	public SNV getSnv() {
		return snv;
	}

	public Transcript getTranscript() {
		return transcript;
	}
	
	public String getReferenceProtein() {
		return referenceProtein;
	}
	public void setReferenceProtein(String referenceProtein) {
		if(referenceProtein == null) {
			this.referenceProtein = "";
		} else {
			this.referenceProtein = referenceProtein;
		}
	}
	public String getMutatedProtein() {
		return mutatedProtein;
	}
	public void setMutatedProtein(String mutatedProtein) {
		if(mutatedProtein==null) {
			this.mutatedProtein = "";
		} else {
			this.mutatedProtein = mutatedProtein;
		}
	}
	
	public int getTranscriptPosition() {
		return transcriptPosition;
	}
	public void setTranscriptPosition(int transcriptPosition) {
		this.transcriptPosition = transcriptPosition;
	}
}
