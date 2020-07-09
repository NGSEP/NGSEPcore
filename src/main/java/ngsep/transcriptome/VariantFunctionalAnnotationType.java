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
package ngsep.transcriptome;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * 
 * @author Jorge Duitama
 *
 */
public class VariantFunctionalAnnotationType {
	
	public static final String ANNOTATION_SPLICE_DONOR="splice_donor_variant";
	public static final String ANNOTATION_SPLICE_ACCEPTOR="splice_acceptor_variant";
	public static final String ANNOTATION_FRAMESHIFT="frameshift_variant";
	public static final String ANNOTATION_NONSENSE="stop_gained";
	public static final String ANNOTATION_START_LOST="start_lost";
	
	public static final String ANNOTATION_EXONIC_SPLICE_REGION="exonic_splice_region_variant";
	public static final String ANNOTATION_SPLICE_REGION="splice_region_variant";
	public static final String ANNOTATION_INFRAME_DEL="inframe_deletion";
	public static final String ANNOTATION_INFRAME_INS="inframe_insertion";
	public static final String ANNOTATION_STOP_LOST="stop_lost";
	public static final String ANNOTATION_MISSENSE="missense_variant";
	public static final String ANNOTATION_SYNONYMOUS="synonymous_variant";
	
	public static final String ANNOTATION_CODING="coding_sequence_variant";
	public static final String ANNOTATION_5P_UTR="5_prime_UTR_variant";
	public static final String ANNOTATION_3P_UTR="3_prime_UTR_variant";
	public static final String ANNOTATION_NONCODINGRNA="non_coding_transcript_exon_variant";
	
	public static final String ANNOTATION_UPSTREAM="upstream_transcript_variant";
	public static final String ANNOTATION_DOWNSTREAM="downstream_transcript_variant";
	public static final String ANNOTATION_INTRON="intron_variant";
	public static final String ANNOTATION_INTERGENIC="intergenic_variant";
	
	//Default name consistent with sequence ontology (http://www.sequenceontology.org)
	private String name;
	private String soAccession;
	private String ngsep2Name;
	private boolean coding;
	
	private static Map<String, Integer> annotationPriorities = null;
	private static Map<String, VariantFunctionalAnnotationType> annotationTypesByName = null;
	private static Map<String, VariantFunctionalAnnotationType> annotationTypesByNgsep2Name = null;
	private static Map<String, VariantFunctionalAnnotationType> annotationTypesBySOAccession = null;
	private VariantFunctionalAnnotationType(String name, String ngsep2Name, String soAccession, boolean coding) {
		super();
		this.name = name;
		this.soAccession = soAccession;
		this.ngsep2Name = ngsep2Name;
		this.coding = coding;
	}
	
	private static void loadTypes() {
		List<VariantFunctionalAnnotationType> types = new ArrayList<>();
		types.add(new VariantFunctionalAnnotationType(ANNOTATION_SPLICE_DONOR,null,"SO:0001575", true));
		types.add(new VariantFunctionalAnnotationType(ANNOTATION_SPLICE_ACCEPTOR,null,"SO:0001574", true));
		types.add(new VariantFunctionalAnnotationType(ANNOTATION_FRAMESHIFT,"Frameshift","SO:0001589",true));
		types.add(new VariantFunctionalAnnotationType(ANNOTATION_NONSENSE,"Nonsense","SO:0001587", true));
		types.add(new VariantFunctionalAnnotationType(ANNOTATION_START_LOST,null,"SO:0002012",true));
		
		types.add(new VariantFunctionalAnnotationType(ANNOTATION_EXONIC_SPLICE_REGION,"ExonJunction","SO:0002084",true));
		types.add(new VariantFunctionalAnnotationType(ANNOTATION_SPLICE_REGION,null,"SO:0001630",true));
		types.add(new VariantFunctionalAnnotationType(ANNOTATION_INFRAME_DEL,null,"SO:0001822",true));
		types.add(new VariantFunctionalAnnotationType(ANNOTATION_INFRAME_INS,null,"SO:0001821",true));
		types.add(new VariantFunctionalAnnotationType(ANNOTATION_STOP_LOST,null,"SO:0001578",true));
		types.add(new VariantFunctionalAnnotationType(ANNOTATION_MISSENSE,"Missense","SO:0001583",true));
		types.add(new VariantFunctionalAnnotationType(ANNOTATION_SYNONYMOUS,"Synonymous","SO:0001819",true));
	
		types.add(new VariantFunctionalAnnotationType(ANNOTATION_CODING,"Coding","SO:0001580",true));
		types.add(new VariantFunctionalAnnotationType(ANNOTATION_5P_UTR,"FivePrimeUTR","SO:0001623",true));
		types.add(new VariantFunctionalAnnotationType(ANNOTATION_3P_UTR,"ThreePrimeUTR","SO:0001624",true));
		types.add(new VariantFunctionalAnnotationType(ANNOTATION_NONCODINGRNA,"NCRNA","SO:0001792",false));
	
		types.add(new VariantFunctionalAnnotationType(ANNOTATION_UPSTREAM,"Upstream","SO:0001986",false));
		types.add(new VariantFunctionalAnnotationType(ANNOTATION_DOWNSTREAM,"Downstream","SO:0001987",false));
		types.add(new VariantFunctionalAnnotationType(ANNOTATION_INTRON,"Intron","SO:0001627",false));
		types.add(new VariantFunctionalAnnotationType(ANNOTATION_INTERGENIC,"Intergenic","SO:0001628",false));
		
		annotationPriorities = new HashMap<>();
		annotationTypesByName = new HashMap<>();
		annotationTypesByNgsep2Name = new HashMap<>();
		annotationTypesBySOAccession = new HashMap<>();
		for(int i=0;i<types.size();i++) {
			VariantFunctionalAnnotationType type = types.get(i);
			annotationPriorities.put(type.getName(), i);
			annotationTypesByName.put(type.getName(), type);
			if(type.getNgsep2Name()!=null) annotationTypesByNgsep2Name.put(type.getNgsep2Name(), type);
			annotationTypesBySOAccession.put(type.getSoAccession(), type);
			
		}
	}
	public static VariantFunctionalAnnotationType getTypeByName(String name) {
		if(annotationTypesByName==null) loadTypes();
		return annotationTypesByName.get(name);
	}
	public static VariantFunctionalAnnotationType getTypeByNgsep2Name(String name) {
		if(annotationTypesByName==null) loadTypes();
		return annotationTypesByNgsep2Name.get(name);
	}
	public static VariantFunctionalAnnotationType getTypeBySOAccession(String accessionId) {
		if(annotationTypesByName==null) loadTypes();
		return annotationTypesByNgsep2Name.get(accessionId);
	}
	
	public static boolean isTypeCoding(String name) {
		if(annotationTypesByName==null) loadTypes();
		if(name == null) throw new NullPointerException("Name can not be null");
		VariantFunctionalAnnotationType type = annotationTypesByName.get(name);
		if(type==null) return false;
		return type.isCoding();
	}
	public static boolean isTypeNonSynonymous(String name) {
		if(!isTypeCoding(name)) return false;
		return !name.equals(ANNOTATION_SYNONYMOUS);
	}
	
	public static VariantFunctionalAnnotationType getTypeBySearchKey(String key) {
		if(annotationTypesByName==null) loadTypes();
		VariantFunctionalAnnotationType answer = getTypeByName(key);
		if(answer == null) answer = getTypeBySOAccession(key);
		if(answer == null) answer = getTypeByNgsep2Name(key);
		return answer;
	}
	public static int getNumberSupportedTypes () {
		if(annotationTypesByName==null) loadTypes();
		return annotationTypesByName.size();
	}

	/**
	 * @return the name
	 */
	public String getName() {
		return name;
	}

	/**
	 * @return the soAccession
	 */
	public String getSoAccession() {
		return soAccession;
	}

	/**
	 * @return the ngsep2Name
	 */
	public String getNgsep2Name() {
		return ngsep2Name;
	}

	/**
	 * @return the coding
	 */
	public boolean isCoding() {
		return coding;
	}
	public static Comparator<VariantFunctionalAnnotation> getPriorityComparator () {
		return new Comparator<VariantFunctionalAnnotation>() {
			@Override
			public int compare(VariantFunctionalAnnotation a1, VariantFunctionalAnnotation a2) {
				if(annotationPriorities==null) loadTypes();
				Integer p1 = annotationPriorities.get(a1.getTypeName());
				Integer p2 = annotationPriorities.get(a2.getTypeName());
				if(p1==null) throw new IllegalArgumentException("Unrecognized annotation type: "+a1.getTypeName());
				if(p2==null) throw new IllegalArgumentException("Unrecognized annotation type: "+a2.getTypeName());
				if(p1!=p2) return p1 - p2;
				Transcript t1 = a1.getTranscript();
				Transcript t2 = a2.getTranscript();
				if(t1!=null && t2!=null) return t1.getId().compareTo(t2.getId());
				return 0;
			}
		};
	
	}
}
