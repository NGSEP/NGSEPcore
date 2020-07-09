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
package ngsep.discovery.rd;

import java.io.IOException;
import java.io.PrintStream;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.logging.Logger;

import JSci.maths.SpecialMath;
import JSci.maths.statistics.TDistribution;
import ngsep.genome.GenomicRegionComparator;
import ngsep.genome.ReferenceGenome;
import ngsep.main.io.ParseUtils;
import ngsep.math.PhredScoreHelper;
import ngsep.variants.CalledCNV;
import ngsep.variants.GenomicVariant;
import ngsep.variants.GenomicVariantImpl;
import ngsep.variants.io.GFFVariantsFileHandler;


public class CNVnatorReadDepthAlgorithm implements SingleSampleReadDepthAlgorithm {
	
	public static final String SOURCE_CNVNATOR = "CNVnator";
	private Logger log = Logger.getLogger(CNVnatorReadDepthAlgorithm.class.getName());
	
	private ReadDepthDistribution rdDistribution;
	private byte normalPloidy = 2;
	
	//Parameters of the distribution of read depth. Can be different than those 
	//calculated in the ReadDepthDistribution object if desired by the user
	private double meanReadDepth=0;
	private double sigmaReadDepth = 1;
	private long genomeSize = 0;
	
	private boolean greedyMergeLevels = true;
	
	private double minGradientValue = 0;
	
	
	public static final double MAX_PVALUE_REGION=0.05;
	
	public static void main(String[] args) throws Exception {
		CNVnatorReadDepthAlgorithm cnvNator = new CNVnatorReadDepthAlgorithm();
		int i=0;
		int fileType=0; 
		while(i<args.length && args[i].charAt(0)=='-') {
			if("-genomeSize".equals(args[i])) {
				i++;
				cnvNator.setGenomeSize(Long.parseLong(args[i]));
			} else if("-ploidy".equals(args[i])) {
				i++;
				cnvNator.setNormalPloidy(Byte.parseByte(args[i]));
			} else if("-meanRD".equals(args[i])) {
				i++;
				cnvNator.meanReadDepth = Double.parseDouble(args[i]);
			} else if("-sdRD".equals(args[i])) {
				i++;
				cnvNator.sigmaReadDepth = Double.parseDouble(args[i]);
			} else if("-minG".equals(args[i])) {
				i++;
				cnvNator.minGradientValue = Double.parseDouble(args[i]);
			} else if("-t".equals(args[i])) {
				i++;
				fileType = Integer.parseInt(args[i]);
			} else if("-noGreedyMerge".equals(args[i])) {
				cnvNator.greedyMergeLevels = false;
			}
			i++;
		}
		String reference = args[i++];
		String inFile = args[i++];
		String outFile = args[i++];
		String outFilePartition = outFile+".partition";
		
		cnvNator.log.info("Loading genome reference");
		ReferenceGenome genome = new ReferenceGenome(reference);
		cnvNator.log.info("Loaded genome reference. Sequences: "+genome.getNumSequences());
		if(fileType == 0) {
			cnvNator.log.info("Loading bins");
			cnvNator.rdDistribution = new ReadDepthDistribution(genome, ReadDepthDistribution.DEFAULT_BIN_SIZE);
			cnvNator.log.info("Loaded bins. Genome size assembly: "+cnvNator.rdDistribution.getGenomeSize());
			cnvNator.log.info("Processing alignments file: "+inFile);
			cnvNator.rdDistribution.processAlignments(inFile);
			cnvNator.log.info("Processed alignments file: "+inFile);
			cnvNator.rdDistribution.correctDepthByGCContent();
			cnvNator.log.info("Corrected GC-Content");
		} else {
			cnvNator.rdDistribution = new ReadDepthDistribution(inFile,fileType==2);
			cnvNator.log.info("Loaded bins from file: "+inFile+". Genome size bins file: "+cnvNator.rdDistribution.getGenomeSize());
		}
		if(cnvNator.genomeSize==0) cnvNator.setGenomeSize(cnvNator.rdDistribution.getGenomeSize());
		if(cnvNator.meanReadDepth==0) {
			cnvNator.log.info("Calculating read depth parameters");
			cnvNator.rdDistribution.calculateReadDepthDistParameters();
			cnvNator.meanReadDepth = cnvNator.rdDistribution.getMeanReadDepth();
			cnvNator.sigmaReadDepth = cnvNator.rdDistribution.getSigmaReadDepth();
			cnvNator.log.info("Calculated read depth parameters. Mean read depth: "+cnvNator.meanReadDepth+". Standard deviation: "+cnvNator.sigmaReadDepth);
		}
		
		List<CalledCNV> cnvs;
		if(fileType!=2) {
			cnvs = cnvNator.callCNVs();
		} else {
			cnvs = cnvNator.makeCNVCallsPartition();
		}
		cnvNator.log.info("Called "+cnvs.size()+" cnvs");
		cnvNator.printPartition(outFilePartition);
		cnvNator.log.info("Saved partition");
		GFFVariantsFileHandler cnvFH = new GFFVariantsFileHandler();
		PrintStream out = new PrintStream(outFile);
		cnvFH.saveVariants(cnvs, out);
		out.flush();
		out.close();
		cnvNator.log.info("Saved cnvs");
	}
	
	public void setReadDepthDistribution(ReadDepthDistribution rdDistribution) {
		this.rdDistribution = rdDistribution;
		this.meanReadDepth = rdDistribution.getMeanReadDepth();
		this.sigmaReadDepth = rdDistribution.getSigmaReadDepth();
	}
	
	public List<CalledCNV> callCNVs () {
		log.info("Calculating partition");
		calculatePartition(128, false);
		log.info("Calculated partition");
		log.info("Calling CNVs");
		return makeCNVCallsPartition();
	}
	
	

	
	private void calculatePartition (int range, boolean skipMasked) {
		List<String> sequenceNames = rdDistribution.getSequences().getNamesStringList();
		for(String seqName:sequenceNames) {
			List<ReadDepthBin> seqBins = rdDistribution.getBins(seqName);
			log.info("Calculating partition for sequence: "+seqName+". Sequence bins: "+seqBins.size());
			boolean [] mask = new boolean[seqBins.size()];
			Arrays.fill(mask, false);
			
			for(int binBand=2;binBand<=range;binBand++) {
				//Restore original RD values
				for(int i=0;i<seqBins.size();i++) {
					if(!mask[i]) {
						ReadDepthBin bin = seqBins.get(i);
						bin.setReadDepthLevel(bin.getCorrectedReadDepth());
					}
				}
				calcLevels(seqBins, binBand, mask, skipMasked);
				calcLevels(seqBins, binBand, mask, skipMasked);
				calcLevels(seqBins, binBand, mask, skipMasked);
				updateMask(seqBins, mask);
				if(binBand>=8) binBand++;
				if(binBand>=16) binBand+=2;
				if(binBand>=32) binBand+=4;
				if(binBand>=64) binBand+=8;
				if(binBand>=128) binBand+=16;
				if(binBand>=256) binBand+=32;
				if(binBand>=512) binBand+=64;
			}
		}
	}
	private void calcLevels(List<ReadDepthBin> seqBins, int binBand, boolean [] mask, boolean skipMasked) {
		double [] gradients = new double[seqBins.size()];
		Arrays.fill(gradients, 0);
		double inv2BinBand = 1.0/(binBand*binBand);
		double mean4 = meanReadDepth/4;
		double sigma2 = 4/(sigmaReadDepth*sigmaReadDepth);
		double ms2 = meanReadDepth/(sigmaReadDepth*sigmaReadDepth);
		int window = 3*binBand;
		double [] expsWindow = new double [window+1];
		for(int i=0;i<expsWindow.length;i++) {
			expsWindow[i] = i*Math.exp(-0.5*i*i*inv2BinBand);
		}
		double [] hrs = new double[seqBins.size()];
		for(int b=0;b<seqBins.size();b++) {
			ReadDepthBin binB = seqBins.get(b);
			hrs[b] = sigma2;
			if(binB.getReadDepthLevel() >= mean4) hrs[b] = ms2/binB.getReadDepthLevel();
		}
		//Calculate gradients
		for(int b=0;b<seqBins.size();b++) {
			ReadDepthBin binB = seqBins.get(b);
			if(mask[b]) continue;
			int distance = 0;
			for(int i=b+1;i<seqBins.size();i++) {
				ReadDepthBin binI = seqBins.get(i);
				if(mask[i]) continue;
				distance++;
				double difLevel = binI.getReadDepthLevel()-binB.getReadDepthLevel();
				double val = -0.5*difLevel*difLevel;
				//if(b==1058132 || b == 1058133)System.out.println("Bin start: "+binB.getFirst()+" BinI: "+binI.getFirst()+" current grad: "+gradients[b]+" Dif level: "+difLevel+" val: "+val+" Level: "+binB.getReadDepthLevel()+" invB: "+hrs[b]+" ms2: "+ms2);
				gradients[b] += expsWindow[distance]*Math.exp(val*hrs[b]);
				gradients[i] -= expsWindow[distance]*Math.exp(val*hrs[i]);
				if(distance == window) break;
			}
		}
		//Average signal on each segment
		for(int b=0;b<seqBins.size();b++) {
			//Determine next segment
			if(mask[b]) continue;
			int bStart = b;
			//Change from original implementation to disallow breaks due to gradients close to zero
			if(skipMasked) {
				while (b < seqBins.size() && (gradients[b] >= -minGradientValue || mask[b])) b++;
				while (b < seqBins.size() && (gradients[b] <  minGradientValue || mask[b])) b++;
			} else {
				while (b < seqBins.size() && gradients[b] >= -minGradientValue && !mask[b]) b++;
				while (b < seqBins.size() && gradients[b] <  minGradientValue && !mask[b]) b++;
			}
			int bStop = --b;
			assert (bStart<=bStop);
			//Calculate average signal for the segment
			double nl = 0;
			int n = 0;
		    for (int i = bStart;i <= bStop;i++) {
		    	ReadDepthBin bin = seqBins.get(i);
		    	if (mask[i]) continue;
		    	nl += bin.getReadDepthLevel();
		    	n++;
		    }
		    //Update with the average signal
		    nl/=n;
			for (int i = bStart;i <= bStop;i++) {
		    	ReadDepthBin bin = seqBins.get(i);
		    	if (!mask[i]) bin.setReadDepthLevel(nl);
		    	//if(i==1058132 || i == 1058133)System.out.println("Bin start: "+bin.getFirst()+" Level: "+bin.getReadDepthLevel()+" mask: "+mask[i]+" gradient: "+gradients[i]);
		    }
		}
	}
	
	private void updateMask(List<ReadDepthBin> seqBins, boolean[] mask) {
		Arrays.fill(mask, false);
		List<LevelRegion> regions = calculateRegions(seqBins);
		for(int i=1;i<regions.size()-1;i++) {
			LevelRegion rL = regions.get(i-1);
			LevelRegion region = regions.get(i);
			LevelRegion rR = regions.get(i+1);
			if(region.nBins<=1) continue;
			boolean smallN = rL.nBins<=15 || region.nBins<=15 || rR.nBins<=15;	
			//System.out.println("Testing region: "+region.first+" - "+region.last+ " Left start: "+rL.first+" right stop: "+ rR.last+ " levels "+rL.level+" "+region.level+" "+region.level+" "+rR.level);
			if(calculatePValue(rL,region,genomeSize,smallN)>=0.01) continue;
			if(calculatePValue(rR,region,genomeSize,smallN)>=0.01) continue;
			//System.out.println("P value region: "+region.first+" - "+region.last+" "+calculatePValue(region,false));
			if(calculatePValue(region,false)>MAX_PVALUE_REGION) continue;
			//System.out.println("Masking region: "+region.first+" - "+region.last);
			for(int j=region.first;j<=region.last;j++) {
				mask[j] = true;
			}
		}
		
		
	}
	private List<LevelRegion> calculateRegions(List<ReadDepthBin> seqBins) {
		List<LevelRegion> regions = new ArrayList<LevelRegion>();
		double lastLevel = -1;
		int first = 0;
		int last = 0;
		int nBins = 0;
		double sum = 0;
		double sum2 = 0;
		for(int i=0;i<seqBins.size();i++) {
			ReadDepthBin bin = seqBins.get(i);
			if(lastLevel != -1 && Math.abs(bin.getReadDepthLevel()-lastLevel)>=0.01) {
				regions.add(new LevelRegion(first, last, nBins, lastLevel, sum, sum2));
				first = i;
				
				sum = sum2 = nBins = 0;
			}
			double signal = bin.getCorrectedReadDepth(); 
			sum+=signal;
			sum2+=(signal*signal);
			nBins++;
			last = i;
			lastLevel = bin.getReadDepthLevel();
		}
		regions.add(new LevelRegion(first, last, nBins, lastLevel, sum, sum2));
		return regions;
	}
	private double calculatePValue(LevelRegion region, boolean correctMT) {
		if(region.nBins<=1) return 0.5;
		TDistribution tDist = new TDistribution(region.nBins-1);
		double stDev = Math.sqrt(region.getVarianceRD());
		if(correctMT && stDev < 0.001) stDev = sigmaReadDepth*Math.sqrt(region.getAverageRD()/meanReadDepth);
		if(stDev < 0.001) stDev = 1;
		double x = (region.getAverageRD() - meanReadDepth)*Math.sqrt(region.nBins)/stDev;
		double p = tDist.cumulative(x);
		if (x > 0) p = 1 - p;
		if(correctMT) p *= (0.99*genomeSize/(rdDistribution.getBinSize()*region.nBins));
		return p;
	}
	private double calculateModifiedPValue(List<ReadDepthBin> seqBins, LevelRegion region, boolean deletion) {
		double max = -1;
		double min = -1;
		for(int i=region.first;i<=region.last;i++) {
			ReadDepthBin bin = seqBins.get(i);
			if(bin.getCorrectedReadDepth()>max) max = bin.getCorrectedReadDepth();
			if(min==-1 || bin.getCorrectedReadDepth()<min) min = bin.getCorrectedReadDepth();
		}
		double p;
		if(deletion) {
			double x = 0.707*(max - meanReadDepth)/sigmaReadDepth;
		    p = 0.5*(1 + SpecialMath.error(x));
		  } else {
		    double x = 0.707*(min - meanReadDepth)/sigmaReadDepth;
		    p = 0.5*(1 - SpecialMath.error(x));
		  }
		return 0.99*genomeSize*Math.pow(p, region.nBins);
	}
	private double calculatePValue(LevelRegion borderRegion, LevelRegion region, double scale, boolean smallN) {
		if(smallN) {
			double nsigma = 1.8*Math.sqrt(borderRegion.level/meanReadDepth)*sigmaReadDepth;
			if(Math.abs(borderRegion.level-region.level)<nsigma) return 0.5;
			else return 0;
		}
		double var1 = borderRegion.getVarianceRD();
		if(var1 == 0 ) var1 = 1;
		double var2 = region.getVarianceRD();
		if(var2 == 0 ) var2 = 1;
		double tmp1 = var1/borderRegion.nBins,tmp2 = var2/region.nBins;
		double s = Math.sqrt(tmp1 + tmp2);
		double t = (borderRegion.getAverageRD() - region.getAverageRD())/s;
		double tmp = (tmp1 + tmp2)*(tmp1 + tmp2)*(borderRegion.nBins - 1)*(region.nBins - 1);
		tmp /= (tmp1*tmp1*(region.nBins - 1) + tmp2*tmp2*(borderRegion.nBins - 1));
		int ndf = (int)(tmp + 0.5);
		TDistribution tDist = new TDistribution(ndf);
		double pValue = tDist.cumulative(t);
		if(t > 0 ) pValue = 1 - pValue;
		//Correction for multiple testing
		pValue *= scale/(rdDistribution.getBinSize()*(borderRegion.nBins+region.nBins));
		return pValue;
	}
	
	private void printPartition(String outFilePartition) throws IOException {
		DecimalFormat fmt = ParseUtils.ENGLISHFMT;
		PrintStream out = new PrintStream(outFilePartition);
		double cut = meanReadDepth / (2*normalPloidy);
		List<String> sequenceNames = rdDistribution.getSequences().getNamesStringList();
		for(String seqName:sequenceNames) {
			List<ReadDepthBin> seqBins = rdDistribution.getBins(seqName);
			List<LevelRegion> regions = calculateFinalRegions(seqBins, cut,true,false);
			for(LevelRegion region:regions) {
				for(int i=region.first;i<=region.last;i++) {
					ReadDepthBin bin = seqBins.get(i);
					out.print(""+seqName+"\t"+bin.getFirst()+"\t"+bin.getLast()+"\t"+fmt.format(100*bin.getGcContent()));
					out.print("\t"+fmt.format(bin.getRawReadDepth())+"\t"+fmt.format(bin.getCorrectedReadDepth())+"\t"+fmt.format(bin.getReadDepthLevel()));
					out.print("\t"+fmt.format(region.getAverageRD())+"\t"+fmt.format(Math.sqrt(region.getVarianceRD()))+"\t"+fmt.format(region.level));
					out.print("\t"+region.pValueMethod);
					out.print("\t"+region.pValueGaussian);
					out.println();
				}
			}
			out.flush();
		}
		out.close();
	}
	
	private List<CalledCNV> makeCNVCallsPartition() {
		List<CalledCNV> cnvs = new ArrayList<CalledCNV>();
		double cut = meanReadDepth / (2*normalPloidy);
		List<String> sequenceNames = rdDistribution.getSequences().getNamesStringList();
		for(String seqName:sequenceNames) {
			if(sequenceNames.size()<100) log.info("Calling CNVs for sequence "+seqName);
			List<ReadDepthBin> seqBins = rdDistribution.getBins(seqName);
			List<LevelRegion> regions = calculateFinalRegions(seqBins, cut,false,greedyMergeLevels);
			
			for(LevelRegion region:regions) {
				double pValue = region.getPValue();
				
				if(pValue <= MAX_PVALUE_REGION) {
					float copies = (float) (region.getAverageRD()*normalPloidy/meanReadDepth);
					GenomicVariantImpl cnv = new GenomicVariantImpl(seqName, seqBins.get(region.first).getFirst(), seqBins.get(region.last).getLast(),GenomicVariant.TYPE_CNV);
					CalledCNV calledCNV = new CalledCNV(cnv,copies);
					calledCNV.setSource(SOURCE_CNVNATOR);
					calledCNV.setTotalReadDepth((int)Math.round(region.getAverageRD()*region.nBins));
					calledCNV.setGenotypeQuality(PhredScoreHelper.calculatePhredScore(pValue));
					cnvs.add(calledCNV);
				}
			}
		}
		Collections.sort(cnvs,new GenomicRegionComparator(rdDistribution.getSequences()));
		return cnvs;
	}
	
	private List<LevelRegion> calculateFinalRegions(List<ReadDepthBin> seqBins, double cut, boolean includeNeutral, boolean mergeLevels) {
		if(mergeLevels)mergeLevels(seqBins,cut);
		List <LevelRegion> regions = calculateAdjustedRegions(seqBins,cut);
		mergeSmallMiddleRegions(regions,seqBins);
		addGaussianCNVs(regions,seqBins,cut);
		if(!includeNeutral) return regions;
		List <LevelRegion> allRegions = new ArrayList<LevelRegion>();
		for(int i=0;i<regions.size();i++) {
			LevelRegion region = regions.get(i);
			LevelRegion lastRegion = null;
			if(allRegions.size()>0) lastRegion = allRegions.get(allRegions.size()-1);
			LevelRegion middleRegion = null;
			if(lastRegion!=null) {
				if(lastRegion.last+1<region.first) {
					middleRegion = new LevelRegion(seqBins, lastRegion.last+1, region.first-1);
					middleRegion.pValueMethod = calculatePValue(middleRegion, true);
				}
			} else if(region.first>0) {
				middleRegion = new LevelRegion(seqBins, 0, region.first-1);
				middleRegion.pValueMethod = calculatePValue(middleRegion, true);
			}
			if(middleRegion!=null) allRegions.add(middleRegion);
			allRegions.add(region);
		}
		return allRegions;
	}
	private void mergeLevels(List<ReadDepthBin> seqBins, double cut) {
		while (true) {
			boolean change = false;
			List<LevelRegion> regions = calculateRegions(seqBins);
			int i=0;
			int i2=1;
			double lastLevel = -1;
			while(i2<regions.size()) {
				LevelRegion r1 = regions.get(i);
				LevelRegion r2 = regions.get(i2);
				//Try to merge region i with region i+1
				double v1 = Math.abs(r1.level-r2.level);
				//System.out.println("Trying to merge levels starting at "+r1.first+" and "+r2.first+". Level1: "+r1.level+". Level2: "+r2.level+" diff: "+v1);
				if(v1 < cut) {
					double v2 = v1+1;
					double v3 = v1 + 1;
					if(lastLevel>0) v2 = Math.abs(lastLevel-r1.level);
					if(i2<regions.size()-1) v3 = Math.abs(r2.level-regions.get(i2+1).level);
					if(v1 < v2 && v1 < v3) {		
						//Merge region i with region i2
						change = true;
						double newLevel = r1.level*r1.nBins;
						newLevel += r2.level*r2.nBins;
						newLevel /= (r1.nBins + r2.nBins);
						//System.out.println("Merging levels starting at "+r1.first+" and "+r2.first+". Level1: "+r1.level+". Level2: "+r2.level+" diff: "+v1+" diff2: "+v2+" diff3: "+v3+" new level: "+newLevel);
						for(int j=r1.first;j<=r2.last;j++) {
							seqBins.get(j).setReadDepthLevel(newLevel);
						}
						r1.last = r2.last;
						r1.nBins +=r2.nBins;
						r1.level = newLevel;
						r1.sum +=r2.sum;
						r1.sumSquares +=r2.sumSquares;
					} else {
						lastLevel = r1.level;
						i=i2;
					}
				} else {
					lastLevel = r1.level;
					i=i2;
				}
				i2++;
			}
			if(!change) {
				return; 
			}
		}
		
	}
	private List<LevelRegion> calculateAdjustedRegions(List<ReadDepthBin> seqBins,double cut) {
		List<LevelRegion> regions = new ArrayList<LevelRegion>();
		int first=0;
		int limitFirst = 0;
		while(first<seqBins.size()) {
			ReadDepthBin bin = seqBins.get(first);
			double level = bin.getReadDepthLevel();
			double diffLevel = level-meanReadDepth; 
			if(Math.abs(diffLevel)<cut) {
				first++;
				continue;
			}
			int last = first+1;
			for(;last < seqBins.size();last++) {
				ReadDepthBin bin2 = seqBins.get(last);
				double diffLevel2 = bin2.getReadDepthLevel()-meanReadDepth;
				if(Math.abs(diffLevel2)<cut || (diffLevel>0 != diffLevel2>0)) {
					break;
				}
			}
			last--;
			if(last==first) {
				first++;
				continue;
			}
			//System.out.println("Calculating adjusted region for "+first+" "+last+". Initial pValue: "+calculatePValue(new LevelRegion(seqBins, first, last), true));
			LevelRegion region = localAdjust(seqBins, first, last, limitFirst);
			//System.out.println("Calculated adjusted region limits "+region.first+" "+region.last+". pValue: "+region.pValueMethod);
			if(region.pValueMethod>MAX_PVALUE_REGION) {
				//first+=3;
				first = last+1;
				continue;
			}
			regions.add(region);
			first = region.last+1;
			limitFirst = first;
		}
		return regions;
	}
	private void mergeSmallMiddleRegions(List<LevelRegion> regions, List<ReadDepthBin> seqBins) {
		List<LevelRegion> newRegionsList = new ArrayList<LevelRegion>();
		while(true) {
			boolean change = false;
			for(int i=0;i<regions.size();i++) {
				LevelRegion region = regions.get(i);
				if(region.pValueMethod > MAX_PVALUE_REGION) System.err.println("ERROR: region to merge starting at "+region.first+" has a high pValue: "+region.pValueMethod);
				LevelRegion lastRegion = null;
				if(newRegionsList.size()>0) {
					lastRegion = newRegionsList.get(newRegionsList.size()-1);
				}
				
				//if(lastRegion!=null) System.out.println("Trying to merge regions with limits "+lastRegion.first+" "+lastRegion.last+" "+region.first+" "+region.last+" pvalues "+lastRegion.pValueMethod+" "+region.pValueMethod+" levels "+lastRegion.level+" "+region.level);
				//Regions to merge are not neutral
				boolean mergeRegions = lastRegion!=null && lastRegion.pValueMethod <= MAX_PVALUE_REGION && region.pValueMethod <= MAX_PVALUE_REGION;
				//Both neighbors have the same sign
				mergeRegions = mergeRegions && ((lastRegion.level>meanReadDepth) == (region.level>meanReadDepth));
				if(mergeRegions && lastRegion.last+1<region.first-1) {
					LevelRegion middleRegion = new LevelRegion(seqBins, lastRegion.last+1, region.first-1);
					//Middle region is smaller than their neighbors
					mergeRegions = middleRegion.nBins < lastRegion.nBins && middleRegion.nBins < region.nBins;
					if(mergeRegions) {
						double pValue1 = calculatePValue(lastRegion, middleRegion,0.01*genomeSize,false);
						double pValue2 = calculatePValue(region, middleRegion,0.01*genomeSize,false);
						//The middle region is not significantly different than the neighbors
						mergeRegions = mergeRegions && (pValue1>=0.01 || pValue2>=0.01);
						//if(lastRegion!=null) System.out.println("Middle region "+middleRegion.first+" "+middleRegion.last+" pvalues "+pValue1+" "+pValue2);
					}
				}
				if(mergeRegions) {					
					LevelRegion combinedRegion = new LevelRegion(seqBins, lastRegion.first, region.last);
					combinedRegion.pValueMethod = calculatePValue(combinedRegion, true);
					if(combinedRegion.pValueMethod <=MAX_PVALUE_REGION) {
						//System.out.println("Merging regions with limits "+lastRegion.first+" "+lastRegion.last+" "+region.first+" "+region.last+" combined pValue: "+combinedRegion.pValueMethod);
						newRegionsList.set(newRegionsList.size()-1, combinedRegion);
						change = true;
					} else {
						mergeRegions = false;
					}
				} 
				if(!mergeRegions){
					newRegionsList.add(region);
				}
			}
			if(change) {
				regions.clear();
				regions.addAll(newRegionsList);
				newRegionsList.clear();
			} else break;
		}
	}

	private LevelRegion localAdjust(List<ReadDepthBin> seqBins, int first, int last, int limitFirst) { 
		double pValue = calculatePValue(new LevelRegion(seqBins, first, last),true);
		for (int i=0;i<1000 && pValue > MAX_PVALUE_REGION && last > first +1;i++) {
			int bestIndex = 0;
			if(first>limitFirst) {
				LevelRegion region2 = new LevelRegion(seqBins, first-1, last);
				double pVal2 = calculatePValue(region2, true);
				if(pVal2<pValue) {
					pValue = pVal2;
					bestIndex = 1;
				}
			}
			if(last+1<seqBins.size()) {
				LevelRegion region2 = new LevelRegion(seqBins, first, last+1);
				double pVal2 = calculatePValue(region2, true);
				if(pVal2<pValue) {
					pValue = pVal2;
					bestIndex = 2;
				}
			}
			if(first+1<last) {
				LevelRegion region2 = new LevelRegion(seqBins, first+1, last);
				double pVal2 = calculatePValue(region2, true);
				if(pVal2<pValue) {
					pValue = pVal2;
					bestIndex = 3;
				}
				region2 = new LevelRegion(seqBins, first, last-1);
				pVal2 = calculatePValue(region2, true);
				if(pVal2<pValue) {
					pValue = pVal2;
					bestIndex = 4;
				}
			}
			if(bestIndex==0) break;
			if(bestIndex ==1) {
				first--;
			} else if (bestIndex == 2) {
				last++;
			} else if (bestIndex == 3) {
				first++;
			} else if (bestIndex == 4) {
				last--;
			}
		}
		LevelRegion region = new LevelRegion(seqBins, first, last);
		region.pValueMethod = pValue;
		return region;
	}
	
	private void addGaussianCNVs(List<LevelRegion> regions,List<ReadDepthBin> seqBins,double cut) {
		List<LevelRegion> finalRegions = new ArrayList<LevelRegion>();
		int lastEnd = -1;
		for(int i=0;i<regions.size();i++) {
			LevelRegion region = regions.get(i);
			int last = region.first-1;
			finalRegions.addAll(lookForGaussian(seqBins,lastEnd+1,last,cut));
			finalRegions.add(region);
			lastEnd = region.last;
		}
		if(lastEnd<seqBins.size()-1)finalRegions.addAll(lookForGaussian(seqBins,lastEnd+1,seqBins.size()-1,cut));
		regions.clear();
		regions.addAll(finalRegions);
	}
	
	private List<LevelRegion> lookForGaussian(List<ReadDepthBin> seqBins, int first, int last,double cut) {
		List<LevelRegion> gaussianDeletions = new ArrayList<LevelRegion>();
		double min = meanReadDepth-cut;
		for(int i=first;i<=last;i++) {
			ReadDepthBin bS = seqBins.get(i);
			if(bS.getReadDepthLevel()<min) {
				int j = i+1;
				while (j<=last && seqBins.get(j).getReadDepthLevel()<min) j++;
				j--;
				if(i<j) {
					LevelRegion r = new LevelRegion(seqBins, i, j);
					r.pValueMethod = calculatePValue(r, true);
					r.pValueGaussian = calculateModifiedPValue(seqBins, r, true);
					if(r.pValueGaussian < MAX_PVALUE_REGION) gaussianDeletions.add(r);
					i=r.last;
				}
				
			}
		}
		return gaussianDeletions;
	}

	public byte getNormalPloidy() {
		return normalPloidy;
	}
	public void setNormalPloidy(byte normalPloidy) {
		this.normalPloidy = normalPloidy;
	}
	public double getMeanReadDepth() {
		return meanReadDepth;
	}
	public double getSigmaReadDepth() {
		return sigmaReadDepth;
	}
	public long getGenomeSize() {
		return genomeSize;
	}
	public void setGenomeSize(long genomeSize) {
		this.genomeSize = genomeSize;
	}
	public Logger getLog() {
		return log;
	}
	public void setLog(Logger log) {
		if (log == null) throw new NullPointerException("Log can not be null");
		this.log = log;
	}
	
	
	
}
class LevelRegion {
	int first;
	int last;
	int nBins;
	double level;
	double sum;
	double sumSquares;
	double pValueMethod=1;
	double pValueGaussian=1;
	public LevelRegion(int first, int last, int nBins, double level, double sum, double sumSquares) {
		super();
		this.first = first;
		this.last = last;
		this.nBins = nBins;
		this.level = level;
		this.sum = sum;
		this.sumSquares = sumSquares;
	}
	public LevelRegion(List<ReadDepthBin> seqBins, int first, int last) {
		this.first = first;
		this.last = last;
		this.sum = 0;
		this.sumSquares = 0;
		this.level = 0;
		nBins = 0;
		for(int i=first;i<=last;i++) {
			ReadDepthBin bin = seqBins.get(i);
			double signal = bin.getCorrectedReadDepth(); 
			sum+=signal;
			sumSquares+=signal*signal;
			nBins++;
			level += bin.getReadDepthLevel(); 
		}
		level/=nBins;
	}
	public double getAverageRD() {
		return sum/nBins;
	}
	public double getVarianceRD() {
		//if(nBins==1) return 0;
		//return (sumSquares-sum*sum/nBins)/(nBins-1);
		//As implemented in CNVNator
		double sum2 = sum*sum;
		double nBins2 = nBins*nBins;
		return sumSquares/nBins-sum2/nBins2; 
	}
	public double getPValue() {
		return Math.min(pValueMethod, pValueGaussian);
	}
	
}
class GaussFunction {
    public double function(double[] p, double[] x){
    	double a = (x[0]-p[0])/p[1];
        return (p[2]/(p[1]*Math.sqrt(2.0D*Math.PI)))*Math.exp(-0.5D*a*a);
    }
}