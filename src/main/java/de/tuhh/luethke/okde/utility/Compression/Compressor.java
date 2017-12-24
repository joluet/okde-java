/*
 * Copyright 2014 Jonas Luethke
 */


package de.tuhh.luethke.okde.utility.Compression;

import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;
import java.util.List;

import org.ejml.simple.SimpleMatrix;

import de.tuhh.luethke.okde.Exceptions.EmptyDistributionException;
import de.tuhh.luethke.okde.model.BaseSampleDistribution;
import de.tuhh.luethke.okde.model.OneComponentDistribution;
import de.tuhh.luethke.okde.model.SampleModel;
import de.tuhh.luethke.okde.model.TwoComponentDistribution;
import de.tuhh.luethke.okde.utility.MomentMatcher;
import de.tuhh.luethke.okde.utility.Matrices.MatrixOps;
import de.tuhh.luethke.okde.utility.Projection.ProjectionData;
import de.tuhh.luethke.okde.utility.Projection.Projector;

public class Compressor {
	private final static double CONST_SMALL_TOLERANCE = 1E-10;
	
	//private final static double D_TH = 0.1;
	
	private static final double MIN_EM_DISTANCE = 2.34d;

	private static final float INC_TH_SCALE = 1.5f;
	private static final float DEC_TH_SCALE = 0.6f;
	private static final float CHECK_IF_DEC_SCALE = 0.5f;
	
	private static void setNoOfComponentsThreshold(SampleModel dist, int noOfCompsBeforeCompression, int noOfCompsAfterCompression) {
		float threshold = dist.getNoOfCompsThreshold();
		if (noOfCompsAfterCompression > threshold)
			threshold = threshold * INC_TH_SCALE;
		else if (noOfCompsAfterCompression <= threshold * CHECK_IF_DEC_SCALE)
			threshold = threshold * DEC_TH_SCALE;
		dist.setNoOfCompsThreshold(threshold);
	}
	
	public static boolean emUpdate(SampleModel dist, List<Integer> updatePoints) {
		ArrayList<BaseSampleDistribution> subDistributions = dist.getSubDistributions();
		ArrayList<SimpleMatrix> means = dist.getSubMeans();
		ArrayList<SimpleMatrix> smoothedCovariances = dist.getSubSmoothedCovariances();
		ArrayList<SimpleMatrix> covariances = dist.getSubCovariances();
		ArrayList<Double> weights = dist.getSubWeights();
		boolean pointMerged = false;
		int count = 0;
		for(int point : updatePoints) {
			for(int i=0; i<means.size() && i!=point && !pointMerged; i++){
				if(subDistributions.get(i).getClass() == TwoComponentDistribution.class) {
					TwoComponentDistribution subComponent = (TwoComponentDistribution)subDistributions.get(i);
					// calculate mahalanobis distance (x-m)L(x-m)' to each mean until one is small enough
					double md = means.get(point).minus( means.get(i) ).transpose().mult( smoothedCovariances.get(i).invert() ).mult( means.get(point).minus(means.get(i)) ).trace();
					if(md < MIN_EM_DISTANCE) {
						// just add the new point to sub model
						
						// which subcomponent is closest?
						OneComponentDistribution[] subSubComponents = subComponent.getSubComponents();
						double distance1 = euclidianDistance(subSubComponents[0].getGlobalMean(), means.get(point));
						double distance2 = euclidianDistance(subSubComponents[1].getGlobalMean(), means.get(point));

						int mergeId = 0;
						if(distance1 < distance2)
							mergeId = 0;
						else
							mergeId = 1;
						OneComponentDistribution componentToMerge = subSubComponents[mergeId];

							
						SimpleMatrix[] meansArray = { componentToMerge.getGlobalMean(), means.get(point) };
						SimpleMatrix[] covarianceArray = { componentToMerge.getGlobalCovariance(), covariances.get(point) };
						
						double subSubweight1 = componentToMerge.getGlobalWeight()*subComponent.getGlobalWeight();
						double subSubweight2 = weights.get(point);
						double globalWeight = subComponent.getGlobalWeight() + subSubweight2;
						double subSubWeightSum = subSubweight1 + subSubweight2;
						subSubweight1 /= subSubWeightSum;
						subSubweight2 /= subSubWeightSum;

						double[] weightsArray = { subSubweight1, subSubweight2 };
						
						
						OneComponentDistribution oneCompDist = null;
						try {
							TwoComponentDistribution twoCompDist = new TwoComponentDistribution(weightsArray, meansArray, covarianceArray, dist.getBandwidthMatrix());
							double subWeight1 = subSubComponents[0].getGlobalWeight()*subComponent.getGlobalWeight();
							double subWeight2 = subSubComponents[1].getGlobalWeight()*subComponent.getGlobalWeight();
							if(mergeId == 0)
								subWeight1 += weights.get(point);
							else
								subWeight2 += weights.get(point);
							double subWeightSum = subWeight1+subWeight2;
							subWeight1 /= subWeightSum;
							subWeight2 /= subWeightSum;
							
							MomentMatcher.matchMoments(twoCompDist);
							oneCompDist = new OneComponentDistribution(twoCompDist);
							subSubComponents[mergeId] = oneCompDist;
							subSubComponents[0].setGlobalWeight(subWeight1);
							subSubComponents[1].setGlobalWeight(subWeight2);
							MomentMatcher.matchMoments(subComponent);
							subComponent.setGlobalWeight(globalWeight);
							
							dist.mEMCount++;
							/*double compressionError = Hellinger.calculateUnscentedHellingerDistance(oneCompDist, twoCompDist);
							dist.mEMError = (dist.mEMError*dist.mEMCount + compressionError)/(dist.mEMCount + 1);
							dist.mEMCount++;*/
						} catch (EmptyDistributionException e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
						} catch (Exception e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
						}
						/*						
						Double[] subWeights = twoCompDist.getSubWeights();
						double newWeight1 = subWeights[0] / (subWeights[0] + subWeights[1]);
						double newWeight2 = subWeights[1] / (subWeights[0] + subWeights[1]);
						twoCompDist.getSubComponents()[0].setGlobalWeight(newWeight1);
						twoCompDist.getSubComponents()[1].setGlobalWeight(newWeight2);
						dist.getSubDistributions().set(i, twoCompDist);
						*/
						
						dist.getSubDistributions().remove(point);
						count++;
						pointMerged = true;
					}
				}
			}
		}
		return (count == updatePoints.size());
	}

	public static void compress(SampleModel dist, ArrayList<Integer> newComponents) throws Exception {
		// check wether compression is necessary using hysteresis rule
		if (dist.getSubMeans().size() <= dist.getNoOfCompsThreshold())
			return;
		// try em update
		boolean successfulEMUpdate = emUpdate(dist, newComponents);
		if(successfulEMUpdate) {
			return;
		}
		ProjectionData projectionData = null;
		try {
			projectionData = Projector.projectSampleDistToSubspace(dist);
		} catch(Exception e) {
			// if projection fails: stop compression
			System.out.println("Projection failed. Aborted Compression");
			return;
		}
		revitalizeComponents(dist);
		//System.out.println("COMPRESS");
		int noOfCompsBeforeCompression = dist.getSubMeans().size();
		SampleModel inputModelCopy = new SampleModel(dist);
		double compressionError = Double.MAX_VALUE;
		if(inputModelCopy.getSubDistributions().size() > 1)
			compressionError = mergeTwoClosestComps(inputModelCopy);
		while (compressionError < dist.mCompressionThreshold) {
			dist.overWirite(inputModelCopy);
			if(inputModelCopy.getSubDistributions().size() > 1)
				compressionError = mergeTwoClosestComps(inputModelCopy);
			else
				compressionError = Double.MAX_VALUE;
		}
		Projector.projectSampleDistToOriginalSpace(dist, projectionData);
		int noOfCompsAfterCompression = dist.getSubMeans().size();
		setNoOfComponentsThreshold(dist, noOfCompsBeforeCompression, noOfCompsAfterCompression);
	}

	private static void revitalizeComponents(SampleModel dist) throws Exception {
		// check which sub distributions have to be revitalized
		for (int i = 0; i < dist.getSubDistributions().size(); i++) {
			if (dist.getSubDistributions().get(i).getClass() == TwoComponentDistribution.class) {
				TwoComponentDistribution subDist = (TwoComponentDistribution) dist.getSubDistributions().get(i);
				double tmpWeight = subDist.getGlobalWeight();
				MomentMatcher.matchMoments(subDist);
				OneComponentDistribution oneCompDist = new OneComponentDistribution(subDist);
				double compressionError = Hellinger.calculateUnscentedHellingerDistance(oneCompDist, subDist);
				subDist.setGlobalWeight(tmpWeight);
				if (compressionError >= dist.mCompressionThreshold) {
					OneComponentDistribution subComp1 = subDist.getSubComponents()[0];
					OneComponentDistribution subComp2 = subDist.getSubComponents()[1];
					BaseSampleDistribution splitDist1 = null, splitDist2 = null;
					// check wether covariance of sub component is zero --> no splitting necessary
					if(subComp1.getGlobalCovariance().elementSum() > CONST_SMALL_TOLERANCE)
						splitDist1 = subComp1.split(tmpWeight);
					else{ 
						subComp1.scaleGlobalWeight(tmpWeight);
						splitDist1 = subComp1;
					}
					// check wether covariance of sub component is zero --> no splitting necessary
					if(subComp2.getGlobalCovariance().elementSum() > CONST_SMALL_TOLERANCE)
						splitDist2 = subComp2.split(tmpWeight);
					else{ 
						subComp2.scaleGlobalWeight(tmpWeight);
						splitDist2 = subComp2;
					}
					dist.getSubDistributions().set(i,splitDist1);
					dist.getSubDistributions().add(splitDist2);
				}
			}
		}
	}
	
	private static double mergeTwoClosestComps(SampleModel dist) throws IllegalArgumentException, InvocationTargetException, NoSuchMethodException,
			SecurityException, InstantiationException, IllegalAccessException {
		double compressionError = 0;
		TwoComponentDistribution twoCompDist = null;
		ArrayList<SimpleMatrix> means = dist.getSubMeans();
		ArrayList<SimpleMatrix> covs = dist.getSubCovariances();
		ArrayList<Double> weights = dist.getSubWeights();
		double distance = -1d;
		int indexComp1 = 0, indexComp2 = 0;
		for (int i = 0; i < means.size(); i++) {
			SimpleMatrix mean1 = means.get(i);
			for (int j = (i + 1); j < means.size(); j++) {
				SimpleMatrix mean2 = means.get(j);
				double tmpDistance = euclidianDistance(mean1, mean2);
				if ((distance == -1) | (tmpDistance < distance)) {
					distance = tmpDistance;
					indexComp1 = i;
					indexComp2 = j;
				}
				if (distance == 0)
					break;
			}
			if (distance == 0)
				break;
		}
		SimpleMatrix[] meansArray = { means.get(indexComp1), means.get(indexComp2) };
		SimpleMatrix[] covarianceArray = { covs.get(indexComp1), covs.get(indexComp2) };
		double[] weightsArray = { weights.get(indexComp1), weights.get(indexComp2) };
		try {
			twoCompDist = new TwoComponentDistribution(weightsArray, meansArray, covarianceArray, dist.getBandwidthMatrix());
			MomentMatcher.matchMoments(twoCompDist);
			OneComponentDistribution oneCompDist = new OneComponentDistribution(twoCompDist);
			compressionError = Hellinger.calculateUnscentedHellingerDistance(oneCompDist, twoCompDist);
		} catch (EmptyDistributionException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		Double[] subWeights = twoCompDist.getSubWeights();
		double newWeight1 = subWeights[0] / (subWeights[0] + subWeights[1]);
		double newWeight2 = subWeights[1] / (subWeights[0] + subWeights[1]);
		twoCompDist.getSubComponents()[0].setGlobalWeight(newWeight1);
		twoCompDist.getSubComponents()[1].setGlobalWeight(newWeight2);
		dist.getSubDistributions().set(indexComp2, twoCompDist);
		dist.getSubDistributions().remove(indexComp1);
		return compressionError;
	}

	public static double euclidianDistance(SimpleMatrix columnVector1, SimpleMatrix columnVector2) {
		double distance = 0;
		SimpleMatrix distVector = columnVector2.minus(columnVector1);
		distance = Math.sqrt(MatrixOps.elemPow(distVector, 2).elementSum());
		return distance;
	}
	
	/*private static double geoDistance(SimpleMatrix columnVector1, SimpleMatrix columnVector2){
		double lat1 = columnVector1.get(0,0);
		double lat2 = columnVector2.get(0,0);
		double lat11 = columnVector1.get(2,0);
		double lat21 = columnVector2.get(2,0);
		double lon1 = columnVector1.get(1,0);
		double lon2 = columnVector1.get(1,0);
		double lon11 = columnVector1.get(3,0);
		double lon21 = columnVector1.get(3,0);
		LatLng point1 = new LatLng(lat1, lon1);
		LatLng point2 = new LatLng(lat2, lon2);
		double distance1 = LatLngTool.distance(point1, point2, LengthUnit.METER);		
		LatLng point3 = new LatLng(lat11, lon11);
		LatLng point4 = new LatLng(lat21, lon21);
		double distance2 = LatLngTool.distance(point3, point4, LengthUnit.METER);		
		
		return (distance1+distance2);
	}*/

}
