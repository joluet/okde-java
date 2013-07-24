package de.tuhh.luethke.oKDE.utility;

import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;

import org.ejml.simple.SimpleMatrix;

import de.tuhh.luethke.oKDE.Exceptions.EmptyDistributionException;
import de.tuhh.luethke.oKDE.model.BaseSampleDistribution;
import de.tuhh.luethke.oKDE.model.OneComponentDistribution;
import de.tuhh.luethke.oKDE.model.SampleModel;
import de.tuhh.luethke.oKDE.model.TwoComponentDistribution;

public class Compressor {
	private final static double CONST_SMALL_TOLERANCE = 1E-10;
	
	private final static double D_TH = 0.02;

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

	public static void compress(SampleModel dist) throws Exception {
		// check wether compression is necessary using hysteresis rule
		if (dist.getSubMeans().size() <= dist.getNoOfCompsThreshold())
			return;
		ProjectionData projectionData = Projector.projectSampleDistToSubspace(dist);
		revitalizeComponents(dist);
		//System.out.println("COMPRESS");
		int noOfCompsBeforeCompression = dist.getSubMeans().size();
		SampleModel inputModelCopy = new SampleModel(dist);
		double compressionError = mergeTwoClosestComps(inputModelCopy);
		while (compressionError < D_TH) {
			dist.overWirite(inputModelCopy);
			compressionError = mergeTwoClosestComps(inputModelCopy);
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
				if (compressionError >= D_TH) {
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

	private static double euclidianDistance(SimpleMatrix columnVector1, SimpleMatrix columnVector2) {
		double distance = 0;
		SimpleMatrix distVector = columnVector2.minus(columnVector1);
		distance = MatrixOps.elemPow(distVector, 2).elementSum();
		return distance;
	}

}
