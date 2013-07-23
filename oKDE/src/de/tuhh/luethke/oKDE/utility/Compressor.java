package de.tuhh.luethke.oKDE.utility;

import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;

import org.ejml.simple.SimpleMatrix;

import de.tuhh.luethke.oKDE.Exceptions.EmptyDistributionException;
import de.tuhh.luethke.oKDE.model.OneComponentDistribution;
import de.tuhh.luethke.oKDE.model.SampleModel;
import de.tuhh.luethke.oKDE.model.TwoComponentDistribution;

public class Compressor {
	private final static double D_TH = 0.02;

	public static SampleModel compress(SampleModel dist) throws IllegalArgumentException, InvocationTargetException, NoSuchMethodException, SecurityException, InstantiationException, IllegalAccessException {
		SampleModel inputModelCopy = new SampleModel(dist);
		double compressionError = mergeTwoClosestComps(dist);
		while (compressionError < D_TH) {
			inputModelCopy = new SampleModel(dist);
			compressionError = mergeTwoClosestComps(dist);
		}
		return inputModelCopy;
	}

	private static double mergeTwoClosestComps(SampleModel dist) throws IllegalArgumentException, InvocationTargetException, NoSuchMethodException, SecurityException, InstantiationException, IllegalAccessException {
		double compressionError = 0;
		TwoComponentDistribution twoCompDist = null;
		ArrayList<SimpleMatrix> means = dist.getSubMeans();
		ArrayList<SimpleMatrix> covs = dist.getSubCovariances();
		ArrayList<Double> weights = dist.getSubWeights();
		double distance = -1d;
		int indexComp1=0, indexComp2=0;
		for (int i = 0; i < means.size(); i++) {
			SimpleMatrix mean1 = means.get(i);
			for (int j = (i + 1); j < means.size(); j++) {
				SimpleMatrix mean2 = means.get(j);
				double tmpDistance = euclidianDistance(mean1,mean2);
				if( (distance == -1) | (tmpDistance < distance) ) {
					distance = tmpDistance;
					indexComp1 = i;
					indexComp2 = j;
				}
				if(distance == 0)
					break;
			}
			if(distance == 0)
				break;
		}
		SimpleMatrix[] meansArray = {means.get(indexComp1), means.get(indexComp2)};
		SimpleMatrix[] covarianceArray = {covs.get(indexComp1), covs.get(indexComp2)};
		double[] weightsArray = {weights.get(indexComp1), weights.get(indexComp2)};
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
		dist.getSubDistributions().set(indexComp2,twoCompDist);
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
