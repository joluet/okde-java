package de.tuhh.luethke.oKDE.utility;

import java.util.ArrayList;
import java.util.List;

import org.ejml.simple.SimpleMatrix;
import org.ejml.simple.SimpleSVD;

import de.tuhh.luethke.oKDE.Exceptions.EmptyDistributionException;
import de.tuhh.luethke.oKDE.model.SampleDist;

public class Hellinger {

	private static final double MIN_TOL = 1e-5;

	public static void calculateUnscentedHellingerDistance(SampleDist dist1, SampleDist dist2) throws Exception {
		SampleDist dist0 = mergeSampleDists(dist1, dist2, 0.5, 0.5);
		double wSum=0;
		for(double w : dist0.getSubWeights()){
			wSum+=w;
		}
		System.out.println("WeightSum: "+wSum);
		List<SigmaPoint> sigmaPoints = getAllSigmaPoints(dist0, 3);
		System.out.println("sigmapoints: "+sigmaPoints.size());
		for(SigmaPoint p : sigmaPoints)
			System.out.println(p.getmPointVecor());
	}
	
	private static SampleDist mergeSampleDists(SampleDist dist1, SampleDist dist2, double w1, double w2) {
		SampleDist dist = new SampleDist();
		ArrayList<SimpleMatrix> means = dist2.getSubMeans();
		means.add(0,dist1.getGlobalMean());
		ArrayList<SimpleMatrix> covs = dist2.getSubSmoothedCovariances();
		covs.add(0,dist1.getGlobalCovariance());
		ArrayList<Double> weights = dist2.getSubWeights();
		weights.add(0,1d*w1);
		for(int i=1; i<weights.size(); i++){
			weights.set(i, weights.get(i)*w2);
		}
		double[] weightArray = new double[weights.size()];
		for(int i=0; i<weightArray.length; i++){
			System.out.println(weights.get(i));
			weightArray[i] = weights.get(i);
		}
		try {
			dist.updateDistribution(means.toArray(new SimpleMatrix[0]), covs.toArray(new SimpleMatrix[0]),
					weightArray);
			dist.setSubWeights(weights);
		} catch (EmptyDistributionException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return dist;
	}

	public static List<SigmaPoint> getAllSigmaPoints(SampleDist distribution, int max) throws Exception {
		ArrayList<SigmaPoint> sigmaPoints = new ArrayList<SigmaPoint>();
		int noOfComponents = distribution.getSubMeans().size();
		int dim = distribution.getSubMeans().get(0).numRows();
		int k = max - dim;
		int noOfSigmaPoints;
		if (k != 0)
			noOfSigmaPoints = 2 * dim + 1;
		else
			noOfSigmaPoints = 2 * dim;
		/*
		 * for(int i=0; i<(noOfSigmaPoints*noOfComponents); i++){ X.add(new
		 * SimpleMatrix(dim,0)); }
		 */
		
		ArrayList<Double> weights = new ArrayList<Double>();
		for (int i = 0; i < (2 * dim); i++) {
			weights.add(new Double(1d / (2 * ((double) dim + k))));
		}
		if (k != 0)
			weights.add(new Double((double) k / (double) (dim + k)));
		double sum = 0;
		for (Double d : weights) {
			sum += d;
		}
		if ((sum - 1) > MIN_TOL)
			throw new Exception("Weights in the unscented transform should sum to one!");

		for (int i = 0; i < noOfComponents; i++) {
			List<SimpleMatrix> x = getSigmaPoints(distribution.getSubMeans().get(i), distribution.getSubCovariances().get(i),
					noOfSigmaPoints, k);
			int count = 0;
			for(SimpleMatrix m : x){
				sigmaPoints.add(new SigmaPoint(m, weights.get(count)));
				count++;
			}
		}

		

		return sigmaPoints;
	}

	/**
	 * Returns 2n+k sigma points starting with mean as the first point
	 * 
	 * @param mean
	 * @param cov
	 * @param no
	 * @param k
	 * @return
	 */
	private static List<SimpleMatrix> getSigmaPoints(SimpleMatrix mean, SimpleMatrix cov, int no, int k) {
		List<SimpleMatrix> resultVectors = new ArrayList<SimpleMatrix>();

		int n = cov.numRows();
		SimpleSVD svd = cov.svd(true);
		SimpleMatrix U = svd.getU();
		SimpleMatrix S = svd.getW();
		SimpleMatrix V = svd.getV();

		S = U.mult(MatrixOps.elemSqrt(S)).scale(Math.sqrt(n + k));

		for (int i = 0; i < S.numCols(); i++) {
			SimpleMatrix columnVector = S.extractVector(false, i);
			SimpleMatrix negColumnVector = S.extractVector(false, i).scale(-1);
			resultVectors.add(columnVector.plus(mean));
			resultVectors.add(negColumnVector.plus(mean));
		}
		if (k != 0)
			resultVectors.add(mean);

		return resultVectors;
	}
	
}
