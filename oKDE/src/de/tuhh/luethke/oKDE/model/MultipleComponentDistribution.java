package de.tuhh.luethke.oKDE.model;

import java.util.ArrayList;

import org.ejml.simple.SimpleMatrix;

abstract class MultipleComponentDistribution extends SampleDist {
	
	public MultipleComponentDistribution(){
	}

	public MultipleComponentDistribution(double w, SimpleMatrix mean, SimpleMatrix covariance) {
		super(w, mean, covariance);
	}
	
	// component distributions
	private ArrayList<SampleDist> mSubDistributions;

	@Override
	public double evaluate(SimpleMatrix pointVector, boolean useSubDists, boolean useSmoothedCov) {
		ArrayList<SimpleMatrix> means = new ArrayList<SimpleMatrix>();
		ArrayList<SimpleMatrix> covs = new ArrayList<SimpleMatrix>();
		ArrayList<Double> weights = new ArrayList<Double>();
		if (useSubDists) {
			means = this.getSubMeans();
			if (useSmoothedCov)
				covs = this.getSubSmoothedCovariances();
			else
				covs = this.getSubCovariances();
			weights = this.getSubWeights();
		} else {
			means.add(this.mGlobalMean);
			covs.add(this.mGlobalCovariance);
			weights.add(this.mWeightSum);
		}

		SimpleMatrix bandwidth = this.mBandwidthMatrix;
		/*
		 * double[][] dxVector = { { x }, { y } }; SimpleMatrix xVector = new
		 * SimpleMatrix(dxVector);
		 */
		double d = 0d;
		double n = means.get(0).numRows();
		double a = Math.pow(Math.sqrt(2 * Math.PI), n);
		for (int i = 0; i < means.size(); i++) {
			SimpleMatrix m = means.get(i);
			SimpleMatrix c = covs.get(i);
			double w = weights.get(i);
			double tmp = (-0.5d) * pointVector.minus(m).transpose().mult(c.invert()).mult(pointVector.minus(m)).trace();
			d += ((1 / (a * Math.sqrt(c.determinant()))) * Math.exp(tmp)) * w;
		}
		return d;
	}

	/**
	 * Evaluates the distribution at the given n-dimensional points and returns
	 * the results in a List of double-values.
	 * 
	 * @param points
	 * @return array of double values
	 */
	@Override
	public ArrayList<Double> evaluate(ArrayList<SimpleMatrix> points, boolean useSubDists, boolean useSmoothedCov) {
		ArrayList<Double> resultPoints = new ArrayList<Double>();
		for (SimpleMatrix point : points) {
			resultPoints.add(evaluate(point, useSubDists, useSmoothedCov));
		}
		return resultPoints;
	}


	public void setSubMeans(ArrayList<SimpleMatrix> means) {
		for (int i = 0; i < mSubDistributions.size(); i++) {
			mSubDistributions.get(i).setGlobalMean(means.get(i));
		}
	}

	public void setSubCovariances(ArrayList<SimpleMatrix> covariances) {
		for (int i = 0; i < mSubDistributions.size(); i++) {
			mSubDistributions.get(i).setGlobalCovariance(covariances.get(i));
		}
	}

	public ArrayList<SimpleMatrix> getSubSmoothedCovariances() {
		ArrayList<SimpleMatrix> covs = new ArrayList<SimpleMatrix>();
		for (SampleDist d : mSubDistributions)
			covs.add(d.getmGlobalCovarianceSmoothed());
		return covs;
	}

	public ArrayList<SimpleMatrix> getSubMeans() {
		ArrayList<SimpleMatrix> means = new ArrayList<SimpleMatrix>();
		for (SampleDist d : mSubDistributions)
			means.add(d.getGlobalMean());
		return means;
	}

	public ArrayList<SimpleMatrix> getSubCovariances() {
		ArrayList<SimpleMatrix> covs = new ArrayList<SimpleMatrix>();
		for (SampleDist d : mSubDistributions)
			covs.add(d.getGlobalCovariance());
		return covs;
	}

	public ArrayList<Double> getSubWeights() {
		ArrayList<Double> weights = new ArrayList<Double>();
		for (SampleDist d : mSubDistributions)
			weights.add(d.getWeightSum());
		return weights;
	}

	public void setSubWeights(ArrayList<Double> weights) {
		for (int i = 0; i < mSubDistributions.size(); i++) {
			mSubDistributions.get(i).setWeightSum(weights.get(i));
		}
	}

}
