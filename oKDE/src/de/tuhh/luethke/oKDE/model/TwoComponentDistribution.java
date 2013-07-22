package de.tuhh.luethke.oKDE.model;

import java.util.ArrayList;

import org.ejml.simple.SimpleMatrix;

import de.tuhh.luethke.oKDE.Exceptions.TooManyComponentsException;

public class TwoComponentDistribution extends SampleDist {

	// component distributions
	private OneComponentDistribution[] mSubDistributions;

	public TwoComponentDistribution(double[] weights, SimpleMatrix[] means, SimpleMatrix[] covariances) throws TooManyComponentsException {
		// check number of components
		if (!(weights.length == 2 & means.length == 2 & covariances.length == 2))
			throw new TooManyComponentsException();
		// add components to distribution
		mSubDistributions = new OneComponentDistribution[2];
		mSubDistributions[0] = new OneComponentDistribution(weights[0], means[0], covariances[0]);
		mSubDistributions[1] = new OneComponentDistribution(weights[1], means[1], covariances[1]);
		mWeightSum = 0;
		for (double w : weights) {
			mWeightSum += w;
		}
		mN_eff = weights.length;
		mForgettingFactor = 1;
	}

	public TwoComponentDistribution(double w, SimpleMatrix mean, SimpleMatrix covariance) {
		mWeightSum = w;
		mGlobalMean = mean;
		mGlobalCovariance = covariance;
		// initialize smoothed covariance to zero
		mGlobalCovarianceSmoothed = new SimpleMatrix(covariance).scale(0);
		// initialize bandwidth matrix to zero
		mBandwidthMatrix = new SimpleMatrix(covariance).scale(0);
		mN_eff = 0;
		mForgettingFactor = 1;
		mSubDistributions = null;
	}

	/**
	 * Copy constructor
	 * 
	 * @param dist
	 */
	TwoComponentDistribution(TwoComponentDistribution dist) {
		OneComponentDistribution[] subDists = dist.getSubComponents();
		OneComponentDistribution[] copy = {new OneComponentDistribution(subDists[0]),
				new OneComponentDistribution(subDists[1])};
		this.mSubDistributions = copy;
		this.mGlobalCovarianceSmoothed = dist.getmGlobalCovarianceSmoothed();
		this.mBandwidthFactor = dist.getmBandwidthFactor();
		this.mBandwidthMatrix = dist.getmBandwidthMatrix();
		this.mGlobalCovariance = dist.getGlobalCovariance();
		this.mGlobalMean = dist.getGlobalMean();
		this.mSubspaceGlobalCovariance = dist.getSubspaceGlobalCovariance();
		this.mSubspaceInverseCovariance = dist.getmSubspaceInverseCovariance();
		this.mWeightSum = dist.getWeightSum();
		this.mN_eff = dist.getNeff();
	}

	public void setSubComponents(OneComponentDistribution[] subComponents) {
		this.mSubDistributions = subComponents;
	}

	public OneComponentDistribution[] getSubComponents() {
		return mSubDistributions;
	}

	public void setSubMeans(SimpleMatrix mean1, SimpleMatrix mean2) {
		mSubDistributions[0].setGlobalMean(mean1);
		mSubDistributions[1].setGlobalMean(mean2);
	}

	public void setSubCovariances(SimpleMatrix covariance1, SimpleMatrix covariance2) {
		mSubDistributions[0].setGlobalCovariance(covariance1);
		mSubDistributions[1].setGlobalCovariance(covariance1);
	}

	public SimpleMatrix[] getSubSmoothedCovariances() {
		SimpleMatrix[] covs = { mSubDistributions[0].getmGlobalCovarianceSmoothed(), mSubDistributions[1].getmGlobalCovarianceSmoothed() };
		return covs;
	}

	public SimpleMatrix[] getSubMeans() {
		SimpleMatrix[] means = { mSubDistributions[0].getGlobalMean(), mSubDistributions[1].getGlobalMean() };
		return means;
	}

	public SimpleMatrix[] getSubCovariances() {
		SimpleMatrix[] covs = { mSubDistributions[0].getGlobalCovariance(), mSubDistributions[1].getGlobalCovariance() };
		return covs;
	}

	public Double[] getSubWeights() {
		Double[] weights = { mSubDistributions[0].getWeightSum(), mSubDistributions[1].getWeightSum() };
		return weights;
	}

	public void setSubWeights(double weight1, double weight2) {
		mSubDistributions[0].setWeightSum(weight1);
		mSubDistributions[1].setWeightSum(weight2);
	}

	@Override
	public double evaluate(SimpleMatrix pointVector, boolean useSubDists, boolean useSmoothedCov) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public ArrayList<Double> evaluate(ArrayList<SimpleMatrix> points, boolean useSubDists, boolean useSmoothedCov) {
		// TODO Auto-generated method stub
		return null;
	}

}
