package de.tuhh.luethke.oKDE.model;

import java.util.ArrayList;

import org.ejml.simple.SimpleMatrix;

import de.tuhh.luethke.oKDE.Exceptions.TooManyComponentsException;

public class TwoComponentDistribution extends MultipleComponentDistribution {
	private final static int NO_OF_COMPONENTS = 2;


	public TwoComponentDistribution(double[] weights, SimpleMatrix[] means, SimpleMatrix[] covariances) throws TooManyComponentsException {
		super(weights, means, covariances);
		// check number of components
		if (!(weights.length == NO_OF_COMPONENTS & means.length == NO_OF_COMPONENTS & covariances.length == NO_OF_COMPONENTS))
			throw new TooManyComponentsException();
	}

	/*public TwoComponentDistribution(double w, SimpleMatrix mean, SimpleMatrix covariance) {
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
	}*/

	/**
	 * Copy constructor
	 * 
	 * @param dist
	 */
	TwoComponentDistribution(TwoComponentDistribution dist) {
		super(dist);
	}

/*	public void setSubComponents(OneComponentDistribution[] subComponents) {
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
	}*/



}
