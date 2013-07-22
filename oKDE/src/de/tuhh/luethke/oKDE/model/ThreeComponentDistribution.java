package de.tuhh.luethke.oKDE.model;

import java.util.ArrayList;

import org.ejml.simple.SimpleMatrix;

import de.tuhh.luethke.oKDE.Exceptions.TooManyComponentsException;

public class ThreeComponentDistribution extends SampleDist{

	// component distributions
	private OneComponentDistribution[] mSubDistributions;

	public ThreeComponentDistribution(OneComponentDistribution[] components) throws TooManyComponentsException {
		// check number of components
		if (!(components.length == 3))
			throw new TooManyComponentsException();
		// add components to distribution
		mSubDistributions = components;
		mWeightSum = 0;
		for (OneComponentDistribution d : components) {
			mWeightSum += d.getWeightSum();
		}
		mN_eff = components.length;
		mForgettingFactor = 1;
	}
	
	public ThreeComponentDistribution(double[] weights, SimpleMatrix[] means, SimpleMatrix[] covariances) throws TooManyComponentsException {
		// check number of components
		if (!(weights.length == 3 & means.length == 3 & covariances.length == 3))
			throw new TooManyComponentsException();
		// add components to distribution
		mSubDistributions = new OneComponentDistribution[3];
		mSubDistributions[0] = new OneComponentDistribution(weights[0], means[0], covariances[0]);
		mSubDistributions[1] = new OneComponentDistribution(weights[1], means[1], covariances[1]);
		mSubDistributions[2] = new OneComponentDistribution(weights[2], means[2], covariances[2]);
		mWeightSum = 0;
		for (double w : weights) {
			mWeightSum += w;
		}
		mN_eff = weights.length;
		mForgettingFactor = 1;
	}

	public void setSubComponents(OneComponentDistribution[] subComponents) {
		this.mSubDistributions = subComponents;
	}

	public OneComponentDistribution[] getSubComponents() {
		return mSubDistributions;
	}
	
	public SimpleMatrix[] getSubMeans() {
		SimpleMatrix[] means = { mSubDistributions[0].getGlobalMean(), mSubDistributions[1].getGlobalMean(), mSubDistributions[2].getGlobalMean() };
		return means;
	}

	public SimpleMatrix[] getSubCovariances() {
		SimpleMatrix[] covs = { mSubDistributions[0].getGlobalCovariance(), mSubDistributions[1].getGlobalCovariance(), mSubDistributions[2].getGlobalCovariance() };
		return covs;
	}
	
	public Double[] getSubWeights() {
		Double[] weights = { mSubDistributions[0].getWeightSum(), mSubDistributions[1].getWeightSum(), mSubDistributions[2].getWeightSum() };
		return weights;
	}

	/*public void setSubMeans(SimpleMatrix mean1, SimpleMatrix mean2, ) {
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
