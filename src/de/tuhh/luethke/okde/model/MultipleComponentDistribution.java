/*
 * Copyright 2014 Jonas Luethke
 */


package de.tuhh.luethke.okde.model;

import java.util.ArrayList;

import org.ejml.simple.SimpleMatrix;

import de.tuhh.luethke.okde.Exceptions.NoOfArgumentsException;

abstract public class MultipleComponentDistribution extends BaseSampleDistribution {

	// component distributions
	private OneComponentDistribution[] mSubDistributions;
	

	public MultipleComponentDistribution(double[] weights, SimpleMatrix[] means, SimpleMatrix[] covariances, SimpleMatrix bandwidth) {
		if(bandwidth == null)
			bandwidth = covariances[0].scale(0);
		mBandwidthMatrix = bandwidth;
		// add components to distribution
		mSubDistributions = new OneComponentDistribution[weights.length];
		for(int i=0; i<mSubDistributions.length; i++){
			mSubDistributions[i] = new OneComponentDistribution(weights[i], means[i], covariances[i], bandwidth);	
		}
		mGlobalWeight = 0;
		for (double w : weights) {
			mGlobalWeight += w;
		}
		mForgettingFactor = 1;
	}

	/**
	 * Copy constructor
	 * 
	 * @param dist
	 */
	public MultipleComponentDistribution(MultipleComponentDistribution dist) {
		OneComponentDistribution[] subDists = dist.getSubComponents();
		OneComponentDistribution[] copy = new OneComponentDistribution[subDists.length];
		for(int i=0; i<subDists.length; i++) {
			copy[i] = new OneComponentDistribution(subDists[i]);
		}
		this.mSubDistributions = copy;
		this.mBandwidthMatrix = dist.getBandwidthMatrix();
		this.mGlobalCovariance = dist.getGlobalCovariance();
		this.mGlobalMean = dist.getGlobalMean();
		this.mSubspaceGlobalCovariance = dist.getSubspaceGlobalCovariance();
		this.mSubspaceInverseCovariance = dist.getSubspaceInverseCovariance();
		this.mGlobalWeight = dist.getGlobalWeight();
	}
	

	@Override
	public double evaluate(SimpleMatrix pointVector) {
		SimpleMatrix[] means = this.getSubMeans();
		SimpleMatrix[] covs = this.getSubCovariances();
		Double[] weights = this.getSubWeights();
		double d = 0d;
		double n = means[0].numRows();
		double a = Math.pow(Math.sqrt(2 * Math.PI), n);
		for (int i = 0; i < means.length; i++) {
			SimpleMatrix m = means[i];
			SimpleMatrix c = covs[i].plus(this.mBandwidthMatrix);
			double w = weights[i];
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
	public ArrayList<Double> evaluate(ArrayList<SimpleMatrix> points) {
		ArrayList<Double> resultPoints = new ArrayList<Double>();
		for (SimpleMatrix point : points) {
			resultPoints.add(evaluate(point));
		}
		return resultPoints;
	}


	
	
	public void setSubComponents(OneComponentDistribution[] subComponents) {
		this.mSubDistributions = subComponents;
	}

	public OneComponentDistribution[] getSubComponents() {
		return mSubDistributions;
	}
	
	public SimpleMatrix[] getSubMeans() {
		SimpleMatrix[] means = new SimpleMatrix[mSubDistributions.length];
		for (int i=0; i<mSubDistributions.length; i++)
			means[i] = mSubDistributions[i].getGlobalMean();
		return means;
	}

	public SimpleMatrix[] getSubCovariances() {
		SimpleMatrix[] covs = new SimpleMatrix[mSubDistributions.length];
		for (int i=0; i<mSubDistributions.length; i++)
			covs[i] = mSubDistributions[i].getGlobalCovariance();
		return covs;
	}
	
	public Double[] getSubWeights() {
		Double[] weights = new Double[mSubDistributions.length];
		for (int i=0; i<mSubDistributions.length; i++)
			weights[i] = mSubDistributions[i].getGlobalWeight();
		return weights;
	}
	
	public void setSubMeans(SimpleMatrix[] means) throws NoOfArgumentsException {
		if(means.length != mSubDistributions.length)
			throw new NoOfArgumentsException();
		for (int i=0; i<mSubDistributions.length; i++)
			mSubDistributions[i].setGlobalMean(means[i]);
	}

	public void setSubCovariances(SimpleMatrix[] covariances) throws NoOfArgumentsException {
		if(covariances.length != mSubDistributions.length)
			throw new NoOfArgumentsException();
		for (int i=0; i<mSubDistributions.length; i++)
			mSubDistributions[i].setGlobalCovariance(covariances[i]);
	}
	
	@Override
	public void setBandwidthMatrix(SimpleMatrix mBandwidthMatrix) {
		this.mBandwidthMatrix = mBandwidthMatrix;
		for(BaseSampleDistribution d : mSubDistributions){
			d.setBandwidthMatrix(mBandwidthMatrix);
		}
	}
}
