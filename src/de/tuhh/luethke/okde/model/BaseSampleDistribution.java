/*
 * Copyright 2014 Jonas Luethke
 */


package de.tuhh.luethke.okde.model;

import java.util.ArrayList;

import org.ejml.simple.SimpleMatrix;

/**
 * This is the abstract base class for sample distributions.
 * 
 * @author Jonas Luethke
 *
 */
public abstract class BaseSampleDistribution {

	// bandwidth matrix
	protected SimpleMatrix mBandwidthMatrix;

	// overall weight sum
	protected double mGlobalWeight;

	// overall mean
	protected SimpleMatrix mGlobalMean;

	// overall covariance
	protected SimpleMatrix mGlobalCovariance;

	// overall covariance in subspace
	protected SimpleMatrix mSubspaceGlobalCovariance;

	// overall inverse covariance in subspace
	protected SimpleMatrix mSubspaceInverseCovariance;

	// forgetting factor, to be used for non-stationary distributions
	protected double mForgettingFactor;

	public double getForgettingFactor() {
		return mForgettingFactor;
	}

	public void setForgettingFactor(double forgettingFactor) {
		this.mForgettingFactor = forgettingFactor;
	}

	/**
	 * Create new SampleDist
	 */
	public BaseSampleDistribution() {
		super();
		mGlobalWeight = 0;
		mForgettingFactor = 0;
	}

	public double getGlobalWeight() {
		return mGlobalWeight;
	}

	public SimpleMatrix getBandwidthMatrix() {
		return mBandwidthMatrix;
	}

	abstract public void setBandwidthMatrix(SimpleMatrix mBandwidthMatrix);

	public SimpleMatrix getGlobalCovariance() {
		return mGlobalCovariance;
	}

	public void setGlobalCovariance(SimpleMatrix globalCovariance) {
		this.mGlobalCovariance = globalCovariance;
	}

	public SimpleMatrix getSubspaceGlobalCovariance() {
		return mSubspaceGlobalCovariance;
	}

	public void setGlobalWeight(double weight) {
		this.mGlobalWeight = weight;
	}
	
	public void scaleGlobalWeight(double scaleFactor) {
		this.mGlobalWeight = this.mGlobalWeight*scaleFactor;
	}

	public void setSubspaceGlobalCovariance(SimpleMatrix subspaceCovariance) {
		this.mSubspaceGlobalCovariance = subspaceCovariance;
	}

	public SimpleMatrix getSubspaceInverseCovariance() {
		return mSubspaceInverseCovariance;
	}

	public void setSubspaceInverseCovariance(SimpleMatrix subspaceInverseCovariance) {
		this.mSubspaceInverseCovariance = subspaceInverseCovariance;
	}

	/**
	 * Evaluates the distribution at the given n-dimensional point and returns
	 * the result as a double-value.
	 * 
	 * @param point point to evaluate
	 * @return result as double value
	 */
	abstract public double evaluate(SimpleMatrix pointVector);

	/**
	 * Evaluates the distribution at the given n-dimensional points and returns
	 * the results in a List of double-values.
	 * 
	 * @param points points to evaluate
	 * @return result as array of double values
	 */
	abstract public ArrayList<Double> evaluate(ArrayList<SimpleMatrix> points);

	public SimpleMatrix getGlobalMean() {
		return mGlobalMean;
	}

	public void setGlobalMean(SimpleMatrix globalMean) {
		this.mGlobalMean = globalMean;
	}

	public SimpleMatrix getmGlobalCovarianceSmoothed() {
		if (mBandwidthMatrix == null)
			mBandwidthMatrix = mGlobalCovariance.scale(0);
		return (mGlobalCovariance.plus(mBandwidthMatrix));
	}

}
