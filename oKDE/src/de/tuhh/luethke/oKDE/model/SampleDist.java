package de.tuhh.luethke.oKDE.model;

import java.util.ArrayList;

import org.ejml.simple.SimpleMatrix;

public abstract class SampleDist {

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

	protected double mForgettingFactor;

	/**
	 * Create new SampleDist
	 */
	public SampleDist() {
		super();
		mGlobalWeight = 0;
		mForgettingFactor = 1;
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

	public void setSubspaceGlobalCovariance(SimpleMatrix subspaceCovariance) {
		this.mSubspaceGlobalCovariance = subspaceCovariance;
	}

	public SimpleMatrix getSubspaceInverseCovariance() {
		return mSubspaceInverseCovariance;
	}

	public void setSubspaceInverseCovariance(SimpleMatrix subspaceInverseCovariance) {
		this.mSubspaceInverseCovariance = subspaceInverseCovariance;
	}

	abstract public double evaluate(SimpleMatrix pointVector);

	/**
	 * Evaluates the distribution at the given n-dimensional points and returns
	 * the results in a List of double-values.
	 * 
	 * @param points
	 * @return array of double values
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
