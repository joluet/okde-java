package de.tuhh.luethke.oKDE.model;

import java.util.ArrayList;
import java.util.List;

import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;
import org.ejml.simple.SimpleMatrix;
import org.ejml.simple.SimpleSVD;

import de.tuhh.luethke.oKDE.Exceptions.EmptyDistributionException;
import de.tuhh.luethke.oKDE.utility.MomentMatcher;

public abstract class SampleDist {
	
	//TODO: remove this field
	public double compressionError;

	// bandwidth factor
	protected double mBandwidthFactor;

	// bandwidth matrix
	protected SimpleMatrix mBandwidthMatrix;

	// overall weight sum
	protected double mWeightSum;

	// overall mean
	protected SimpleMatrix mGlobalMean;

	// overall covariance
	protected SimpleMatrix mGlobalCovariance;

	// overall covariance plus bandwidth
	protected SimpleMatrix mGlobalCovarianceSmoothed;

	// overall covariance in subspace
	protected SimpleMatrix mSubspaceGlobalCovariance;

	// overall inverse covariance in subspace
	protected SimpleMatrix mSubspaceInverseCovariance;



	protected double mForgettingFactor;

	// effective number of observed samples
	protected double mN_eff;

	/**
	 * Create new SampleDist
	 */
	/*public SampleDist() {
		super();
		mWeightSum = 0;
		mN_eff = 0;
		mForgettingFactor = 1;
		mSubDistributions = new ArrayList<SampleDist>();
	}*/

	/*public SampleDist(double[] weights, SimpleMatrix[] means, SimpleMatrix[] covariances) {
		super();
		mSubDistributions = new ArrayList<SampleDist>();
		for(int i=0; i< means.length; i++)
			mSubDistributions.add(new SampleDist(weights[i], means[i], covariances[i]));
		mWeightSum = 0;
		for (double w : weights) {
			mWeightSum += w;
		}
		mN_eff = weights.length;
		mForgettingFactor = 1;
	}

	public SampleDist(double w, SimpleMatrix mean, SimpleMatrix covariance) {
		super();
		mWeightSum = w;
		mGlobalMean = mean;
		mGlobalCovariance = covariance;
		// initialize smoothed covariance to zero
		mGlobalCovarianceSmoothed = new SimpleMatrix(covariance).scale(0);
		// initialize bandwidth matrix to zero
		mBandwidthMatrix = new SimpleMatrix(covariance).scale(0);
		mN_eff = 0;
		mForgettingFactor = 1;
		mSubDistributions = new ArrayList<SampleDist>();
	}*/

	/**
	 * Copy constructor
	 * 
	 * @param dist
	 */
	/* SampleDist(SampleDist dist) {
		List<SampleDist> subDists = dist.getSubDistributions();
		ArrayList<SampleDist> copy = new ArrayList<SampleDist>();
		for (SampleDist d : subDists) {
			copy.add(new SampleDist(d));
		}
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
	}*/

	public double getWeightSum() {
		return mWeightSum;
	}

	public double getmBandwidthFactor() {
		return mBandwidthFactor;
	}

	public void setmBandwidthFactor(double mBandwidthFactor) {
		this.mBandwidthFactor = mBandwidthFactor;
	}

	public SimpleMatrix getmBandwidthMatrix() {
		return mBandwidthMatrix;
	}

	public void setmBandwidthMatrix(SimpleMatrix mBandwidthMatrix) {
		if (this.mGlobalCovariance == null)
			this.mGlobalCovarianceSmoothed = new SimpleMatrix(mBandwidthMatrix);
		else
			this.mGlobalCovarianceSmoothed = this.mGlobalCovariance.plus(mBandwidthMatrix);
		this.mBandwidthMatrix = mBandwidthMatrix;
	}

	public SimpleMatrix getGlobalCovariance() {
		return mGlobalCovariance;
	}

	public void setGlobalCovariance(SimpleMatrix globalCovariance) {
		this.mGlobalCovariance = globalCovariance;
	}

	public SimpleMatrix getSubspaceGlobalCovariance() {
		return mSubspaceGlobalCovariance;
	}

	public void setWeightSum(double weightSum) {
		this.mWeightSum = weightSum;
	}

	public void setSubspaceGlobalCovariance(SimpleMatrix mSubspaceCovariance) {
		this.mSubspaceGlobalCovariance = mSubspaceCovariance;
	}

	public SimpleMatrix getmSubspaceInverseCovariance() {
		return mSubspaceInverseCovariance;
	}

	public void setmSubspaceInverseCovariance(SimpleMatrix mSubspaceInverseCovariance) {
		this.mSubspaceInverseCovariance = mSubspaceInverseCovariance;
	}


	abstract public double evaluate(SimpleMatrix pointVector, boolean useSubDists, boolean useSmoothedCov);
	
	/**
	 * Evaluates the distribution at the given n-dimensional points and returns the results in a List of double-values.
	 * @param points 
	 * @return array of double values
	 */
	abstract public ArrayList<Double> evaluate(ArrayList<SimpleMatrix> points, boolean useSubDists, boolean useSmoothedCov);


	public SimpleMatrix getGlobalMean() {
		return mGlobalMean;
	}

	public void setGlobalMean(SimpleMatrix globalMean) {
		this.mGlobalMean = globalMean;
	}

	public double getNeff() {
		return this.mN_eff;
	}

	public void setGlobalMean(double nEff) {
		this.mN_eff = nEff;
	}

	public SimpleMatrix getmGlobalCovarianceSmoothed() {
		return mGlobalCovarianceSmoothed;
	}

	public void setmGlobalCovarianceSmoothed(SimpleMatrix mGlobalCovarianceSmoothed) {
		this.mGlobalCovarianceSmoothed = mGlobalCovarianceSmoothed;
	}

}
