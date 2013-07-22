package de.tuhh.luethke.oKDE.model;

import java.util.ArrayList;

import org.ejml.simple.SimpleMatrix;

public class OneComponentDistribution extends SampleDist {

	public OneComponentDistribution(double w, SimpleMatrix mean, SimpleMatrix covariance) {
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
	}

	public OneComponentDistribution(OneComponentDistribution oneComponentDistribution) {
		this.mGlobalCovarianceSmoothed = oneComponentDistribution.getmGlobalCovarianceSmoothed();
		this.mBandwidthFactor = oneComponentDistribution.getmBandwidthFactor();
		this.mBandwidthMatrix = oneComponentDistribution.getmBandwidthMatrix();
		this.mGlobalCovariance = oneComponentDistribution.getGlobalCovariance();
		this.mGlobalMean = oneComponentDistribution.getGlobalMean();
		this.mSubspaceGlobalCovariance = oneComponentDistribution.getSubspaceGlobalCovariance();
		this.mSubspaceInverseCovariance = oneComponentDistribution.getmSubspaceInverseCovariance();
		this.mWeightSum = oneComponentDistribution.getWeightSum();
		this.mN_eff = oneComponentDistribution.getNeff();
	}
	
	public OneComponentDistribution(TwoComponentDistribution twoComponentDistribution) {
		this.mGlobalCovarianceSmoothed = twoComponentDistribution.getmGlobalCovarianceSmoothed();
		this.mBandwidthFactor = twoComponentDistribution.getmBandwidthFactor();
		this.mBandwidthMatrix = twoComponentDistribution.getmBandwidthMatrix();
		this.mGlobalCovariance = twoComponentDistribution.getGlobalCovariance();
		this.mGlobalMean = twoComponentDistribution.getGlobalMean();
		this.mSubspaceGlobalCovariance = twoComponentDistribution.getSubspaceGlobalCovariance();
		this.mSubspaceInverseCovariance = twoComponentDistribution.getmSubspaceInverseCovariance();
		this.mWeightSum = twoComponentDistribution.getWeightSum();
		this.mN_eff = twoComponentDistribution.getNeff();
	}

	@Override
	public double evaluate(SimpleMatrix pointVector, boolean useSubDists, boolean useSmoothedCov) {
		SimpleMatrix bandwidth = this.mBandwidthMatrix;
		/*
		 * double[][] dxVector = { { x }, { y } }; SimpleMatrix xVector = new
		 * SimpleMatrix(dxVector);
		 */
		double d = 0d;
		double n = mGlobalMean.numRows();
		double a = Math.pow(Math.sqrt(2 * Math.PI), n);
		double tmp = (-0.5d) * pointVector.minus(mGlobalMean).transpose().mult(mGlobalCovariance.invert()).mult(pointVector.minus(mGlobalMean)).trace();
		d += ((1 / (a * Math.sqrt(mGlobalCovariance.determinant()))) * Math.exp(tmp)) * mWeightSum;

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
}
