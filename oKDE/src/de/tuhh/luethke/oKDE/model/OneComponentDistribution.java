package de.tuhh.luethke.oKDE.model;

import java.util.ArrayList;

import org.ejml.simple.SimpleMatrix;

public class OneComponentDistribution extends BaseSampleDistribution {

	public OneComponentDistribution(double w, SimpleMatrix mean, SimpleMatrix covariance) {
		super();
		mGlobalWeight = w;
		mGlobalMean = mean;
		mGlobalCovariance = covariance;
		// initialize bandwidth matrix to zero
		mBandwidthMatrix = new SimpleMatrix(covariance).scale(0);
		mForgettingFactor = 1;
	}

	public OneComponentDistribution(OneComponentDistribution oneComponentDistribution) {
		this.mBandwidthMatrix = oneComponentDistribution.getBandwidthMatrix();
		this.mGlobalCovariance = oneComponentDistribution.getGlobalCovariance();
		this.mGlobalMean = oneComponentDistribution.getGlobalMean();
		this.mSubspaceGlobalCovariance = oneComponentDistribution.getSubspaceGlobalCovariance();
		this.mSubspaceInverseCovariance = oneComponentDistribution.getSubspaceInverseCovariance();
		this.mGlobalWeight = oneComponentDistribution.getGlobalWeight();
	}
	
	public OneComponentDistribution(TwoComponentDistribution twoComponentDistribution) {
		this.mBandwidthMatrix = twoComponentDistribution.getBandwidthMatrix();
		this.mGlobalCovariance = twoComponentDistribution.getGlobalCovariance();
		this.mGlobalMean = twoComponentDistribution.getGlobalMean();
		this.mSubspaceGlobalCovariance = twoComponentDistribution.getSubspaceGlobalCovariance();
		this.mSubspaceInverseCovariance = twoComponentDistribution.getSubspaceInverseCovariance();
		this.mGlobalWeight = twoComponentDistribution.getGlobalWeight();
	}

	/**
	 * @see de.tuhh.luethke.oKDE.model.BaseSampleDistribution#evaluate(SimpleMatrix pointVector)
	 */
	@Override
	public double evaluate(SimpleMatrix pointVector) {
		double d = 0d;
		double n = mGlobalMean.numRows();
		double a = Math.pow(Math.sqrt(2 * Math.PI), n);
		double tmp = (-0.5d) * pointVector.minus(mGlobalMean).transpose().mult(mGlobalCovariance.invert()).mult(pointVector.minus(mGlobalMean)).trace();
		d += ((1 / (a * Math.sqrt(mGlobalCovariance.determinant()))) * Math.exp(tmp)) * mGlobalWeight;

		return d;
	}

	/**
	 * @see de.tuhh.luethke.oKDE.model.BaseSampleDistribution#evaluate(ArrayList<SimpleMatrix> points)
	 */
	@Override
	public ArrayList<Double> evaluate(ArrayList<SimpleMatrix> points) {
		ArrayList<Double> resultPoints = new ArrayList<Double>();
		for (SimpleMatrix point : points) {
			resultPoints.add(evaluate(point));
		}
		return resultPoints;
	}
	
	@Override
	public void setBandwidthMatrix(SimpleMatrix mBandwidthMatrix) {
		this.mBandwidthMatrix = mBandwidthMatrix;
	}
}
