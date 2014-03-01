/*
 * Copyright 2014 Jonas Luethke
 */


package de.tuhh.luethke.okde.model;

import java.util.ArrayList;

import org.ejml.simple.SimpleMatrix;
import org.ejml.simple.SimpleSVD;

import de.tuhh.luethke.okde.Exceptions.TooManyComponentsException;
import de.tuhh.luethke.okde.utility.Matrices.MatrixOps;

public class OneComponentDistribution extends BaseSampleDistribution {

	public OneComponentDistribution(double w, SimpleMatrix mean, SimpleMatrix covariance, SimpleMatrix bandwidth) {
		super();
		mGlobalWeight = w;
		mGlobalMean = mean;
		mGlobalCovariance = covariance;
		mBandwidthMatrix = bandwidth;
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
	 * Splits a single component distribution into two components as described in the oKDE-paper.
	 * @return a TwoComponentDistribution
	 */
	public TwoComponentDistribution split(double parentWeight){
		SimpleSVD<?> svd = mGlobalCovariance.svd(true);
		SimpleMatrix S = svd.getW();
		SimpleMatrix V = svd.getV();
		SimpleMatrix d = S.extractDiag();
		double max = MatrixOps.maxVectorElement(d);
		int maxIndex = MatrixOps.maxVectorElementIndex(d);
		int len = mGlobalCovariance.numRows();
		SimpleMatrix M = new SimpleMatrix(len,1);
		M.set(maxIndex, 0, 1.0d);
		SimpleMatrix dMean = V.mult(M).scale(0.5*Math.sqrt(max));
		SimpleMatrix meanSplit1 = mGlobalMean.plus(dMean);
		SimpleMatrix meanSplit2 = mGlobalMean.minus(dMean);
		
		SimpleMatrix dyadMean = mGlobalMean.mult(mGlobalMean.transpose());
		SimpleMatrix dyadMeanSplit1 = meanSplit1.mult(meanSplit1.transpose());
		SimpleMatrix dyadMeanSplit2 = meanSplit2.mult(meanSplit2.transpose());
		SimpleMatrix covSplit = mGlobalCovariance.plus(dyadMean).minus(dyadMeanSplit1.plus(dyadMeanSplit2).scale(0.5));
		
		SimpleMatrix[] means = {meanSplit1, meanSplit2};
		SimpleMatrix[] covariances = {covSplit, covSplit};
		double[] weights = {0.5, 0.5};
		TwoComponentDistribution splitDist = null;
		try {
			splitDist = new TwoComponentDistribution(weights, means, covariances, mBandwidthMatrix);
			splitDist.setGlobalWeight(parentWeight*mGlobalWeight);
			splitDist.setGlobalCovariance(mGlobalCovariance);
			splitDist.setGlobalMean(mGlobalMean);
		} catch (TooManyComponentsException e) {
			// cant be thrown
		}
		return splitDist;
	}

	/**
	 * @see de.tuhh.luethke.okde.model.BaseSampleDistribution#evaluate(SimpleMatrix pointVector)
	 */
	@Override
	public double evaluate(SimpleMatrix pointVector) {
		SimpleMatrix smoothedCov = mGlobalCovariance.plus(mBandwidthMatrix);
		double d = 0d;
		double n = mGlobalMean.numRows();
		double a = Math.pow(Math.sqrt(2 * Math.PI), n);
		double tmp = (-0.5d) * pointVector.minus(mGlobalMean).transpose().mult(smoothedCov.invert()).mult(pointVector.minus(mGlobalMean)).trace();
		d += ((1 / (a * Math.sqrt(smoothedCov.determinant()))) * Math.exp(tmp)) * mGlobalWeight;

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
	
	/**
	 * @see de.tuhh.luethke.okde.model.BaseSampleDistribution#setBandwidthMatrix(SimpleMatrix mBandwidthMatrix)
	 */
	@Override
	public void setBandwidthMatrix(SimpleMatrix mBandwidthMatrix) {
		this.mBandwidthMatrix = mBandwidthMatrix;
	}
}
