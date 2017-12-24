/*
 * Copyright 2014 Jonas Luethke
 */


package de.tuhh.luethke.okde.utility;

import java.util.ArrayList;

import org.ejml.simple.SimpleMatrix;

import de.tuhh.luethke.okde.Exceptions.EmptyDistributionException;
import de.tuhh.luethke.okde.model.SampleModel;
import de.tuhh.luethke.okde.model.TwoComponentDistribution;

public class MomentMatcher {
	/**
	 * This method takes a sample distribution (mixture of Gaussians) with
	 * multiple mixture components and matches the first two moments (mean,
	 * covariance) of this distribution to generate a single component
	 * distribution.
	 * 
	 * The results are written to the the given distribution using
	 * setGlobalCovariance, setGlobalMean, setWeightSum.
	 * 
	 * @param distribution
	 * @throws EmptyDistributionException
	 */

	public static void matchMoments(SampleModel distribution) throws EmptyDistributionException {
		// Array of covariance matrices of components
		ArrayList<SimpleMatrix> smCovariances = distribution.getSubCovariances();
		// Array of mean vectors of components
		ArrayList<SimpleMatrix> smMeans = distribution.getSubMeans();

		// Array of component weights of components
		ArrayList<Double> smWeights = distribution.getSubWeights();

		// if the given distribution has only one component
		// just return empty covariance
		if (smWeights.size() == 0) {
			return;
		}
		if (smWeights.size() == 1) {
			SimpleMatrix newMean = smMeans.get(0);
			SimpleMatrix newCovariance = null;
			if (smCovariances.size() > 0)
				newCovariance = smCovariances.get(0);
			distribution.setGlobalCovariance(newCovariance);
			distribution.setGlobalMean(newMean);
			distribution.setGlobalWeight(smWeights.get(0));
			return;
		}

		// calculate new weight
		double newWeight = 0;
		for (double d : smWeights) {
			newWeight += d;
		}
		// calculate new mean vector
		SimpleMatrix newMean = new SimpleMatrix(smMeans.get(0).numRows(), smMeans.get(0).numCols());
		for (int i = 0; i < smMeans.size(); i++) {
			newMean = newMean.plus((smMeans.get(i).scale(smWeights.get(i))));
		}
		newMean = newMean.scale(1 / newWeight);
		// calculate new covariance matrix
		SimpleMatrix newCovariance = new SimpleMatrix(smCovariances.get(0).numRows(), smCovariances.get(0).numCols());
		for (int i = 0; i < smCovariances.size(); i++) {
			SimpleMatrix dyadSmMean = smMeans.get(i).mult(smMeans.get(i).transpose());
			SimpleMatrix S = smCovariances.get(i).plus(dyadSmMean);
			newCovariance = newCovariance.plus(S.scale(smWeights.get(i)));
		}
		newCovariance = newCovariance.scale(1 / newWeight);
		SimpleMatrix dyadNewMean = newMean.mult(newMean.transpose());
		newCovariance = newCovariance.minus(dyadNewMean);
		//System.out.println("matching moments");
		//System.out.println(newCovariance);
		//System.out.println(newMean);
		// set calculated parameters to distribution
		distribution.setGlobalCovariance(newCovariance);
		distribution.setGlobalMean(newMean);
		distribution.setGlobalWeight(newWeight);
	}
	
	public static void matchMoments(TwoComponentDistribution distribution) throws EmptyDistributionException {
		// Array of covariance matrices of components
		SimpleMatrix[] smCovariances = distribution.getSubCovariances();
		// Array of mean vectors of components
		SimpleMatrix[] smMeans = distribution.getSubMeans();

		// Array of component weights of components
		Double[] smWeights = distribution.getSubWeights();

		// if the given distribution has only one component
		// just return empty covariance
		if (smWeights.length == 0) {
			return;
		}
		if (smWeights.length == 1) {
			SimpleMatrix newMean = smMeans[0];
			SimpleMatrix newCovariance = null;
			if (smCovariances.length > 0)
				newCovariance = smCovariances[0];
			distribution.setGlobalCovariance(newCovariance);
			distribution.setGlobalMean(newMean);
			distribution.setGlobalWeight(smWeights[0]);
			return;
		}

		// calculate new weight
		double newWeight = 0;
		for (double d : smWeights) {
			newWeight += d;
		}
		// calculate new mean vector
		SimpleMatrix newMean = new SimpleMatrix(smMeans[0].numRows(), smMeans[0].numCols());
		for (int i = 0; i < smMeans.length; i++) {
			newMean = newMean.plus((smMeans[i].scale(smWeights[i])));
		}
		newMean = newMean.scale(1 / newWeight);
		// calculate new covariance matrix
		SimpleMatrix newCovariance = new SimpleMatrix(smCovariances[0].numRows(), smCovariances[0].numCols());
		for (int i = 0; i < smCovariances.length; i++) {
			SimpleMatrix dyadSmMean = smMeans[i].mult(smMeans[i].transpose());
			SimpleMatrix S = smCovariances[i].plus(dyadSmMean);
			newCovariance = newCovariance.plus(S.scale(smWeights[i]));
		}
		newCovariance = newCovariance.scale(1 / newWeight);
		SimpleMatrix dyadNewMean = newMean.mult(newMean.transpose());
		newCovariance = newCovariance.minus(dyadNewMean);
		//System.out.println("matching moments");
		//System.out.println(newCovariance);
		//System.out.println(newMean);
		// set calculated parameters to distribution
		distribution.setGlobalCovariance(newCovariance);
		distribution.setGlobalMean(newMean);
		distribution.setGlobalWeight(newWeight);
	}
}
