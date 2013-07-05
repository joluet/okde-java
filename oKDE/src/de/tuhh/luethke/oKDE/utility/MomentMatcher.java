package de.tuhh.luethke.oKDE.utility;

import java.util.ArrayList;

import org.ejml.simple.SimpleMatrix;

import de.tuhh.luethke.oKDE.model.SampleDist;
import de.tuhh.luethke.oKDE.Exceptions.EmptyDistributionException;

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
    public static void matchMoments(SampleDist distribution)
	    throws EmptyDistributionException {
	// Array of covariance matrices of components
	ArrayList<SimpleMatrix> smCovariances = distribution.getCovariances();

	// Array of mean vectors of components
	ArrayList<SimpleMatrix> smMeans = distribution.getMeans();

	// Array of component weights of components
	ArrayList<Double> smWeights = distribution.getWeights();

	// if the given distribution has only one component
	// just return empty covariance
	if (smWeights.size() == 1) {
	    SimpleMatrix newMean = smMeans.get(0);
	    SimpleMatrix newCovariance = null;
	    if (smCovariances.size() > 0)
		newCovariance = smCovariances.get(0);
	    distribution.setGlobalCovariance(newCovariance);
	    distribution.setGlobalMean(newMean);
	    distribution.setWeightSum(smWeights.get(0)); 
	    return;
	}

	// calculate new weight
	double newWeight = 0;
	for (double d : smWeights) {
	    newWeight += d;
	}
	// calculate new mean vector
	SimpleMatrix newMean = new SimpleMatrix(smMeans.get(0).numRows(),
		smMeans.get(0).numCols());
	for (int i = 0; i < smMeans.size(); i++) {
	    newMean = newMean.plus((smMeans.get(i).scale(smWeights.get(i))));
	}
	newMean = newMean.scale(1 / newWeight);
	// calculate new covariance matrix
	SimpleMatrix newCovariance = new SimpleMatrix(smCovariances.get(0)
		.numRows(), smCovariances.get(0).numCols());
	for (int i = 0; i < smCovariances.size(); i++) {
	    SimpleMatrix dyadSmMean = smMeans.get(i).mult(
		    smMeans.get(i).transpose());
	    SimpleMatrix S = smCovariances.get(i).plus(dyadSmMean);
	    newCovariance = newCovariance.plus(S.scale(smWeights.get(i)));
	}
	newCovariance = newCovariance.scale(1 / newWeight);
	SimpleMatrix dyadNewMean = newMean.mult(newMean.transpose());
	newCovariance = newCovariance.minus(dyadNewMean);

	// set calculated parameters to distribution
	distribution.setGlobalCovariance(newCovariance);
	distribution.setGlobalMean(newMean);
	distribution.setWeightSum(newWeight); 
    }
}
