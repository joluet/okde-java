package de.tuhh.luethke.oKDE.utility;

import org.ejml.simple.SimpleMatrix;

import de.tuhh.luethke.oKDE.model.SampleDist;


public class MomentMatcher {
    public SampleDist matchMoments(SampleDist subMixture) {
	SimpleMatrix[] smCovariances = subMixture.getCovariance();
	SimpleMatrix[] smSampleData = subMixture.getSampleData();
	double[] smWeights = subMixture.getWeights();

	// calculate new weight
	double weightNew = 0;
	for (double d : smWeights) {
	    weightNew += d;
	}

	// calculate new sample data
	SimpleMatrix sampleDataNew = new SimpleMatrix(smSampleData[0].numRows(), smSampleData[0].numCols());
	for (int i = 0; i < smSampleData.length; i++) {
	    sampleDataNew = sampleDataNew.plus((smSampleData[i]
		    .scale(smWeights[i])));
	}
	sampleDataNew.scale(1 / weightNew);

	// calculate new covariance
	SimpleMatrix covarianceNew = new SimpleMatrix(smCovariances[0].numRows(), smCovariances[0].numCols());
	for (int i = 0; i < smCovariances.length; i++) {
	    
	    SimpleMatrix dyadSmSample = smSampleData[i].mult(smSampleData[i].transpose());		
	    	
	    SimpleMatrix S = smCovariances[i].plus(dyadSmSample);
	    
	    covarianceNew = S.scale(smWeights[i]);
	}
	covarianceNew = covarianceNew.scale(1 / weightNew);
	SimpleMatrix dyadNewSample = sampleDataNew.mult(sampleDataNew.transpose());
	covarianceNew = covarianceNew.minus(dyadNewSample);
	
	// put parameters into sample dist
	SampleDist compressed = new SampleDist(weightNew, sampleDataNew, covarianceNew);
	
	return compressed;
    }
}
