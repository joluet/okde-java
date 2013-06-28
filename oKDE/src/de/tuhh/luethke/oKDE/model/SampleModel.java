package de.tuhh.luethke.oKDE.model;

import java.util.ArrayList;

import org.ejml.simple.SimpleMatrix;

public class SampleModel {

    private SampleDist mSampleDist;

    public SampleDist getSampleDist() {
	return mSampleDist;
    }

    public void setSampleDist(SampleDist sampleDist) {
	this.mSampleDist = sampleDist;
    }

    /**
     * Updates the given sample model with a single new observed data point or
     * several data points (batch mode)
     * 
     * @param model
     * @param forgettingFactor
     */
    private void updateSampleModel(float forgettingFactor,
	    SimpleMatrix observedData, double[] observedDataWeights) {
	ArrayList<Double> sampleWeights = mSampleDist.getWeights();
	double weightSum = mSampleDist.getWeightSum();
	

	mSampleDist.addWeights(observedDataWeights);
    }

}
