package de.tuhh.luethke.oKDE;

import java.util.ArrayList;

import org.ejml.simple.SimpleMatrix;

import de.tuhh.luethke.oKDE.model.SampleDist;
import de.tuhh.luethke.oKDE.model.SampleModel;

public class Algorithm {
    /**
     * Updates the given sample model with a single new observed data point or
     * several data points (batch mode)
     * 
     * @param model
     * @param forgettingFactor
     */
    private void updateSampleModel(SampleModel model, float forgettingFactor, SimpleMatrix observedData, double[] observedDataWeights) {
	SampleDist sampleDist = model.getSampleDist();
	ArrayList<Double> sampleWeights = sampleDist.getWeights();
	double weightSum = sampleDist.getWeightSum();
	
    }
}
