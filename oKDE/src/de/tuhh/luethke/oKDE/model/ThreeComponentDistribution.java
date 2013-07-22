package de.tuhh.luethke.oKDE.model;

import org.ejml.simple.SimpleMatrix;

import de.tuhh.luethke.oKDE.Exceptions.TooManyComponentsException;

public class ThreeComponentDistribution extends MultipleComponentDistribution{

	public ThreeComponentDistribution(double[] weights, SimpleMatrix[] means, SimpleMatrix[] covariances) throws TooManyComponentsException {
		super(weights, means, covariances);
		// check number of components
		if (!(weights.length == 3 & means.length == 3 & covariances.length == 3))
			throw new TooManyComponentsException();
	}
	
}
