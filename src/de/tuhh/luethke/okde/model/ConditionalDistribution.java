/*
 * Copyright 2014 Jonas Luethke
 */


package de.tuhh.luethke.okde.model;

import java.util.ArrayList;

import org.ejml.simple.SimpleMatrix;

public class ConditionalDistribution {

	public ArrayList<SimpleMatrix> conditionalMeans;

	public ArrayList<SimpleMatrix> conditionalCovs;
	
	public ArrayList<Double> conditionalWeights;

	public ConditionalDistribution(ArrayList<SimpleMatrix> conditionalMeans, ArrayList<SimpleMatrix> conditionalCovs, ArrayList<Double> conditionalWeights) {
		super();
		this.conditionalMeans = conditionalMeans;
		this.conditionalCovs = conditionalCovs;
		this.conditionalWeights = conditionalWeights;
	}
	
	

}
