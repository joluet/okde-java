/*
 * Copyright 2014 Jonas Luethke
 */


package de.tuhh.luethke.okde.utility.Optimization;

import org.ejml.simple.SimpleMatrix;

public class SearchResult {

	public SimpleMatrix point;
	public double probability;
	public SearchResult(SimpleMatrix point, double probability) {
		super();
		this.point = point;
		this.probability = probability;
	}
	
	

}
