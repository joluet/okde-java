package de.tuhh.luethke.oKDE.utility;

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