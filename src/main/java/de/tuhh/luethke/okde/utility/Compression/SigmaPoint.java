/*
 * Copyright 2014 Jonas Luethke
 */


package de.tuhh.luethke.okde.utility.Compression;

import org.ejml.simple.SimpleMatrix;

public class SigmaPoint{
	private SimpleMatrix mPointVecor;
	// weight of the sigma point 
	private double mWeightInComponent;
	// weight of the sigma point multiplied with weight of corresponding component
	private double mWeight;
	public double getmWeightInComponent() {
		return mWeightInComponent;
	}
	public void setmWeightInComponent(double mWeightInComponent) {
		this.mWeightInComponent = mWeightInComponent;
	}
	public SigmaPoint(SimpleMatrix mPointVecor, double weight, double weightInComponent) {
		super();
		this.mPointVecor = mPointVecor;
		this.mWeight = weight;
		this.mWeightInComponent = weightInComponent;
	}
	public SimpleMatrix getmPointVecor() {
		return mPointVecor;
	}
	public void setmPointVecor(SimpleMatrix mPointVecor) {
		this.mPointVecor = mPointVecor;
	}
	public double getmWeight() {
		return mWeight;
	}
	public void setmWeight(double mWeight) {
		this.mWeight = mWeight;
	}
}
