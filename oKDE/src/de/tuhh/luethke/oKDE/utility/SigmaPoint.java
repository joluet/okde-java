package de.tuhh.luethke.oKDE.utility;

import org.ejml.simple.SimpleMatrix;

public class SigmaPoint{
	private SimpleMatrix mPointVecor;
	private double mWeight;
	public SigmaPoint(SimpleMatrix mPointVecor, double mWeight) {
		super();
		this.mPointVecor = mPointVecor;
		this.mWeight = mWeight;
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