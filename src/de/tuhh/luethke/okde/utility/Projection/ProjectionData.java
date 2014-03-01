/*
 * Copyright 2014 Jonas Luethke
 */


package de.tuhh.luethke.okde.utility.Projection;

import javax.swing.text.StyleContext.SmallAttributeSet;

import org.ejml.simple.SimpleMatrix;
import org.ejml.simple.SimpleSVD;

public class ProjectionData {

	public SimpleSVD<?> mSVD;
	public Double[] mValidElements;
	public int mCountValidElements;
	public SimpleMatrix mBandwidthMatrix;
	public SimpleMatrix mGlobalMean;
}
