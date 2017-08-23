/*
 * Copyright 2014 Jonas Luethke
 */


package de.tuhh.luethke.okde.utility.Matrices;

import org.ejml.simple.SimpleMatrix;

public class HashableSimpleMatrix extends SimpleMatrix {

	public HashableSimpleMatrix(SimpleMatrix m) {
		super(m.getMatrix());
	}
	
	@Override
	public boolean equals(Object obj) {
		SimpleMatrix m = (SimpleMatrix) obj;
	    if(m.isIdentical(this, 1E-30))
	    	return true;
	    else
	    	return false;
	}
	
	@Override
    public int hashCode() {
		return (int) this.get(this.numRows()-1,0);
	}
}
