package de.tuhh.luethke.oKDE.utility;

import org.ejml.simple.SimpleMatrix;

public class MatrixOps {

    public static SimpleMatrix abs(SimpleMatrix matrix) {
	for(int i=0; i<matrix.numRows(); i++){
	    for(int j=0; j<matrix.numCols(); j++){
		matrix.set(i,j,Math.abs(matrix.get(i, j)));
	    }
	}
	return matrix;
    }
    
    public static SimpleMatrix elemSqrt(SimpleMatrix matrix) {
	for(int i=0; i<matrix.numRows(); i++){
	    for(int j=0; j<matrix.numCols(); j++){
		matrix.set(i,j,Math.sqrt(matrix.get(i, j)));
	    }
	}
	return matrix;
    }
    
    public static SimpleMatrix elemPow(SimpleMatrix matrix, double p) {
	for(int i=0; i<matrix.numRows(); i++){
	    for(int j=0; j<matrix.numCols(); j++){
		matrix.set(i,j,Math.pow(matrix.get(i, j),p));
	    }
	}
	return matrix;
    }
}
