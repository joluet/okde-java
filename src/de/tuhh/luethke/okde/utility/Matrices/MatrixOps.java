/*
 * Copyright 2014 Jonas Luethke
 */


package de.tuhh.luethke.okde.utility.Matrices;

import java.util.List;

import org.ejml.simple.SimpleMatrix;

public class MatrixOps {

	public static SimpleMatrix abs(SimpleMatrix matrix) {
		for (int i = 0; i < matrix.numRows(); i++) {
			for (int j = 0; j < matrix.numCols(); j++) {
				matrix.set(i, j, Math.abs(matrix.get(i, j)));
			}
		}
		return matrix;
	}

	public static SimpleMatrix elemSqrt(SimpleMatrix matrix) {
		for (int i = 0; i < matrix.numRows(); i++) {
			for (int j = 0; j < matrix.numCols(); j++) {
				matrix.set(i, j, Math.sqrt(matrix.get(i, j)));
			}
		}
		return matrix;
	}

	public static SimpleMatrix elemPow(SimpleMatrix matrix, double p) {
		for (int i = 0; i < matrix.numRows(); i++) {
			for (int j = 0; j < matrix.numCols(); j++) {
				matrix.set(i, j, Math.pow(matrix.get(i, j), p));
			}
		}
		return matrix;
	}

	public static SimpleMatrix deleteElementsFromVector(SimpleMatrix vector, List<Double> elements, int vectorSize) {
		SimpleMatrix newVector = new SimpleMatrix(vectorSize, 1);
		int j = 0;
		for (int i = 0; i < vector.numRows(); i++)
			if (elements.get(i) == 1)
				newVector.set(j++, 0, vector.get(i));
		return newVector;
	}

	public static SimpleMatrix ones(int rows, int cols) {
		SimpleMatrix matrix = new SimpleMatrix(rows, cols);
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < rows; j++) {
				matrix.set(i, j, 1);
			}
		}
		return matrix;
	}

	public static SimpleMatrix doubleListToMatrix(List<Double> valueList) {
		SimpleMatrix m = new SimpleMatrix(1, valueList.size());
		for (int i = 0; i < valueList.size(); i++)
			m.set(0, i, valueList.get(i));
		return m;
	}

	public static List<Double> setNegativeValuesToZero(List<Double> valueList) {
		for (int i = 0; i < valueList.size(); i++) {
			if (valueList.get(i) < 0)
				valueList.set(i, 0d);
		}
		return valueList;
	}
	
	public static double maxVectorElement(SimpleMatrix matrix){
		double d = Double.MIN_VALUE;
		for (int i = 0; i < matrix.numRows(); i++) {
				if(matrix.get(i,0)>d)
					d = matrix.get(i,0);
		}
		return d;
	}
	public static int maxVectorElementIndex(SimpleMatrix matrix){
		double d = Double.MIN_VALUE;
		int row = 0;
		for (int i = 0; i < matrix.numRows(); i++) {
				if(matrix.get(i,0)>d){
					d = matrix.get(i,0);
					row = i;
				}
		}
		return row;
	}
}
