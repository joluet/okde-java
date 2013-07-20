package de.tuhh.luethke.oKDE.utility;

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

	public static SimpleMatrix deleteElementsFromVector(SimpleMatrix vector, List<Double> elements) {
		SimpleMatrix newVector = new SimpleMatrix(elements.size(), 1);
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
				matrix.set(i,j,1);
			}
		}
		return matrix;
	}
}
