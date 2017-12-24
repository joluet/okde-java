/*
 * Copyright 2014 Jonas Luethke
 */


package de.tuhh.luethke.okde.utility.Optimization;

import java.util.ArrayList;
import java.util.List;

import org.ejml.simple.SimpleEVD;
import org.ejml.simple.SimpleMatrix;

import de.tuhh.luethke.okde.model.SampleModel;

/**
 * This class provides an optimization method for distributions obtained by kernel density estimation.
 * It contains a method for searching local maxima in a given distribution.
 * 
 * @author Jonas Luethke
 *
 */
public class Optimization {
	private static final double MAX_MAHALANOBIS_DIST = 25;
	
	private static final double START_STEP_SIZE = 5;
	
	private static final double STOP_STEP_SIZE = 1E-8;
	
	/**
	 * This method searches a local maximum by gradient-quadratic search. First a direct leap to the maximum by 
	 * quadratic optimization is tried. Then gradient search is used to refine the result in case of an overshoot.
	 * Uses means, covariances and component weights given as parameters.
	 * 
	 * This algorithm was motivated by this paper: 
	 * Miguel A. Carreira-Perpinan (2000): "Mode-finding for mixtures of
	 * Gaussian distributions", IEEE Trans. on Pattern Analysis and
	 * Machine Intelligence 22(11): 1318-1323.
	 * 
	 * @param start Defines the starting point for the search.
	 * @return The serach result containing the point and the probability value at that point.
	 */
	public static SearchResult gradQuadrSearch(SimpleMatrix start, ArrayList<SimpleMatrix> means, ArrayList<SimpleMatrix> covs, ArrayList<Double> weights, SampleModel model){
	
		SimpleMatrix gradient = new SimpleMatrix(2,1);
		SimpleMatrix hessian = new SimpleMatrix(2,2);
		double n = means.get(0).numRows();
		double a = Math.pow(Math.sqrt(2 * Math.PI), n);
		
		SimpleMatrix x = new SimpleMatrix(2,1);
		x.set(0,0,start.get(start.numRows()-2,0));
		x.set(1,0,start.get(start.numRows()-1,0));
		ArrayList<Double> mahalanobisDistances;
		double step = START_STEP_SIZE;
		double probability = 0;
		SimpleMatrix gradStep = null;
		int count =0;
		do {
			mahalanobisDistances = mahalanobis(x, means, covs);
			double prob = 0;
			// this loop calculates gradient and hessian as well as probability at x
			for (int i = 0; i < means.size(); i++) {
				// check whether the component actually contributes to to the density at given point by mahalanobis distance
				if(mahalanobisDistances.get(i) < MAX_MAHALANOBIS_DIST) {
					SimpleMatrix m = means.get(i);
					SimpleMatrix dm = m.minus(x);
					SimpleMatrix c = covs.get(i);
					SimpleMatrix invC = c.invert();
					double w = weights.get(i);
					//probability p(x,m) under component m
					double p = ((1 / (a * Math.sqrt(c.determinant()))) * Math.exp((-0.5d) * mahalanobisDistances.get(i))) * w;
					prob += p; 
					// gradient at x
					gradient = gradient.plus( invC.mult(dm).scale(p) );
					// hessian at x
					hessian = hessian.plus( invC.mult( dm.mult(dm.transpose()).minus(c) ).mult(invC).scale(p) );
				}


			}
			// save x
			SimpleMatrix xOld = new SimpleMatrix(x);
			double tst = evaluate(xOld, means, covs, weights);
			// check if hessian is negative definite
			SimpleEVD hessianEVD = hessian.eig();
			int maxEVIndex = hessianEVD.getIndexMax();
			// try a direct leap by quadratic optimization
			if(hessianEVD.getEigenvalue(maxEVIndex).getReal() < 0){
				gradStep = hessian.invert().mult(gradient);
				x = xOld.minus(gradStep);
			}
			double prob1 = 	evaluate(x, means, covs, weights);
			// if quadratic optimization did not work try gradient ascent
			if( prob1 <= prob | hessianEVD.getEigenvalue(maxEVIndex).getReal() >= 0) {
				gradStep = gradient.scale(step);
				x = xOld.plus(gradStep);
				// if still not ok decrease step size iteratively
				while(evaluate(x, means, covs, weights) < prob){
					step = step/2;
					gradStep = gradient.scale(step);
					x = xOld.plus(gradStep);
				}
			}
			probability =	model.evaluate(x, means, covs, weights);
			count++;
			// continue until the last step is sufficiently small or
			// a predefined amount of steps was performed
		}while(gradStep.elementMaxAbs() > STOP_STEP_SIZE && count<10);

		// return results
		return new SearchResult(x, probability);
	}
	
	/**
	 * Evaluate a gaussian mixture defined by the given means, covariances and component weights at a given point x.
	 * 
	 * @param x The point where the mixture shall be evaluated.
	 * @param means The component means.
	 * @param covs The component covariances.
	 * @param weights The component weights.
	 * @return The probability at point x.
	 */
	private static double evaluate(SimpleMatrix x, List<SimpleMatrix> means, List<SimpleMatrix> covs, List<Double> weights){
		ArrayList<Double> mahalanobisDistances = mahalanobis(x, means, covs);
		double n = x.numRows();
		double a = Math.pow(Math.sqrt(2 * Math.PI), n);

		double prob = 0;
		for (int i = 0; i < means.size(); i++) {
			// check wether the component actually contributes to to the density at given point 
			if(mahalanobisDistances.get(i) < MAX_MAHALANOBIS_DIST) {
				SimpleMatrix m = means.get(i);
				SimpleMatrix c = covs.get(i);
				double w = weights.get(i);
				//probability p(x,m) under component m
				double p = ((1 / (a * Math.sqrt(c.determinant()))) * Math.exp((-0.5d) * mahalanobisDistances.get(i))) * w;
				prob += p; 
			}
		}
		return prob;
	}
	
	/**
	 * Calculcates Mahalanobis distance between given point x and all given components defined by their means and covariances.
	 * 
	 * @param x The reference point.
	 * @param means The component means.
	 * @param covs The component covariances.
	 * @return A list with Mahalanobis distances between x and the given components.
	 */
	private static ArrayList<Double> mahalanobis(SimpleMatrix x, List<SimpleMatrix> means, List<SimpleMatrix> covs) {
		ArrayList<Double> mahalanobisDistances = new java.util.ArrayList<Double>();
		for (int i = 0; i < means.size(); i++) {
			SimpleMatrix m = means.get(i);
			SimpleMatrix c = covs.get(i);
			// calculate Mahalanobis distance
			double distance = x.minus(m).transpose().mult(c.invert()).mult(x.minus(m)).trace();
			mahalanobisDistances.add(distance);
		}
		return mahalanobisDistances;
	}
}
