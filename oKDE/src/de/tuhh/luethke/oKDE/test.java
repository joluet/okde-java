package de.tuhh.luethke.oKDE;

import java.util.ArrayList;

import org.ejml.factory.SingularMatrixException;
import org.ejml.simple.SimpleMatrix;

import de.tuhh.luethke.oKDE.model.SampleDist;

public class test {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		System.out.println("Start Testing!");
		ArrayList<SimpleMatrix> sampleData = new ArrayList<SimpleMatrix>();
		ArrayList<Double> weights = new ArrayList<Double>();
		SimpleMatrix[] covariances = new SimpleMatrix[1];
		SampleDist dist = new SampleDist(sampleData, weights, covariances);
		double d[] = {1d}; 
		dist.addWeights(d);
		dist.addWeights(d);
		dist.addWeights(d);
		dist.addWeights(d);
		dist.addWeights(d);
		
		
		System.out.println("Test intsqrd");
		double[][] dMu = {{1.18835013453899,0.0698038241915721,-1.25815395873056},{0.766696783442554,-1.41248979682270,0.645793013380147}};
		
		
		
		
		SimpleMatrix mu = new SimpleMatrix(dMu);
		double[] w= {0.333333333333333,0.333333333333333,0.333333333333333};
		double[][] c = {{0.0,0.0}, {0.0,0.0}};
		SimpleMatrix[] cov = {new SimpleMatrix(c), new SimpleMatrix(c), new SimpleMatrix(c)};
		//double[][] dG = {{0.693361274350635,0.0}, {0.0,0.693361274350635}};
		//SimpleMatrix g = new SimpleMatrix(dG);
		System.out.println(mu);
		System.out.println(cov[0]);
		//System.out.println(g);
		double[][] cov_smp = {{1,0.0}, {0.0,1}};
		SimpleMatrix Cov_smp = new SimpleMatrix(cov_smp);
		System.out.println(Cov_smp);

		SampleDist.estimateBandwidth(mu, cov, w, Cov_smp, 3);
		
		//double I = SampleDist.getIntSquaredHessian(mu, w, cov, g);
		
		//System.out.println(I);
	}
}
