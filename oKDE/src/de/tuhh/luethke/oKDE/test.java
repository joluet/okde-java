package de.tuhh.luethke.oKDE;

import java.util.ArrayList;

import org.ejml.factory.SingularMatrixException;
import org.ejml.simple.SimpleMatrix;

import de.tuhh.luethke.oKDE.Exceptions.EmptyDistributionException;
import de.tuhh.luethke.oKDE.model.SampleDist;

public class test {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		System.out.println("Start Testing!");
		SampleDist dist = new SampleDist();
		
		System.out.println("Test intsqrd");
		double[][] dMu = {{1.18835013453899,0.0698038241915721,-1.25815395873056},{0.766696783442554,-1.41248979682270,0.645793013380147}};
		
		
		
		/*
		SimpleMatrix mu = new SimpleMatrix(dMu);
		SimpleMatrix[] mus = {mu};
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

		try {
		    dist.updateDistribution(mus, cov, w);
		} catch (EmptyDistributionException e) {
		    // TODO Auto-generated catch block
		    e.printStackTrace();
		}*/
		System.out.println("---------------------------------------");
		dist = new SampleDist();
		ArrayList<SimpleMatrix> cov1 = new ArrayList<SimpleMatrix>();
		double[][] c = {{0.0,0.0}, {0.0,0.0}};
		cov1.add(new SimpleMatrix(c));
		cov1.add(new SimpleMatrix(c));
		cov1.add(new SimpleMatrix(c));
		cov1.add(new SimpleMatrix(c));
		cov1.add(new SimpleMatrix(c));
		double[][] mean1 = {{1},{0}};
		double[][] mean2 = {{1.5},{1.5}};
		double[][] mean3 = {{9},{8}};
		double[][] mean4 = {{9},{8}};
		double[][] mean5 = {{9},{8}};
		ArrayList<SimpleMatrix> means = new ArrayList<SimpleMatrix>();
		means.add(new SimpleMatrix(mean1));
		means.add(new SimpleMatrix(mean2));
		means.add(new SimpleMatrix(mean3));
		means.add(new SimpleMatrix(mean4));
		means.add(new SimpleMatrix(mean5));
		//dist.setMeans(means);
		//dist.setCovariances(cov1);
		double[] weights = new double[5];
		weights[0] = 1d;
		weights[1] = 1d;
		weights[2] = 1d;
		weights[3] = 1d;
		weights[4] = 1d;
		//dist.setWeights(weights);
		/*try {
		    SampleDist.projectToSubspace(dist);
		} catch (EmptyDistributionException e) {
		    // TODO Auto-generated catch block
		    e.printStackTrace();
		}*/
		try {
		    dist.updateDistribution(means.toArray(new SimpleMatrix[0]), cov1.toArray(new SimpleMatrix[0]), weights);
		} catch (EmptyDistributionException e) {
		    // TODO Auto-generated catch block
		    e.printStackTrace();
		}
		
		//double I = SampleDist.getIntSquaredHessian(mu, w, cov, g);
		
		//System.out.println(I);
	}
}
