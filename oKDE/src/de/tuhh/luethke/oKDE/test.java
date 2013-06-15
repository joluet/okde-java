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
		}
		
		//double I = SampleDist.getIntSquaredHessian(mu, w, cov, g);
		
		//System.out.println(I);
	}
}
