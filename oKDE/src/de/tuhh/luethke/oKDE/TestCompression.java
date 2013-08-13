package de.tuhh.luethke.oKDE;

import java.util.ArrayList;

import org.ejml.simple.SimpleMatrix;

import de.tuhh.luethke.oKDE.Exceptions.EmptyDistributionException;
import de.tuhh.luethke.oKDE.model.BaseSampleDistribution;
import de.tuhh.luethke.oKDE.utility.Compression.Hellinger;

public class TestCompression {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		BaseSampleDistribution distribution = new BaseSampleDistribution();
		double[][] mean1 = {{0},{0}};
		double[][] mean2 = {{0.810488643593117},{0.658276986625753}};
		double[][] mean3 = {{0.810488643593117},{0.658276986625753}};
		double[][] mean4 = {{0.700923094831665},{-1.96374467859908}};
		double[][] mean5 = {{0.587934445337502},{0.578840296979537}};
		double[][] mean6 = {{-0.969944942451800},{0.0227834694560141}};
		double[][] mean7 = {{-1.19249914070741},{-0.0566532201902028}};
		double[][] mean8 = {{0.747390744196186},{0.102220159102233}};
		ArrayList<SimpleMatrix> means = new ArrayList<SimpleMatrix>();
		means.add(new SimpleMatrix(mean1));
		means.add(new SimpleMatrix(mean2));
		means.add(new SimpleMatrix(mean3));
		means.add(new SimpleMatrix(mean4));
		means.add(new SimpleMatrix(mean5));
		means.add(new SimpleMatrix(mean6));
		means.add(new SimpleMatrix(mean7));
		means.add(new SimpleMatrix(mean8));
		double[] w= {1,1,1,1,1,1,1,1};
		double[][] c = {{0.0,0.0}, {0.0,0.0}};
		SimpleMatrix[] cov = {new SimpleMatrix(c), new SimpleMatrix(c),new SimpleMatrix(c),new SimpleMatrix(c),new SimpleMatrix(c),new SimpleMatrix(c),new SimpleMatrix(c),new SimpleMatrix(c)};
		try {
			distribution.updateDistribution(means.toArray(new SimpleMatrix[8]), cov, w);
		} catch (EmptyDistributionException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		try {
			Hellinger.getAllSigmaPoints(distribution, 3);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

}
