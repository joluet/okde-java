package de.tuhh.luethke.oKDE.model;

import static org.junit.Assert.*;

import java.util.ArrayList;

import org.ejml.simple.SimpleMatrix;
import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

import de.tuhh.luethke.oKDE.Exceptions.EmptyDistributionException;

public class SampleDistTest {
    
    @BeforeClass
    public static void setUpBeforeClass() throws Exception {
	
    }

    @AfterClass
    public static void tearDownAfterClass() throws Exception {
    }

    @Test
    public void testUpdateDistribution() {
	fail("Not yet implemented");
    }

    @Test
    public void testProjectToSubspace() {
	SampleDist testDist = new SampleDist();
	//First test set
	ArrayList<SimpleMatrix> cov1 = new ArrayList<SimpleMatrix>();
	double[][] c = {{0.0,0.0}, {0.0,0.0}};
	cov1.add(new SimpleMatrix(c));
	cov1.add(new SimpleMatrix(c));
	cov1.add(new SimpleMatrix(c));
	double[][] mean1 = {{1},{1}};
	double[][] mean2 = {{1.5},{1.5}};
	double[][] mean3 = {{2},{2}};
	ArrayList<SimpleMatrix> means = new ArrayList<SimpleMatrix>();
	means.add(new SimpleMatrix(mean1));
	means.add(new SimpleMatrix(mean2));
	means.add(new SimpleMatrix(mean3));
	ArrayList<Double> weights = new ArrayList<Double>();
	weights.add(1d);
	weights.add(1d);
	weights.add(1d);
	
	testDist.setMeans(means);
	testDist.setCovariances(cov1);
	testDist.setWeights(weights);
	try {
	    SampleDist.projectToSubspace(testDist);
	} catch (EmptyDistributionException e) {
	    // TODO Auto-generated catch block
	    e.printStackTrace();
	}
    }

}
