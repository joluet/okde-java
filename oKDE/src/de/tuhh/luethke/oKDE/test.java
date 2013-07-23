package de.tuhh.luethke.oKDE;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Writer;
import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;

import org.ejml.simple.SimpleMatrix;

import de.tuhh.luethke.oKDE.Exceptions.EmptyDistributionException;
import de.tuhh.luethke.oKDE.model.BaseSampleDistribution;
import de.tuhh.luethke.oKDE.model.SampleModel;
import de.tuhh.luethke.oKDE.utility.Compressor;
import de.tuhh.luethke.oKDE.utility.MomentMatcher;
import de.tuhh.luethke.oKDE.utility.ProjectionData;
import de.tuhh.luethke.oKDE.utility.Projector;

public class test {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		System.out.println("Start Testing!");
		SampleModel dist = new SampleModel();
		
		System.out.println("Test intsqrd");
		
		//  1 1.5 1.5 5 5.5 4.5 5 4.5;  1 1 1.5 5 5.5 4.5 4.5 5];
		//double[][] dMu = {{1.18835013453899,0.0698038241915721,-1.25815395873056},{0.766696783442554,-1.41248979682270,0.645793013380147}};
		double[][] mean1 = {{1},{1}};
		double[][] mean2 = {{1},{1}};
		double[][] mean3 = {{1.523},{1}};
		double[][] mean4 = {{1.5},{1.5}};
		double[][] mean5 = {{5},{5}};
		double[][] mean6 = {{5.5},{5.5}};
		double[][] mean7 = {{4.5},{4.5}};
		double[][] mean8 = {{5},{4.5}};
		double[][] mean9 = {{4.5},{5}};
		double[][] mean10 = {{1},{1}};
		double[][] mean11 = {{1.5},{1}};
		double[][] mean12 = {{1.5},{1.5}};
		double[][] mean13 = {{5},{5}};
		double[][] mean14 = {{5.5},{5.5}};
		double[][] mean15 = {{4.5},{4.5}};
		double[][] mean16 = {{5},{4.5}};
		double[][] mean17 = {{4.5},{5.1222}};
		/*
		double[][] mean7 = {{5},{4.5}};
		double[][] mean8 = {{4.5},{5}};
		double[][] mean9 = {{5},{5}};
		double[][] mean10 = {{5},{5}};*/
		
		
		//SimpleMatrix mu = new SimpleMatrix(dMu);
		//SimpleMatrix[] mus = {mu};
		double[] w= {1,1,1};
		double[][] c = {{0.0,0.0}, {0.0,0.0}};
		SimpleMatrix[] cov = {new SimpleMatrix(c),new SimpleMatrix(c),new SimpleMatrix(c)};
		//double[][] dG = {{0.693361274350635,0.0}, {0.0,0.693361274350635}};
		//SimpleMatrix g = new SimpleMatrix(dG);
		//System.out.println(mu);

		try {
			ArrayList<SimpleMatrix>  means = new ArrayList<SimpleMatrix>();
			means.add(new SimpleMatrix(mean1));
			means.add(new SimpleMatrix(mean2));
			means.add(new SimpleMatrix(mean3));
		    dist.updateDistribution(means.toArray(new SimpleMatrix[3]), cov, w);
		    
		    means = new ArrayList<SimpleMatrix>();
			means.add(new SimpleMatrix(mean4));
			means.add(new SimpleMatrix(mean5));
			means.add(new SimpleMatrix(mean6));
			means.add(new SimpleMatrix(mean7));
			means.add(new SimpleMatrix(mean8));
			means.add(new SimpleMatrix(mean9));
			means.add(new SimpleMatrix(mean10));
			means.add(new SimpleMatrix(mean11));
			means.add(new SimpleMatrix(mean12));
			means.add(new SimpleMatrix(mean13));
			means.add(new SimpleMatrix(mean14));
			means.add(new SimpleMatrix(mean15));
			means.add(new SimpleMatrix(mean16));
			means.add(new SimpleMatrix(mean17));
		    
		    for(int i=0; i<14; i++) {
		    	dist.updateDistribution(means.get(i), new SimpleMatrix(c), 1.0d);
		    }
		    
		    

		} catch (EmptyDistributionException e) {
		    // TODO Auto-generated catch block
		    e.printStackTrace();
		} catch (IllegalArgumentException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (InvocationTargetException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (NoSuchMethodException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (SecurityException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (InstantiationException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IllegalAccessException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		System.out.println("---------------------------------------");
		/*dist = new SampleDist();
		double[][] c = {{0.1,1.1}, {2.0,1.0}};
		SimpleMatrix cov = new SimpleMatrix(c);
		dist.setGlobalCovariance(cov);
		dist.setmBandwidthMatrix(cov);*/
		System.out.println("BW"+dist.getBandwidthMatrix());

		/*for(int i=0; i<dist.getSubCovariances().size(); i++){
			dist.setGlobalCovariance(new SimpleMatrix(dist.getmBandwidthMatrix()));
		}*/

		
		for(BaseSampleDistribution d : dist.getSubDistributions())
			System.out.println("subdist: "+ d.getGlobalMean() +  "weight " + d.getGlobalWeight());
		System.out.println("finished");

		/*ArrayList<SimpleMatrix> cov1 = new ArrayList<SimpleMatrix>();
		double[][] c = {{0.0,0.0}, {0.0,0.0}};
		cov1.add(new SimpleMatrix(c));
		cov1.add(new SimpleMatrix(c));
		cov1.add(new SimpleMatrix(c));
		cov1.add(new SimpleMatrix(c));
		cov1.add(new SimpleMatrix(c));
		cov1.add(new SimpleMatrix(c));
		cov1.add(new SimpleMatrix(c));
		cov1.add(new SimpleMatrix(c));
		cov1.add(new SimpleMatrix(c));
		cov1.add(new SimpleMatrix(c));
		
		double[][] mean1 = {{1},{1}};
		double[][] mean2 = {{1.5},{1}};
		double[][] mean3 = {{1.5},{1.5}};
		double[][] mean4 = {{5},{5}};
		double[][] mean5 = {{5.5},{5.5}};
		double[][] mean6 = {{4.5},{4.5}};
		double[][] mean7 = {{5},{4.5}};
		double[][] mean8 = {{4.5},{5}};
		double[][] mean9 = {{5},{5}};
		double[][] mean10 = {{5},{5}};
		ArrayList<SimpleMatrix> means = new ArrayList<SimpleMatrix>();
		means.add(new SimpleMatrix(mean1));
		means.add(new SimpleMatrix(mean2));
		means.add(new SimpleMatrix(mean3));
		means.add(new SimpleMatrix(mean4));
		means.add(new SimpleMatrix(mean5));
		means.add(new SimpleMatrix(mean6));
		means.add(new SimpleMatrix(mean7));
		means.add(new SimpleMatrix(mean8));
		means.add(new SimpleMatrix(mean9));
		means.add(new SimpleMatrix(mean10));
		//dist.setMeans(means);
		//dist.setCovariances(cov1);
		double[] weights = new double[10];
		weights[0] = 1d;
		weights[1] = 1d;
		weights[2] = 1d;
		weights[3] = 1d;
		weights[4] = 1d;
		weights[5] = 1d;
		weights[6] = 1d;
		weights[7] = 1d;
		weights[8] = 1d;
		weights[9] = 1d;*/
		
		/*
		double[] weights = new double[15];
		ArrayList<SimpleMatrix> means = new ArrayList<SimpleMatrix>();
		ArrayList<SimpleMatrix> cov1 = new ArrayList<SimpleMatrix>();
		
		for(int i=0; i<15; i++) {
		    if(i<100){
			double d1 = StdRandom.gaussian(2, 1);
			double d2 = StdRandom.gaussian(2, 1);
			double[][] mean = {{d1},{d2}};
			means.add(new SimpleMatrix(mean));
			double[][] c = {{0.0,0.0}, {0.0,0.0}};
			cov1.add(new SimpleMatrix(c));
			weights[i] = 1d;
		    }else if(i<10){
			double d1 = StdRandom.gaussian(2, 1);
			double d2 = StdRandom.gaussian(7, 1);
			double[][] mean = {{d1},{d2}};
			means.add(new SimpleMatrix(mean));
			double[][] c = {{0.0,0.0}, {0.0,0.0}};
			cov1.add(new SimpleMatrix(c));
			weights[i] = 1d;
		    }else{
			double d1 = StdRandom.gaussian(7, 1);
			double d2 = StdRandom.gaussian(2, 1);
			double[][] mean = {{d1},{d2}};
			means.add(new SimpleMatrix(mean));
			double[][] c = {{0.0,0.0}, {0.0,0.0}};
			cov1.add(new SimpleMatrix(c));
			weights[i] = 1d;
		    }
		}
		dataToFile(means);*/
		
		
		/*means = readFromFile();
		for(int i=0; i<1000; i++) {
		    double[][] c = {{0.0,0.0}, {0.0,0.0}};
		    cov1.add(new SimpleMatrix(c));
		    weights[i] = 1d;
		 }*/
		//dist.setWeights(weights);
		/*try {
		    SampleDist.projectToSubspace(dist);
		} catch (EmptyDistributionException e) {
		    // TODO Auto-generated catch block
		    e.printStackTrace();
		}*/
		/*try {
		    dist.updateDistribution(means.toArray(new SimpleMatrix[0]), cov1.toArray(new SimpleMatrix[0]), weights);
		} catch (EmptyDistributionException e) {
		    // TODO Auto-generated catch block
		    e.printStackTrace();
		}*/
		
        	// define your data
        	/*double[] x = new double[100];
        	double[] y = new double[100];
        	double coord = 0;
        	for(int i=0; i<100; i++) {
        	    coord += 0.1;
        	    x[i] = coord;
        	}
        	coord = 0;
        	for(int i=0; i<100; i++) {
        	    coord += 0.1;
        	    y[i] = coord;
        	}
        	double[][] z1 = new double[300][300];
        	z1 = calculateY(x ,y, dist, means);
        	// create your PlotPanel (you can use it as a JPanel) with a legend at
        	// SOUTH
        	Plot3DPanel plot = new Plot3DPanel("SOUTH");
        
        	// add grid plot to the PlotPanel
        	plot.addGridPlot("kernel", x, y, z1);
        
        	// put the PlotPanel in a JFrame like a JPanel
        	JFrame frame = new JFrame("a plot panel");
        	frame.setSize(600, 600);
        	frame.setContentPane(plot);
        	frame.setVisible(true);*/
        	
        
        	// double I = SampleDist.getIntSquaredHessian(mu, w, cov, g);
		
		//System.out.println(I);
	}
	private static double calculateY(double x, double y, BaseSampleDistribution dist, ArrayList<SimpleMatrix> means){
		//ArrayList<SimpleMatrix> means = dist.getMeans();
		SimpleMatrix bandwidth = dist.getBandwidthMatrix();
	    	//double[][] c = {{1.381236356118853, 1.375256977953836},{1.375256977953837, 1.387215734283870}};
	    	//SimpleMatrix bandwidth = new SimpleMatrix(c);
		double[][] dxVector = {{x},{y}};
		SimpleMatrix xVector = new SimpleMatrix(dxVector);
		double d = 0d;
		double n = means.size();
		for(SimpleMatrix m : means) {
		    double tmp = (-0.5d)*xVector.minus(m).transpose().mult(bandwidth.invert()).mult(xVector.minus(m)).trace();
		    d += ( ( 1 / Math.sqrt(4*Math.PI*Math.PI*bandwidth.determinant()) ) * Math.exp(tmp) );
		}
		return d/n;
	}
    
        /*public static double calculateY(double x, double y, SampleDist dist, ArrayList<SimpleMatrix> means) {
        	double z = cos(x * PI) * sin(y * PI);
        	return z;
        }*/
    
        public static double[][] calculateY(double[] x, double[] y, BaseSampleDistribution dist, ArrayList<SimpleMatrix> means) {
        	double[][] z = new double[y.length][x.length];
        	for (int i = 0; i < x.length; i++)
        	    for (int j = 0; j < y.length; j++)
        		z[j][i] = calculateY(x[i], y[j], dist, means);
        	return z;
        }
        
        private static void dataToFile(ArrayList<SimpleMatrix> data) {
            PrintWriter pw = null;
            try
            {
              Writer fw = new FileWriter( "data.txt" );
              Writer bw = new BufferedWriter( fw );
              pw = new PrintWriter( bw );
              pw.print("[");
              for(int i=0; i<data.get(0).numRows(); i++) {
                  for (SimpleMatrix m : data)
                      pw.print(m.get(i, 0)+" ");
                  if(i<(data.get(0).numRows()-1))
                      pw.print(";");
              }
              pw.print("]");
            }
              
            catch ( IOException e ) {
              System.err.println( "Error creating file!" );
            }
            finally {
              if ( pw != null )
                pw.close();
            }
        }
        
        private static ArrayList<SimpleMatrix> readFromFile() {
        	BufferedReader br = null;
        	String[] matrix;
        	String[] matrixRow1 = null;
        	String[] matrixRow2 = null;

        	try {
        	    String data = "";
        	    String sCurrentLine;
        
        	    br = new BufferedReader(new FileReader("data.txt"));
        
        	    while ((sCurrentLine = br.readLine()) != null) {
        		data += sCurrentLine;
        	    }
        	    data = data.substring(1);
        	    data = data.substring(0, data.length()-1);
        	    matrix = data.split(";");
        	    matrixRow1 = matrix[0].split(" ");
        	    matrixRow2 = matrix[1].split(" ");   
        	    
        	} catch (IOException e) {
        	    e.printStackTrace();
        	} finally {
        	    try {
        		if (br != null)
        		    br.close();
        	    } catch (IOException ex) {
        		ex.printStackTrace();
        	    }
        	}
        	ArrayList<SimpleMatrix> dataArray = new ArrayList<SimpleMatrix>();
        	for(int i=0; i<matrixRow1.length; i++){
        	    double[][] d = {{Double.parseDouble(matrixRow1[i])},{Double.parseDouble(matrixRow2[i])}};
        	    dataArray.add(new SimpleMatrix(d));
        	}
        	return dataArray;
        	
        }
}
