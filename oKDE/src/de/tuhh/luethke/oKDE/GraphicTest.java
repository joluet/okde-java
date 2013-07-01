package de.tuhh.luethke.oKDE;

import java.util.ArrayList;

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.ejml.simple.SimpleMatrix;
import org.math.plot.Plot2DPanel;
import org.math.plot.Plot3DPanel;

import de.tuhh.luethke.oKDE.Exceptions.EmptyDistributionException;
import de.tuhh.luethke.oKDE.model.SampleDist;

public class GraphicTest extends JFrame {

    public GraphicTest() {

        initUI();
    }

    private void initUI() {
        
        setTitle("Points");
        setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);


        setSize(900, 500);
        setLocationRelativeTo(null);
        
        
        SampleDist dist = new SampleDist();
	ArrayList<SimpleMatrix> cov1 = new ArrayList<SimpleMatrix>();
	double[][] c = {{0.0}};
	cov1.add(new SimpleMatrix(c));
	cov1.add(new SimpleMatrix(c));
	cov1.add(new SimpleMatrix(c));
	cov1.add(new SimpleMatrix(c));
	cov1.add(new SimpleMatrix(c));
	cov1.add(new SimpleMatrix(c));
	cov1.add(new SimpleMatrix(c));
	double[][] mean1 = {{1}};
	double[][] mean2 = {{1.5}};
	double[][] mean3 = {{0.5}};
	double[][] mean4 = {{9}};
	double[][] mean5 = {{9}};
	double[][] mean6 = {{9}};
	double[][] mean7 = {{11}};
	ArrayList<SimpleMatrix> means = new ArrayList<SimpleMatrix>();
	means.add(new SimpleMatrix(mean1));
	means.add(new SimpleMatrix(mean2));
	means.add(new SimpleMatrix(mean3));
	means.add(new SimpleMatrix(mean4));
	means.add(new SimpleMatrix(mean5));
	means.add(new SimpleMatrix(mean6));
	means.add(new SimpleMatrix(mean7));
	//dist.setMeans(means);
	//dist.setCovariances(cov1);
	double[] weights = new double[7];
	weights[0] = 1d;
	weights[1] = 1d;
	weights[2] = 1d;
	weights[3] = 1d;
	weights[4] = 1d;
	weights[5] = 1d;
	weights[6] = 1d;
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
        
        double x = 0;
        double[] xArray = new double[1500];
        double[] yArray = new double[1500];
        for (int i = 0; i < 1500; i++) {
            x += 0.01;
            xArray[i] = x;
            double y = calculateY(x, dist, means);
            yArray[i] = y;
            //Point2D.Double point = new Point2D.Double(x, y);
            //System.out.println("x "+x+" | y "+y);
            //g2d.draw(new Line2D.Double(x*50d,y,x*50d,y));
            //System.out.println("x"+x*10d+" y"+y);
        }
	// create your PlotPanel (you can use it as a JPanel)
	Plot3DPanel plot = new Plot3DPanel();
	 
	// add a line plot to the PlotPanel
	plot.add.addLinePlot("my plot", xArray, yArray);
	
	setContentPane(plot);
    }
    
    private double calculateY(double x, SampleDist dist, ArrayList<SimpleMatrix> means){
	//ArrayList<SimpleMatrix> means = dist.getMeans();
	SimpleMatrix bandwidth = dist.getmBandwidthMatrix();
	double h =bandwidth.get(0,0);
	double y = 0d;
	int n = means.size();
	for(SimpleMatrix m : means) {
	    double mValue = m.get(0, 0);
	    y += (1d/(n*h))*(Gaussian.phi(x, mValue, h));
	}
	return y;
    }

    public static void main(String[] args) {

        SwingUtilities.invokeLater(new Runnable() {
            @Override
            public void run() {

        	GraphicTest ps = new GraphicTest();
                ps.setVisible(true);
            }
        });
    }

}


    
    
