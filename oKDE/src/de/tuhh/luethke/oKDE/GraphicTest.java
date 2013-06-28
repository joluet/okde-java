package de.tuhh.luethke.oKDE;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Insets;
import java.awt.geom.Line2D;
import java.awt.geom.Point2D;
import java.util.ArrayList;
import java.util.Random;

import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.SwingUtilities;

import org.ejml.simple.SimpleMatrix;

import de.tuhh.luethke.oKDE.Exceptions.EmptyDistributionException;
import de.tuhh.luethke.oKDE.model.SampleDist;

public class GraphicTest extends JFrame {

    public GraphicTest() {

        initUI();
    }

    private void initUI() {
        
        setTitle("Points");
        setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

        add(new Surface());

        setSize(900, 500);
        setLocationRelativeTo(null);
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

class Surface extends JPanel {

    private void doDrawing(Graphics g) {
	
	
	SampleDist dist = new SampleDist();
	ArrayList<SimpleMatrix> cov1 = new ArrayList<SimpleMatrix>();
	double[][] c = {{0.0}};
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
	ArrayList<SimpleMatrix> means = new ArrayList<SimpleMatrix>();
	means.add(new SimpleMatrix(mean1));
	means.add(new SimpleMatrix(mean2));
	means.add(new SimpleMatrix(mean3));
	means.add(new SimpleMatrix(mean4));
	means.add(new SimpleMatrix(mean5));
	means.add(new SimpleMatrix(mean6));
	//dist.setMeans(means);
	//dist.setCovariances(cov1);
	double[] weights = new double[6];
	weights[0] = 1d;
	weights[1] = 1d;
	weights[2] = 1d;
	weights[3] = 1d;
	weights[4] = 1d;
	weights[5] = 1d;
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
	
	

        Graphics2D g2d = (Graphics2D) g;

        g2d.setColor(Color.blue);

        Dimension size = getSize();
        Insets insets = getInsets();

        int w = size.width - insets.left - insets.right;
        int h = size.height - insets.top - insets.bottom;

        Random r = new Random();
        double x = 0;
        for (int i = 0; i < 1500; i++) {
            x += 0.01;
            double y = calculateY(x, dist, means)*500d;
            //Point2D.Double point = new Point2D.Double(x, y);
            //System.out.println("x "+x+" | y "+y);
            g2d.draw(new Line2D.Double(x*50d,y,x*50d,y));
            System.out.println("x"+x*10d+" y"+y);
        }
        for (int i = 0; i < 75; i+=15) {
            //Point2D.Double point = new Point2D.Double(x, y);
            //System.out.println("x "+x+" | y "+y);
           // g2d.draw(new Line2D.Double(i,1,i,1));
            g2d.drawString(String.valueOf(i), i, 20);
        }
    }

    @Override
    public void paintComponent(Graphics g) {

        super.paintComponent(g);
        doDrawing(g);
    }
    
    private double calculateY(double x, SampleDist dist, ArrayList<SimpleMatrix> means){
	//ArrayList<SimpleMatrix> means = dist.getMeans();
	SimpleMatrix bandwidth = dist.getmBandwidthMatrix();
	double h = bandwidth.get(0,0);
	double y = 0d;
	int n = means.size();
	for(SimpleMatrix m : means) {
	    double mValue = m.get(0, 0);
	    y += (1d/(n*h))*(Gaussian.phi(x, mValue, h));
	}
	return y;
    }
}
