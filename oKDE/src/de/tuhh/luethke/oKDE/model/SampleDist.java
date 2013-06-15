package de.tuhh.luethke.oKDE.model;

import java.util.ArrayList;

import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;
import org.ejml.simple.SimpleMatrix;

import de.tuhh.luethke.oKDE.Exceptions.EmptyDistributionException;

public class SampleDist {

    // component means
    private ArrayList<SimpleMatrix> mMeans;

    // component covariances
    private ArrayList<SimpleMatrix> mCovariances;

    // component weights
    private ArrayList<Double> mWeights;

    private double mWeightSum;

    private double mForgettingFactor;

    // effective number of obesverved samples
    private double mN_eff;
    

    public SampleDist() {
	super();
	mMeans = new ArrayList<SimpleMatrix>();
	mCovariances = new ArrayList<SimpleMatrix>();
	mWeightSum = 0;
	mWeights = new ArrayList<Double>();
	mN_eff = 0;
	mForgettingFactor = 1;
    }
    
    /*public SampleDist(ArrayList<Double> weights, ArrayList<SimpleMatrix> means,
	    ArrayList<SimpleMatrix> covariances) {
	super();
	mMeans = means;
	mCovariances = covariances;
	mWeightSum = 0;
	for (double w : weights) {
	    mWeightSum += w;
	}
	mWeights = weights;
	mN_eff = weights.size();
	mForgettingFactor = 1;
    }*/

    public ArrayList<Double> getWeights() {
	return mWeights;
    }

    /**
     * Update the sample distribution adding a single component or a batch of
     * components. This method executes two main steps: 1. Augment sample
     * distribution and recalculate weights 2. Reestimate the bandwidth
     * 
     * @param means
     * @param covariances
     * @param weights
     * @throws EmptyDistributionException
     *             Exception is thrown when one parameter is null or empty
     */
    public void updateDistribution(SimpleMatrix[] means,
	    SimpleMatrix[] covariances, double[] weights)
	    throws EmptyDistributionException {
	// at first check input parameters!
	checkInputParams(means, covariances, weights);

	// augment distribution
	addWeights(weights);
	for (SimpleMatrix m : means)
	    mMeans.add(m);
	for (SimpleMatrix c : covariances)
	    mCovariances.add(c);

	// get empirical sample covariance by moment matching
	SimpleMatrix covSmpl = projectToSubspace(means, covariances, weights,
		mN_eff);
	// reestimate bandwidth as explained in oKDE paper
	reestimateBandwidth(means, covariances, weights, covSmpl, mN_eff);
    }

    private void checkInputParams(SimpleMatrix[] means,
	    SimpleMatrix[] covariances, double[] weights)
	    throws EmptyDistributionException {
	if (weights == null || weights.length == 0)
	    throw new EmptyDistributionException();
	if (means == null || means.length == 0)
	    throw new EmptyDistributionException();
	if (covariances == null || covariances.length == 0)
	    throw new EmptyDistributionException();
    }

    /**
     * Takes new incoming sample weights and updates this distribution using a
     * forgetting factor.
     * 
     * @param weights
     */
    private void addWeights(double[] weights) {
	double sumOfNewWeights = 0;
	for (double w : weights) {
	    sumOfNewWeights += w;
	    mWeights.add(w);
	}

	mN_eff = mN_eff * mForgettingFactor + weights.length;

	// calculate mixing weights for old and new weights
	double mixWeightOld = mWeightSum
		/ (mWeightSum * mForgettingFactor + sumOfNewWeights);
	double mixWeightNew = sumOfNewWeights
		/ (mWeightSum * mForgettingFactor + sumOfNewWeights);

	mWeightSum = mWeightSum * mForgettingFactor + sumOfNewWeights;

	for (int i = 0; i < mWeights.size() - weights.length; i++) {
	    mWeights.set(i, mWeights.get(i) * mixWeightOld);
	}
	for (int i = mWeights.size() - weights.length; i < mWeights.size(); i++) {
	    mWeights.set(i, mWeights.get(i) * mixWeightNew);
	}

	System.out.println(mixWeightOld + "-" + mixWeightNew + " " + mWeightSum
		+ " " + mWeights);
    }

    public ArrayList<SimpleMatrix> getMeans() {
	return mMeans;
    }

    public ArrayList<SimpleMatrix> getCovariances() {
	return mCovariances;
    }

    public double getWeightSum() {
	return mWeightSum;
    }

    private SimpleMatrix projectToSubspace(SimpleMatrix[] means,
	    SimpleMatrix[] covariances, double[] weights, double N_eff) {
	double minBW = 1e-7;
	/*
	 * SampleDist subSpaceDist = momentMatchPdf(dist) ;
	 * subSpaceDist.getCovariance().svd(false); [u,s,v] = svd( gC ) ; s =
	 * diag(s) ;
	 */
	return null;
    }

    private static void reestimateBandwidth(SimpleMatrix[] means,
	    SimpleMatrix[] covariance, double[] weights, SimpleMatrix Cov_smp,
	    double N_eff) {

	double d = means[0].numRows();

	SimpleMatrix G = Cov_smp.scale(Math.pow((4 / ((d + 2) * N_eff)),
		(2 / (d + 4))));

	float alphaScale = 1;
	SimpleMatrix F = Cov_smp.scale(alphaScale);

	double Rf2 = getIntSquaredHessian(means, weights, covariance, F, G);
	double hAmise = Math
		.pow((Math.pow(N_eff, (-1))
			* Math.pow(F.determinant(), (-1 / 2)) / (Math.pow(
			Math.sqrt(4 * Math.PI), d)
			* Rf2 * d)), (1 / (d + 4)));
	
	System.out.println("hAmise: " + hAmise);
    }

    private static double getIntSquaredHessian(SimpleMatrix[] means,
	    double[] weights, SimpleMatrix[] covariance, SimpleMatrix F,
	    SimpleMatrix g) {

	long d = means[0].numRows();
	long N = means.length;
	System.out.println("d:" + d);
	// normalizer
	double constNorm = Math.pow((1 / (2 * Math.PI)), (d / 2));

	double w1, w2, m, I = 0, eta, f_t, c;
	SimpleMatrix s1, s2, mu1, mu2, dm, ds, B, b, C;
	for (int i1 = 0; i1 < N; i1++) {
	    s1 = covariance[i1].plus(g);
	    mu1 = means[i1];
	    w1 = weights[i1];
	    for (int i2 = i1; i2 < N; i2++) {
		s2 = covariance[i2];
		mu2 = means[i2];
		w2 = weights[i2];
		SimpleMatrix A = s1.plus(s2).invert();
		dm = mu1.minus(mu2);

		// test if F is identity for speedup
		SimpleMatrix Id = SimpleMatrix.identity(F.numCols());
		double deltaF = F.minus(Id).elementSum();
		// if F is not identity
		if (deltaF > 1e-3) {
		    ds = dm.transpose().mult(A);
		    b = ds.transpose().mult(ds);
		    B = A.minus(b.scale(2));
		    C = A.minus(b);
		    f_t = constNorm * Math.sqrt(A.determinant())
			    * Math.exp(-0.5 * ds.mult(dm).trace());
		    c = 2 * F.mult(A).mult(F).mult(B).trace()
			    + Math.pow(F.mult(C).trace(), 2);
		} else {
		    m = dm.transpose().mult(A).mult(dm).get(0);
		    System.out.println("m: " + m);
		    System.out.println("A:" + A);
		    f_t = constNorm * Math.sqrt(A.determinant())
			    * Math.exp(-0.5 * m);

		    DenseMatrix64F A_sqr = new DenseMatrix64F(A.numRows(),
			    A.numCols());
		    CommonOps.elementMult(A.getMatrix(), A.transpose()
			    .getMatrix(), A_sqr);
		    double sum = CommonOps.elementSum(A_sqr);
		    c = 2 * sum * (1 - 2 * m) + Math.pow((1 - m), 2)
			    * Math.pow(A.trace(), 2);
		}

		// determine the weight of the current term
		if (i1 == i2)
		    eta = 1;
		else
		    eta = 2;
		System.out.println("I: " + I + " f_t: " + f_t + " c: " + c
			+ " w2*w1:" + (w2 * w1));
		I = I + f_t * c * w2 * w1 * eta;
		System.out.println("I: " + I);
	    }
	}

	return I;
    }
    /*
     * FASTER VERSION? private double getIntSquaredHessian(double[] muValues,
     * double[] wValues, double[] cValues, double[] gValue) { //double
     * *muValues, *cValues, *wValues, *gValue, tmp, *ptr_out ;
     * 
     * double A, dm, m, trFA_in, trFA_out, detA ; int i, j, rowLen, colLen, k_i,
     * idx_i, idx_j ; double eta, I, const_norm ;
     * 
     * if ( nrhs != 4 ) mexErrMsgTxt("To few input parameters! (Mu, w, Cov)") ;
     * 
     * 
     * //Get matrix mu colLen = muValues.numCols(); rowLen =
     * muValues.numRows();;
     * 
     * if ( mxGetM(prhs[0]) != mxGetM(prhs[2]) )
     * mexErrMsgTxt("rows of Mu and Cov not equal!") ; if ( mxGetN(prhs[0]) !=
     * mxGetN(prhs[1]) ) mexErrMsgTxt("columns of Mu and Cov not equal!") ; if (
     * mxGetN(prhs[0]) != mxGetN(prhs[2]) )
     * mexErrMsgTxt("columns of Mu and w not equal!") ;
     * 
     * 
     * 
     * 
     * const_norm = Math.pow(1.0/(2.0*Math.PI), (((double)rowLen)/2.0)) ; I =
     * 0.0 ; for( i = 0 ; i < colLen ; i++ ) { idx_i = i*rowLen ;
     * 
     * for ( j = i ; j < colLen ; j++ ) { idx_j = j*rowLen ; if ( i == j ) { eta
     * = 1.0 ; } else { eta = 2.0 ; }
     * 
     * m = 0.0 ; trFA_in = 0.0 ; trFA_out = 0.0 ; detA = 1.0 ; for ( k_i = 0 ;
     * k_i < rowLen ; ++k_i ) { A = 1.0 / (gValue[0] + cValues[idx_i + k_i] +
     * cValues[idx_j + k_i]) ;
     * 
     * dm = muValues[idx_i + k_i] - muValues[idx_j + k_i] ; m += dm*dm*A ;
     * 
     * trFA_in += A*A ; trFA_out += A ; detA *= A ; } trFA_out *= trFA_out ;
     * 
     * 
     * I += wValues[i] * wValues[j] * eta * (
     * const_norm*Math.sqrt(detA)*Math.exp(-0.5*m) ) * ( 2*trFA_in*(1.0-2.0*m) +
     * (1.0-m)*(1.0-m)*trFA_out ) ;
     * 
     * 
     * } } }
     */
}
