package de.tuhh.luethke.oKDE.model;

import java.util.ArrayList;

import org.ejml.data.DenseMatrix64F;
import org.ejml.factory.SingularValueDecomposition;
import org.ejml.ops.CommonOps;
import org.ejml.simple.SimpleMatrix;
import org.ejml.simple.SimpleSVD;

import de.tuhh.luethke.oKDE.Exceptions.EmptyDistributionException;
import de.tuhh.luethke.oKDE.utility.MomentMatcher;

public class SampleDist {

    // bandwidth factor
    private double mBandwidthFactor;

    // bandwidth matrix
    private SimpleMatrix mBandwidthMatrix;

    // component means
    private ArrayList<SimpleMatrix> mMeans;

    // component covariances
    private ArrayList<SimpleMatrix> mCovariances;

    // component weights
    private ArrayList<Double> mWeights;

    // overall covariance
    private SimpleMatrix mCovariance;

    // overall covariance in subspace
    private SimpleMatrix mSubspaceCovariance;

    // overall inverse covariance in subspace
    private SimpleMatrix mSubspaceInverseCovariance;

    // subspace: row/column ids
    private ArrayList<Integer> mSubspace;

    private double mWeightSum;

    private double mForgettingFactor;

    // effective number of obesverved samples
    private double mN_eff;

    /**
     * Create new SampleDist
     */
    public SampleDist() {
	super();
	mMeans = new ArrayList<SimpleMatrix>();
	mCovariances = new ArrayList<SimpleMatrix>();
	mWeightSum = 0;
	mWeights = new ArrayList<Double>();
	mN_eff = 0;
	mForgettingFactor = 1;
    }

    public SampleDist(ArrayList<Double> weights, ArrayList<SimpleMatrix> means,
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
    }

    /**
     * Returns the component weights of the distribution.
     * 
     * @return weights component weights
     */
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
     * @param doubles
     * @throws EmptyDistributionException
     *             Exception is thrown when one parameter is null or empty
     */
    public void updateDistribution(SimpleMatrix[] means,
	    SimpleMatrix[] covariances, double[] doubles)
	    throws EmptyDistributionException {
	// at first check input parameters!
	checkInputParams(means, covariances, doubles);

	// augment distribution
	addWeights(doubles);
	for (SimpleMatrix m : means)
	    mMeans.add(m);
	for (SimpleMatrix c : covariances)
	    mCovariances.add(c);
	
	Double[] weights = getWeights().toArray(new Double[0]);

	// get empirical sample covariance by moment matching
	SimpleMatrix covSmpl = projectToSubspace(this);
	// reestimate bandwidth as explained in oKDE paper
	double bandwidth = reestimateBandwidth(this.getMeans().toArray(new SimpleMatrix[0]), this.getCovariances().toArray(new SimpleMatrix[0]), weights,
		covSmpl, mN_eff);
	this.setmBandwidthFactor(bandwidth);
	// project Bandwidth into original space
	SimpleMatrix bandwidthMatrix = projectBandwidthToOriginalSpace(this);
	this.setmBandwidthMatrix(bandwidthMatrix);
	System.out.println("BW: "+bandwidthMatrix);
	System.out.println(bandwidthMatrix.get(0,0)+" "+bandwidthMatrix.get(1,1));
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
	    mWeights.set(i, mWeights.get(i) * mixWeightNew * (1d/weights.length));
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

    public static SimpleMatrix projectBandwidthToOriginalSpace(
	    SampleDist distribution) {
	SimpleMatrix bandwidth = SimpleMatrix.identity(distribution
		.getmCovariance().numCols());
	// SimpleMatrix distribution
	SimpleMatrix subSpaceBandwidth = distribution.getmSubspaceCovariance()
		.scale(Math.pow(distribution.getmBandwidthFactor(), 2));
	ArrayList<Integer> subspace = distribution.getmSubspace();
	for (int i = 0; i < subSpaceBandwidth.numRows(); i++) {
	    for (int j = 0; j < subSpaceBandwidth.numCols(); j++) {
		if(subspace.contains(new Integer(i)) && subspace.contains(new Integer(j)))
		bandwidth.set(i, j,
			subSpaceBandwidth.get(i, j));
	    }
	}
	// H(valid, valid) = H2 ;
	// H2t = iF'*H*iF ;,
	SimpleMatrix invSubspaceCov = distribution.getmSubspaceInverseCovariance();
	bandwidth = invSubspaceCov.transpose().mult(bandwidth).mult(invSubspaceCov);
	//System.out.println(bandwidth.get(1,1));
	return bandwidth;
    }

    public double getmBandwidthFactor() {
	return mBandwidthFactor;
    }

    public void setmBandwidthFactor(double mBandwidthFactor) {
	this.mBandwidthFactor = mBandwidthFactor;
    }

    public SimpleMatrix getmBandwidthMatrix() {
	return mBandwidthMatrix;
    }

    public void setmBandwidthMatrix(SimpleMatrix mBandwidthMatrix) {
	this.mBandwidthMatrix = mBandwidthMatrix;
    }

    public static SimpleMatrix projectToSubspace(SampleDist distribution)
	    throws EmptyDistributionException {
	double minBW = 1e-7;
	ArrayList<Integer> subSpace = new ArrayList<Integer>();
	SampleDist subSpaceDist = MomentMatcher.matchMoments(distribution);
	System.out.println(subSpaceDist.getMeans().get(0));
	SimpleMatrix overallCovariance = subSpaceDist.getCovariances().get(0);
	System.out.println("cov: " + overallCovariance);
	distribution.setmCovariance(overallCovariance);
	SimpleSVD svd = overallCovariance.svd(true);
	SimpleMatrix U = svd.getU();
	SimpleMatrix S = svd.getW();
	SimpleMatrix V = svd.getV();
	System.out.println("u" + U);
	System.out.println("s" + S);
	System.out.println("v" + V);
	System.out.println("testSVD: " + U.mult(S).mult(V.transpose()));
	S = S.extractDiag();

	SimpleMatrix F = new SimpleMatrix(0, 0);
	double count = 0, mean = 0;
	for (int i = 0; i < U.numRows(); i++) {
	    if (S.get(i, 0) > minBW) {
		subSpace.add(i);
		SimpleMatrix colU = U.extractVector(false, i);
		double rowW = Math.pow(S.get(i, 0), -0.5);
		colU = colU.scale(rowW);
		//System.out.println("col:" + colU);
		F = F.combine(0, F.numCols(), colU);
		mean += S.get(i, 0);
		count++;
	    }
	}
	//System.out.println("F: " + F);
	mean = (mean / count) * 1e-2;
	//System.out.println("mean" + mean);
	for (int i = 0; i < S.numRows(); i++) {
	    //System.out.println(S.get(i, 0) + " - " + minBW);
	    if (S.get(i, 0) < minBW) {
		S.set(i, 0, mean);
		//System.out.println("jaaaa");
	    }
	}
	//System.out.println("S new: " + S.get(1, 0));
	SimpleMatrix iF = new SimpleMatrix(0, 0);
	for (int i = 0; i < U.numCols(); i++) {
	    SimpleMatrix coliF = U.extractVector(false, i);
	    double rowW = Math.pow(S.get(i, 0), 0.5);
	    coliF = coliF.scale(rowW).transpose();
	    iF = iF.combine(iF.numRows(), 0, coliF);
	    //System.out.println("col:" + coliF);
	}
	//System.out.println(iF);
	SimpleMatrix subspaceCov = F.transpose().mult(overallCovariance).mult(F);
	distribution.setmSubspaceCovariance(subspaceCov);
	
	ArrayList<SimpleMatrix> originalMeans = distribution.getMeans();
	SimpleMatrix subspaceMean = subSpaceDist.getMeans().get(0);
	for(int i=0; i<originalMeans.size(); i++){
	    originalMeans.set(i, originalMeans.get(i).minus(subspaceMean));
	}
	ArrayList<SimpleMatrix> covariances = distribution.getCovariances();
	for(int i=0; i<originalMeans.size(); i++){
	    originalMeans.set(i, F.transpose().mult(originalMeans.get(i)));
	    covariances.set(i, F.transpose().mult(covariances.get(i)).mult(F));
	}
	distribution.setCovariances(covariances);
	distribution.setMeans(originalMeans);
	
	distribution.setmSubspaceInverseCovariance(iF);
	distribution.setmSubspace(subSpace);
	//System.out.println("Array: "+distribution.getmSubspace());

	return subspaceCov;
    }

    public ArrayList<Integer> getmSubspace() {
	return mSubspace;
    }

    public void setmSubspace(ArrayList<Integer> mSubspace) {
	this.mSubspace = mSubspace;
    }

    public SimpleMatrix getmCovariance() {
	return mCovariance;
    }

    public void setmCovariance(SimpleMatrix mCovariance) {
	this.mCovariance = mCovariance;
    }

    public SimpleMatrix getmSubspaceCovariance() {
	return mSubspaceCovariance;
    }

    public void setmSubspaceCovariance(SimpleMatrix mSubspaceCovariance) {
	this.mSubspaceCovariance = mSubspaceCovariance;
    }

    public SimpleMatrix getmSubspaceInverseCovariance() {
	return mSubspaceInverseCovariance;
    }

    public void setmSubspaceInverseCovariance(
	    SimpleMatrix mSubspaceInverseCovariance) {
	this.mSubspaceInverseCovariance = mSubspaceInverseCovariance;
    }

    private double reestimateBandwidth(SimpleMatrix[] means,
	    SimpleMatrix[] covariance, Double[] weights, SimpleMatrix Cov_smp,
	    double N_eff) {

	double d = means[0].numRows();

	//Silverman
	SimpleMatrix G = Cov_smp.scale(Math.pow((4 / ((d + 2) * N_eff)),
		(2 / (d + 4))));
	
	//other
	// Cov_smp *(2/(2+d))^(2/(4+d)) * 4 *N_eff^(-2/(4+d))
	//SimpleMatrix G = Cov_smp.scale(Math.pow((2d / (d + 2d) ),
	//		(2d / (d + 4d)) ) * 4 * Math.pow(N_eff,-2d/(4d+d)) );
	
	float alphaScale = 1;
	SimpleMatrix F = Cov_smp.scale(alphaScale);

	double Rf2 = getIntSquaredHessian(means, weights, covariance, F, G);
	double hAmise = Math
		.pow((Math.pow(N_eff, (-1))
			* Math.pow(F.determinant(), (-1 / 2)) / (Math.pow(
			Math.sqrt(4 * Math.PI), d)
			* Rf2 * d)), (1 / (d + 4)));
	System.out.println("hAmise: " + hAmise);
	return hAmise;
    }

    private double getIntSquaredHessian(SimpleMatrix[] means, Double[] weights,
	    SimpleMatrix[] covariance, SimpleMatrix F, SimpleMatrix g) {

	long d = means[0].numRows();
	long N = means.length;
	System.out.println("d:" + d);
	// normalizer
	double constNorm = Math.pow((1d / (2d * Math.PI)), (d / 2d));

	// test if F is identity for speedup
	SimpleMatrix Id = SimpleMatrix.identity(F.numCols());
	double deltaF = F.minus(Id).elementSum();
	
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
		    //System.out.println("m: " + m);
		    //System.out.println("A:" + A);
		    f_t = constNorm * Math.sqrt(A.determinant())
			    * Math.exp(-0.5 * m);

		    DenseMatrix64F A_sqr = new DenseMatrix64F(A.numRows(),
			    A.numCols());
		    CommonOps.elementMult(A.getMatrix(), A.transpose()
			    .getMatrix(), A_sqr);
		    double sum = CommonOps.elementSum(A_sqr);
		    c = 2d * sum * (1d - 2d * m) + Math.pow((1d - m), 2d)
			    * Math.pow(A.trace(), 2);
		}

		// determine the weight of the current term
		if (i1 == i2)
		    eta = 1;
		else
		    eta = 2;
		//System.out.println(" f_t: " + f_t + " c: " + c
		//	+ " w2*w1:" + (w2 * w1));
		I = I + f_t * c * w2 * w1 * eta;
		//System.out.println("I: " + I);
	    }
	}

	return I;
    }
    
    public double calculateY(double[] x, SampleDist dist,
	    ArrayList<SimpleMatrix> means) {
	SimpleMatrix bandwidth = dist.getmBandwidthMatrix();
	double[][] dxVector = { x };
	SimpleMatrix xVector = new SimpleMatrix(dxVector);
	double d = 0d;
	double n = means.size();
	for (SimpleMatrix m : means) {
	    double tmp = (-0.5d)
		    * xVector.minus(m).transpose().mult(bandwidth.invert())
			    .mult(xVector.minus(m)).trace();
	    double powPi = Math.pow(4*Math.PI, xVector.numRows());
	    d += ((1 / Math.sqrt(powPi
		    * bandwidth.determinant())) * Math.exp(tmp));
	}
	return d / n;
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

    public void setMeans(ArrayList<SimpleMatrix> means) {
	this.mMeans = means;
    }

    public void setCovariances(ArrayList<SimpleMatrix> covariances) {
	this.mCovariances = covariances;
    }

    public void setWeights(ArrayList<Double> weights) {
	this.mWeights = weights;
    }
}
