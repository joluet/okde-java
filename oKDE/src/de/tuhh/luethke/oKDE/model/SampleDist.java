package de.tuhh.luethke.oKDE.model;

import java.util.ArrayList;
import java.util.List;

import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;
import org.ejml.simple.SimpleMatrix;
import org.ejml.simple.SimpleSVD;

import de.tuhh.luethke.oKDE.Exceptions.EmptyDistributionException;
import de.tuhh.luethke.oKDE.utility.MomentMatcher;

public class SampleDist {
	
	//TODO: remove this field
	public double compressionError;

	// bandwidth factor
	private double mBandwidthFactor;

	// bandwidth matrix
	private SimpleMatrix mBandwidthMatrix;

	// overall weight sum
	private double mWeightSum;

	// overall mean
	private SimpleMatrix mGlobalMean;

	// overall covariance
	private SimpleMatrix mGlobalCovariance;

	// overall covariance plus bandwidth
	private SimpleMatrix mGlobalCovarianceSmoothed;

	// overall covariance in subspace
	private SimpleMatrix mSubspaceGlobalCovariance;

	// overall inverse covariance in subspace
	private SimpleMatrix mSubspaceInverseCovariance;

	// subspace: row/column ids
	private ArrayList<Integer> mSubspace;

	// component distributions
	private ArrayList<SampleDist> mSubDistributions;

	public ArrayList<SampleDist> getSubDistributions() {
		return mSubDistributions;
	}

	private double mForgettingFactor;

	// effective number of obesverved samples
	private double mN_eff;

	/**
	 * Create new SampleDist
	 */
	public SampleDist() {
		super();
		mWeightSum = 0;
		mN_eff = 0;
		mForgettingFactor = 1;
		mSubDistributions = new ArrayList<SampleDist>();
	}

	public SampleDist(double[] weights, SimpleMatrix[] means, SimpleMatrix[] covariances) {
		super();
		mSubDistributions = new ArrayList<SampleDist>();
		for(int i=0; i< means.length; i++)
			mSubDistributions.add(new SampleDist(weights[i], means[i], covariances[i]));
		mWeightSum = 0;
		for (double w : weights) {
			mWeightSum += w;
		}
		mN_eff = weights.length;
		mForgettingFactor = 1;
	}

	public SampleDist(double w, SimpleMatrix mean, SimpleMatrix covariance) {
		super();
		mWeightSum = w;
		mGlobalMean = mean;
		mGlobalCovariance = covariance;
		// initialize smoothed covariance to zero
		mGlobalCovarianceSmoothed = new SimpleMatrix(covariance).scale(0);
		// initialize bandwidth matrix to zero
		mBandwidthMatrix = new SimpleMatrix(covariance).scale(0);
		mN_eff = 0;
		mForgettingFactor = 1;
		mSubDistributions = new ArrayList<SampleDist>();
	}

	/**
	 * Copy constructor
	 * 
	 * @param dist
	 */
	public SampleDist(SampleDist dist) {
		List<SampleDist> subDists = dist.getSubDistributions();
		ArrayList<SampleDist> copy = new ArrayList<SampleDist>();
		for (SampleDist d : subDists) {
			copy.add(new SampleDist(d));
		}
		this.mSubDistributions = copy;
		this.mGlobalCovarianceSmoothed = dist.getmGlobalCovarianceSmoothed();
		this.mBandwidthFactor = dist.getmBandwidthFactor();
		this.mBandwidthMatrix = dist.getmBandwidthMatrix();
		this.mGlobalCovariance = dist.getGlobalCovariance();
		this.mGlobalMean = dist.getGlobalMean();
		this.mSubspace = dist.getmSubspace();
		this.mSubspaceGlobalCovariance = dist.getSubspaceGlobalCovariance();
		this.mSubspaceInverseCovariance = dist.getmSubspaceInverseCovariance();
		this.mWeightSum = dist.getWeightSum();
		this.mN_eff = dist.getNeff();
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
	public void updateDistribution(SimpleMatrix[] means, SimpleMatrix[] covariances, double[] doubles) throws EmptyDistributionException {
		// at first check input parameters!
		checkInputParams(means, covariances, doubles);

		// augment distribution
		addDistributions(doubles, means, covariances);

		List<SampleDist> subDists = getSubDistributions();
		Double[] weights = new Double[subDists.size()];
		for (int i = 0; i < subDists.size(); i++)
			weights[i] = subDists.get(i).getWeightSum();

		SampleDist subSpaceDist = projectToSubspace(this);
		// reestimate bandwidth as explained in oKDE paper

		double bandwidth = reestimateBandwidth(subSpaceDist.getSubMeans().toArray(new SimpleMatrix[0]),
				subSpaceDist.getSubCovariances().toArray(new SimpleMatrix[0]), weights, subSpaceDist.getSubspaceGlobalCovariance(), mN_eff);
		System.out.println("BANDW" + bandwidth);
		subSpaceDist.setmBandwidthFactor(bandwidth);
		// project Bandwidth into original space
		SimpleMatrix bandwidthMatrix = projectBandwidthToOriginalSpace(subSpaceDist);
		this.mBandwidthMatrix = bandwidthMatrix;
		for (int i = 0; i < this.getSubDistributions().size(); i++) {
			this.getSubDistributions().get(i).setmBandwidthMatrix(bandwidthMatrix);
		}
		System.out.println("BW: " + bandwidthMatrix);
		System.out.println(bandwidthMatrix.get(0, 0) + " " + bandwidthMatrix.get(1, 1));
		if (mGlobalCovariance == null) {
			mGlobalCovariance = new SimpleMatrix(2, 2);
			System.out.println("globcov null");
		}
		/*
		 * mGlobalCovariance = mGlobalCovariance.plus(mBandwidthMatrix);
		 * System.out.println("GLOBCOV"+mGlobalCovariance);
		 */
	}

	private void checkInputParams(SimpleMatrix[] means, SimpleMatrix[] covariances, double[] weights) throws EmptyDistributionException {
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
	private void addDistributions(double[] weights, SimpleMatrix[] means, SimpleMatrix[] covariances) {
		double sumOfNewWeights = 0;
		for (int i = 0; i < weights.length; i++) {
			sumOfNewWeights += weights[i];
			mSubDistributions.add(new SampleDist(weights[i], means[i], covariances[i]));
		}

		mN_eff = mN_eff * mForgettingFactor + weights.length;

		// calculate mixing weights for old and new weights
		double mixWeightOld = mWeightSum / (mWeightSum * mForgettingFactor + sumOfNewWeights);
		double mixWeightNew = sumOfNewWeights / (mWeightSum * mForgettingFactor + sumOfNewWeights);

		mWeightSum = mWeightSum * mForgettingFactor + sumOfNewWeights;

		for (int i = 0; i < mSubDistributions.size() - weights.length; i++) {
			double tmpWeight = mSubDistributions.get(i).getWeightSum();
			mSubDistributions.get(i).setWeightSum(tmpWeight * mixWeightOld);
		}
		for (int i = mSubDistributions.size() - weights.length; i < mSubDistributions.size(); i++) {
			double tmpWeight = mSubDistributions.get(i).getWeightSum();
			mSubDistributions.get(i).setWeightSum(tmpWeight * mixWeightNew * (1d / weights.length));
		}

		// system.out.println(mixWeightOld + "-" + mixWeightNew + " " +
		// mWeightSum
		// + " " + mWeights);
	}

	public double getWeightSum() {
		return mWeightSum;
	}

	public static SimpleMatrix projectBandwidthToOriginalSpace(SampleDist distribution) {
		SimpleMatrix bandwidth = SimpleMatrix.identity(distribution.getGlobalCovariance().numCols());
		// SimpleMatrix distribution
		SimpleMatrix subSpaceBandwidth = distribution.getSubspaceGlobalCovariance().scale(Math.pow(distribution.getmBandwidthFactor(), 2));
		ArrayList<Integer> subspace = distribution.getmSubspace();
		for (int i = 0; i < subSpaceBandwidth.numRows(); i++) {
			for (int j = 0; j < subSpaceBandwidth.numCols(); j++) {
				if (subspace.contains(new Integer(i)) && subspace.contains(new Integer(j)))
					bandwidth.set(i, j, subSpaceBandwidth.get(i, j));
			}
		}
		// H(valid, valid) = H2 ;
		// H2t = iF'*H*iF ;,
		SimpleMatrix invSubspaceCov = distribution.getmSubspaceInverseCovariance();
		bandwidth = invSubspaceCov.transpose().mult(bandwidth).mult(invSubspaceCov);
		// //system.out.println(bandwidth.get(1,1));
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
		if (this.mGlobalCovariance == null)
			this.mGlobalCovarianceSmoothed = new SimpleMatrix(mBandwidthMatrix);
		else
			this.mGlobalCovarianceSmoothed = this.mGlobalCovariance.plus(mBandwidthMatrix);
		this.mBandwidthMatrix = mBandwidthMatrix;
	}

	public static SampleDist projectToSubspace(SampleDist dist) throws EmptyDistributionException {
		double minBW = 1e-7;
		SampleDist distribution = new SampleDist(dist);
		ArrayList<Integer> subSpace = new ArrayList<Integer>();
		MomentMatcher.matchMoments(distribution, false);
		// system.out.println(subSpaceDist.getMeans().get(0));
		SimpleMatrix overallCovariance = distribution.getGlobalCovariance();
		// system.out.println("cov: " + overallCovariance);
		SimpleSVD svd = overallCovariance.svd(true);
		SimpleMatrix U = svd.getU();
		SimpleMatrix S = svd.getW();
		SimpleMatrix V = svd.getV();
		// system.out.println("u" + U);
		// system.out.println("s" + S);
		// system.out.println("v" + V);
		// system.out.println("testSVD: " + U.mult(S).mult(V.transpose()));
		S = S.extractDiag();

		SimpleMatrix F = new SimpleMatrix(0, 0);
		double count = 0, mean = 0;
		for (int i = 0; i < U.numRows(); i++) {
			if (S.get(i, 0) > minBW) {
				subSpace.add(i);
				SimpleMatrix colU = U.extractVector(false, i);
				double rowW = Math.pow(S.get(i, 0), -0.5);
				colU = colU.scale(rowW);
				// //system.out.println("col:" + colU);
				F = F.combine(0, F.numCols(), colU);
				mean += S.get(i, 0);
				count++;
			}
		}
		// //system.out.println("F: " + F);
		mean = (mean / count) * 1e-2;
		// //system.out.println("mean" + mean);
		for (int i = 0; i < S.numRows(); i++) {
			// //system.out.println(S.get(i, 0) + " - " + minBW);
			if (S.get(i, 0) < minBW) {
				S.set(i, 0, mean);
				// //system.out.println("jaaaa");
			}
		}
		// //system.out.println("S new: " + S.get(1, 0));
		SimpleMatrix iF = new SimpleMatrix(0, 0);
		for (int i = 0; i < U.numCols(); i++) {
			SimpleMatrix coliF = U.extractVector(false, i);
			double rowW = Math.pow(S.get(i, 0), 0.5);
			coliF = coliF.scale(rowW).transpose();
			iF = iF.combine(iF.numRows(), 0, coliF);
			// //system.out.println("col:" + coliF);
		}
		// //system.out.println(iF);
		SimpleMatrix subspaceCov = F.transpose().mult(overallCovariance).mult(F);
		distribution.setSubspaceGlobalCovariance(subspaceCov);

		ArrayList<SimpleMatrix> originalMeans = distribution.getSubMeans();
		SimpleMatrix subspaceMean = distribution.getGlobalMean();
		for (int i = 0; i < originalMeans.size(); i++) {
			originalMeans.set(i, originalMeans.get(i).minus(subspaceMean));
		}
		ArrayList<SimpleMatrix> covariances = distribution.getSubCovariances();
		for (int i = 0; i < originalMeans.size(); i++) {
			originalMeans.set(i, F.transpose().mult(originalMeans.get(i)));
			covariances.set(i, F.transpose().mult(covariances.get(i)).mult(F));
		}
		distribution.setSubCovariances(covariances);
		distribution.setSubMeans(originalMeans);

		distribution.setmSubspaceInverseCovariance(iF);
		distribution.setmSubspace(subSpace);
		// //system.out.println("Array: "+distribution.getmSubspace());
		return distribution;
	}

	public void setSubMeans(ArrayList<SimpleMatrix> means) {
		for (int i = 0; i < mSubDistributions.size(); i++) {
			mSubDistributions.get(i).setGlobalMean(means.get(i));
		}
	}

	public void setSubCovariances(ArrayList<SimpleMatrix> covariances) {
		for (int i = 0; i < mSubDistributions.size(); i++) {
			mSubDistributions.get(i).setGlobalCovariance(covariances.get(i));
		}
	}

	public ArrayList<SimpleMatrix> getSubSmoothedCovariances() {
		ArrayList<SimpleMatrix> covs = new ArrayList<SimpleMatrix>();
		for (SampleDist d : mSubDistributions)
			covs.add(d.getmGlobalCovarianceSmoothed());
		return covs;
	}

	public SimpleMatrix getmGlobalCovarianceSmoothed() {
		return mGlobalCovarianceSmoothed;
	}

	public void setmGlobalCovarianceSmoothed(SimpleMatrix mGlobalCovarianceSmoothed) {
		this.mGlobalCovarianceSmoothed = mGlobalCovarianceSmoothed;
	}

	public ArrayList<SimpleMatrix> getSubMeans() {
		ArrayList<SimpleMatrix> means = new ArrayList<SimpleMatrix>();
		for (SampleDist d : mSubDistributions)
			means.add(d.getGlobalMean());
		return means;
	}

	public ArrayList<SimpleMatrix> getSubCovariances() {
		ArrayList<SimpleMatrix> covs = new ArrayList<SimpleMatrix>();
		for (SampleDist d : mSubDistributions)
			covs.add(d.getGlobalCovariance());
		return covs;
	}

	public ArrayList<Double> getSubWeights() {
		ArrayList<Double> weights = new ArrayList<Double>();
		for (SampleDist d : mSubDistributions)
			weights.add(d.getWeightSum());
		return weights;
	}

	public void setSubWeights(ArrayList<Double> weights) {
		for (int i = 0; i < mSubDistributions.size(); i++) {
			mSubDistributions.get(i).setWeightSum(weights.get(i));
		}
	}

	public ArrayList<Integer> getmSubspace() {
		return mSubspace;
	}

	public void setmSubspace(ArrayList<Integer> mSubspace) {
		this.mSubspace = mSubspace;
	}

	public SimpleMatrix getGlobalCovariance() {
		return mGlobalCovariance;
	}

	public void setGlobalCovariance(SimpleMatrix globalCovariance) {
		this.mGlobalCovariance = globalCovariance;
	}

	public SimpleMatrix getSubspaceGlobalCovariance() {
		return mSubspaceGlobalCovariance;
	}

	public void setWeightSum(double weightSum) {
		this.mWeightSum = weightSum;
	}

	public void setSubspaceGlobalCovariance(SimpleMatrix mSubspaceCovariance) {
		this.mSubspaceGlobalCovariance = mSubspaceCovariance;
	}

	public SimpleMatrix getmSubspaceInverseCovariance() {
		return mSubspaceInverseCovariance;
	}

	public void setmSubspaceInverseCovariance(SimpleMatrix mSubspaceInverseCovariance) {
		this.mSubspaceInverseCovariance = mSubspaceInverseCovariance;
	}

	private double reestimateBandwidth(SimpleMatrix[] means, SimpleMatrix[] covariance, Double[] weights, SimpleMatrix Cov_smp, double N_eff) {

		double d = means[0].numRows();

		// Silverman
		// SimpleMatrix G = Cov_smp.scale(Math.pow((4 / ((d + 2) * N_eff)), (2 /
		// (d + 4))));

		// other
		// Cov_smp *(2/(2+d))^(2/(4+d)) * 4 *N_eff^(-2/(4+d))
		SimpleMatrix G = Cov_smp.scale(Math.pow((2d / (d + 2d)), (2d / (d + 4d))) * 4 * Math.pow(N_eff, -2d / (4d + d)));

		float alphaScale = 1;
		SimpleMatrix F = Cov_smp.scale(alphaScale);

		double Rf2 = getIntSquaredHessian(means, weights, covariance, F, G);
		double hAmise = Math.pow((Math.pow(N_eff, (-1)) * Math.pow(F.determinant(), (-1 / 2)) / (Math.pow(Math.sqrt(4 * Math.PI), d) * Rf2 * d)),
				(1 / (d + 4)));
		// system.out.println("hAmise: " + hAmise);
		return hAmise;
	}

	private double getIntSquaredHessian(SimpleMatrix[] means, Double[] weights, SimpleMatrix[] covariance, SimpleMatrix F, SimpleMatrix g) {
		long d = means[0].numRows();
		long N = means.length;
		// system.out.println("d:" + d);
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
					f_t = constNorm * Math.sqrt(A.determinant()) * Math.exp(-0.5 * ds.mult(dm).trace());
					c = 2 * F.mult(A).mult(F).mult(B).trace() + Math.pow(F.mult(C).trace(), 2);
				} else {
					m = dm.transpose().mult(A).mult(dm).get(0);
					// //system.out.println("m: " + m);
					// //system.out.println("A:" + A);
					f_t = constNorm * Math.sqrt(A.determinant()) * Math.exp(-0.5 * m);

					DenseMatrix64F A_sqr = new DenseMatrix64F(A.numRows(), A.numCols());
					CommonOps.elementMult(A.getMatrix(), A.transpose().getMatrix(), A_sqr);
					double sum = CommonOps.elementSum(A_sqr);
					c = 2d * sum * (1d - 2d * m) + Math.pow((1d - m), 2d) * Math.pow(A.trace(), 2);
				}

				// determine the weight of the current term
				if (i1 == i2)
					eta = 1;
				else
					eta = 2;
				// //system.out.println(" f_t: " + f_t + " c: " + c
				// + " w2*w1:" + (w2 * w1));
				I = I + f_t * c * w2 * w1 * eta;
				// //system.out.println("I: " + I);
			}
		}

		return I;
	}

	public double evaluate(SimpleMatrix pointVector, boolean useSubDists, boolean useSmoothedCov) {
		ArrayList<SimpleMatrix> means = new ArrayList<SimpleMatrix>();
		ArrayList<SimpleMatrix> covs = new ArrayList<SimpleMatrix>();
		ArrayList<Double> weights = new ArrayList<Double>();
		if(useSubDists) {
			means = this.getSubMeans();
			if(useSmoothedCov)
				covs = this.getSubSmoothedCovariances();
			else
				covs = this.getSubCovariances();
			weights = this.getSubWeights();
		}else{
			means.add(this.mGlobalMean);
			covs.add(this.mGlobalCovariance);
			weights.add(this.mWeightSum);
		}
			
		SimpleMatrix bandwidth = this.mBandwidthMatrix;
		/*double[][] dxVector = { { x }, { y } };
		SimpleMatrix xVector = new SimpleMatrix(dxVector);*/
		double d = 0d;
		double n = means.get(0).numRows();
		double a = Math.pow(Math.sqrt(2*Math.PI),n);
		for (int i=0; i< means.size(); i++) {
			SimpleMatrix m = means.get(i);
			SimpleMatrix c = covs.get(i);
			double w = weights.get(i);
			double tmp = (-0.5d) * pointVector.minus(m).transpose().mult(c.invert()).mult(pointVector.minus(m)).trace();
			d += ((1 / (a*Math.sqrt(c.determinant()) )) * Math.exp(tmp)) * w;
		}
		return d;
	}
	
	/**
	 * Evaluates the distribution at the given n-dimensional points and returns the results in a List of double-values.
	 * @param points 
	 * @return array of double values
	 */
	public ArrayList<Double> evaluate(ArrayList<SimpleMatrix> points, boolean useSubDists, boolean useSmoothedCov){
		ArrayList<Double> resultPoints = new ArrayList<Double>();
		for(SimpleMatrix point : points) {
			resultPoints.add( evaluate(point, useSubDists, useSmoothedCov) );
		}
		return resultPoints;
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

	public SimpleMatrix getGlobalMean() {
		return mGlobalMean;
	}

	public void setGlobalMean(SimpleMatrix globalMean) {
		this.mGlobalMean = globalMean;
	}

	public double getNeff() {
		return this.mN_eff;
	}

	public void setGlobalMean(double nEff) {
		this.mN_eff = nEff;
	}

	public void setSubDistributions(ArrayList<SampleDist> subDistributions) {
		this.mSubDistributions = subDistributions;
	}

}
