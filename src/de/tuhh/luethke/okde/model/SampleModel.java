/*
 * Copyright 2014 Jonas Luethke
 */


package de.tuhh.luethke.okde.model;

import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;
import org.ejml.simple.SimpleEVD;
import org.ejml.simple.SimpleMatrix;
import org.ejml.simple.SimpleSVD;

import de.tuhh.luethke.okde.Exceptions.EmptyDistributionException;
import de.tuhh.luethke.okde.utility.MomentMatcher;
import de.tuhh.luethke.okde.utility.Compression.Compressor;
import de.tuhh.luethke.okde.utility.Matrices.HashableSimpleMatrix;
import de.tuhh.luethke.okde.utility.Optimization.SearchResult;

public class SampleModel extends BaseSampleDistribution {

	private static final float DEFAULT_NO_OF_COMPS_THRES = 6;
	
	private HashMap<HashableSimpleMatrix, Double> mProbabilityCache = new HashMap<HashableSimpleMatrix, Double>();
	
	// When mahalanobis distance gets bigger than this the component does not contribute
	// to the density that should be calculated
	// exp(-40) ~ 5E-18
	private static final double MAX_MAHALANOBIS_DIST = 40;

	// compression threshold (maximal hellinger distance)
	public double mCompressionThreshold;
	
	// effective number of observed samples
	protected double mEffectiveNoOfSamples;

	// component distributions
	protected ArrayList<BaseSampleDistribution> mSubDistributions;

	// threshold to determine when compression is necessary
	public float mNoOfCompsThreshold;
	
	public double mEMError=0;
	public double mEMCount=0;
	
	public long bwtime = 0;
	
	public SampleModel(double forgettingFactor, double compressionThreshold) {
		super();
		this.mSubDistributions = new ArrayList<BaseSampleDistribution>();
		this.mBandwidthMatrix = null;
		this.mGlobalCovariance = null;
		this.mGlobalMean = null;
		this.mSubspace = null;
		this.mSubspaceGlobalCovariance = null;
		this.mSubspaceInverseCovariance = null;
		this.mGlobalWeight = 0;
		this.mEffectiveNoOfSamples = 0;
		this.mForgettingFactor = forgettingFactor;
		this.mCompressionThreshold = compressionThreshold;
		mNoOfCompsThreshold = DEFAULT_NO_OF_COMPS_THRES;
	}

	/**
	 * Copy constructor
	 * 
	 * @param dist
	 * @throws IllegalAccessException
	 * @throws InstantiationException
	 * @throws SecurityException
	 * @throws NoSuchMethodException
	 * @throws InvocationTargetException
	 * @throws IllegalArgumentException
	 */
	public SampleModel(SampleModel dist) throws IllegalArgumentException, InvocationTargetException, NoSuchMethodException, SecurityException,
			InstantiationException, IllegalAccessException {
		List<BaseSampleDistribution> subDists = dist.getSubDistributions();
		ArrayList<BaseSampleDistribution> copy = new ArrayList<BaseSampleDistribution>();
		// copy the list of sub components
		for (BaseSampleDistribution d : subDists) {
			// We don't know if a sub component is of type
			// OneComponentDistribution or TwoComponentDistribution.
			// Thus we have to call the constructor using reflection api. This
			// way it is generic and works in both cases.
			Constructor<? extends BaseSampleDistribution> ctor = d.getClass().getDeclaredConstructor(d.getClass());
			ctor.setAccessible(true);
			BaseSampleDistribution tmp = (BaseSampleDistribution) ctor.newInstance(d);
			copy.add(tmp);
		}
		this.mSubDistributions = copy;
		this.mBandwidthMatrix = dist.getBandwidthMatrix();
		this.mGlobalCovariance = dist.getGlobalCovariance();
		this.mGlobalMean = dist.getGlobalMean();
		this.mSubspace = dist.getmSubspace();
		this.mSubspaceGlobalCovariance = dist.getSubspaceGlobalCovariance();
		this.mSubspaceInverseCovariance = dist.getSubspaceInverseCovariance();
		this.mGlobalWeight = dist.getGlobalWeight();
		this.mEffectiveNoOfSamples = dist.mEffectiveNoOfSamples;
	}

	public void overWirite(SampleModel dist) throws IllegalArgumentException, InvocationTargetException, NoSuchMethodException, SecurityException,
			InstantiationException, IllegalAccessException {
		List<BaseSampleDistribution> subDists = dist.getSubDistributions();
		ArrayList<BaseSampleDistribution> copy = new ArrayList<BaseSampleDistribution>();
		// copy the list of sub components
		for (BaseSampleDistribution d : subDists) {
			// We don't know if a sub component is of type
			// OneComponentDistribution or TwoComponentDistribution.
			// Thus we have to call the constructor using reflection api. This
			// way it is generic and works in both cases.
			Constructor<? extends BaseSampleDistribution> ctor = d.getClass().getDeclaredConstructor(d.getClass());
			ctor.setAccessible(true);
			BaseSampleDistribution tmp = (BaseSampleDistribution) ctor.newInstance(d);
			copy.add(tmp);
		}
		this.mSubDistributions = copy;
		this.mBandwidthMatrix = dist.getBandwidthMatrix();
		this.mGlobalCovariance = dist.getGlobalCovariance();
		this.mGlobalMean = dist.getGlobalMean();
		this.mSubspace = dist.getmSubspace();
		this.mSubspaceGlobalCovariance = dist.getSubspaceGlobalCovariance();
		this.mSubspaceInverseCovariance = dist.getSubspaceInverseCovariance();
		this.mGlobalWeight = dist.getGlobalWeight();
		this.mEffectiveNoOfSamples = dist.mEffectiveNoOfSamples;
	}

	// subspace: row/column ids
	private ArrayList<Integer> mSubspace;

	public ArrayList<Integer> getmSubspace() {
		return mSubspace;
	}

	public void setmSubspace(ArrayList<Integer> mSubspace) {
		this.mSubspace = mSubspace;
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
	 * @throws IllegalAccessException
	 * @throws InstantiationException
	 * @throws SecurityException
	 * @throws NoSuchMethodException
	 * @throws InvocationTargetException
	 * @throws IllegalArgumentException
	 */
	public void updateDistribution(SimpleMatrix[] means, SimpleMatrix[] covariances, double[] weights) throws EmptyDistributionException,
			IllegalArgumentException, InvocationTargetException, NoSuchMethodException, SecurityException, InstantiationException,
			IllegalAccessException {
		// at first check input parameters!
		checkInputParams(means, covariances, weights);

		// augment distribution
		addDistributions(weights, means, covariances);
		
		// save indices of new components for em updates
		ArrayList<Integer> newPoints = new ArrayList<Integer>();
		for(int i=(mSubDistributions.size()-means.length); i<mSubDistributions.size(); i++) {
			newPoints.add(i);
		}
		updateDistribution(newPoints);
	}

	public void updateDistribution(SimpleMatrix mean, SimpleMatrix covariance, double weight) throws EmptyDistributionException,
			IllegalArgumentException, InvocationTargetException, NoSuchMethodException, SecurityException, InstantiationException,
			IllegalAccessException {

		// augment distribution
		addDistribution(weight, mean, covariance);
		ArrayList<Integer> newPoints = new ArrayList<Integer>();
		// only one component was added, save index for em updates
		newPoints.add(mSubDistributions.size()-1);
		updateDistribution(newPoints);
	}

	private void updateDistribution(ArrayList<Integer> newPoints) {
		List<BaseSampleDistribution> subDists = getSubDistributions();
		Double[] weights = new Double[subDists.size()];
		for (int i = 0; i < subDists.size(); i++)
			weights[i] = subDists.get(i).getGlobalWeight();

		SampleModel subSpaceDist = null;
		try {
			subSpaceDist = projectToSubspace(this);
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
		} catch (EmptyDistributionException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		// reestimate bandwidth as explained in oKDE paper
		long time = System.currentTimeMillis();
		double bandwidthFactor = reestimateBandwidth(subSpaceDist.getSubMeans().toArray(new SimpleMatrix[0]), subSpaceDist.getSubCovariances()
				.toArray(new SimpleMatrix[0]), weights, subSpaceDist.getSubspaceGlobalCovariance(), mEffectiveNoOfSamples);
		bwtime += (System.currentTimeMillis()-time);
		//System.out.println("BANDW" + bandwidthFactor);
		// project Bandwidth into original space
		SimpleMatrix bandwidthMatrix = projectBandwidthToOriginalSpace(subSpaceDist, bandwidthFactor);
		this.mBandwidthMatrix = bandwidthMatrix;
		for (int i = 0; i < this.getSubDistributions().size(); i++) {
			this.getSubDistributions().get(i).setBandwidthMatrix(bandwidthMatrix);
		}
		//System.out.println("BW: " + bandwidthMatrix);
		//System.out.println(bandwidthMatrix.get(0, 0) + " " + bandwidthMatrix.get(1, 1));
		if (mGlobalCovariance == null) {
			mGlobalCovariance = new SimpleMatrix(2, 2);
			//System.out.println("globcov null");
		}
		try {
			Compressor.compress(this, newPoints);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
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
			mSubDistributions.add(new OneComponentDistribution(weights[i], means[i], covariances[i], mBandwidthMatrix));
		}

		// calculate mixing weights for old and new weights
		double mixWeightOld = (mEffectiveNoOfSamples* mForgettingFactor) / (mEffectiveNoOfSamples * mForgettingFactor + sumOfNewWeights);
		double mixWeightNew = sumOfNewWeights / (mEffectiveNoOfSamples * mForgettingFactor + sumOfNewWeights);

		mEffectiveNoOfSamples = mEffectiveNoOfSamples * mForgettingFactor + weights.length;

		// mGlobalWeight = mGlobalWeight * mForgettingFactor + sumOfNewWeights;
		mGlobalWeight = mixWeightOld + mixWeightNew;

		for (int i = 0; i < mSubDistributions.size() - weights.length; i++) {
			double tmpWeight = mSubDistributions.get(i).getGlobalWeight();
			mSubDistributions.get(i).setGlobalWeight(tmpWeight * mixWeightOld);
		}
		for (int i = mSubDistributions.size() - weights.length; i < mSubDistributions.size(); i++) {
			double tmpWeight = mSubDistributions.get(i).getGlobalWeight();
			mSubDistributions.get(i).setGlobalWeight(tmpWeight * mixWeightNew * (1d / weights.length));
		}
		// mEffectiveNoOfSamples = mSubDistributions.size();
	}

	/**
	 * Takes new incoming sample weights and updates this distribution using a
	 * forgetting factor.
	 * 
	 * @param weights
	 */
	private void addDistribution(double weight, SimpleMatrix mean, SimpleMatrix covariance) {
		double sumOfNewWeights = 0;
		sumOfNewWeights += weight;
		mSubDistributions.add(new OneComponentDistribution(weight, mean, covariance, mBandwidthMatrix));

		// calculate mixing weights for old and new weights
		double mixWeightOld = (mEffectiveNoOfSamples* mForgettingFactor) / (mEffectiveNoOfSamples * mForgettingFactor + sumOfNewWeights);
		double mixWeightNew = sumOfNewWeights / (mEffectiveNoOfSamples * mForgettingFactor + sumOfNewWeights);

		mEffectiveNoOfSamples = mEffectiveNoOfSamples * mForgettingFactor + 1;

		// mGlobalWeight = mGlobalWeight * mForgettingFactor + sumOfNewWeights;
		mGlobalWeight = mixWeightOld + mixWeightNew;

		for (int i = 0; i < mSubDistributions.size() - 1; i++) {
			double tmpWeight = mSubDistributions.get(i).getGlobalWeight();
			mSubDistributions.get(i).setGlobalWeight(tmpWeight * mixWeightOld);
		}
		for (int i = mSubDistributions.size() - 1; i < mSubDistributions.size(); i++) {
			double tmpWeight = mSubDistributions.get(i).getGlobalWeight();
			mSubDistributions.get(i).setGlobalWeight(tmpWeight * mixWeightNew * (1d / 1));
		}
		// mEffectiveNoOfSamples = mSubDistributions.size();
	}

	private static SampleModel projectToSubspace(SampleModel dist) throws EmptyDistributionException, IllegalArgumentException,
			InvocationTargetException, NoSuchMethodException, SecurityException, InstantiationException, IllegalAccessException {
		double minBW = 1e-7;
		SampleModel distribution = new SampleModel(dist);
		ArrayList<Integer> subSpace = new ArrayList<Integer>();
		MomentMatcher.matchMoments(distribution);
		SimpleMatrix overallCovariance = distribution.getGlobalCovariance();
		//System.out.println("cov: " + overallCovariance);
		SimpleSVD<?> svd = overallCovariance.svd(true);
		SimpleMatrix U = svd.getU();
		SimpleMatrix S = svd.getW();
		S = S.extractDiag();

		SimpleMatrix F = new SimpleMatrix(0, 0);
		double count = 0, mean = 0;
		for (int i = 0; i < U.numRows(); i++) {
			if (S.get(i, 0) > minBW) {
				subSpace.add(i);
				SimpleMatrix colU = U.extractVector(false, i);
				double rowW = Math.pow(S.get(i, 0), -0.5);
				colU = colU.scale(rowW);
				F = F.combine(0, F.numCols(), colU);
				mean += S.get(i, 0);
				count++;
			}
		}
		mean = (mean / count) * 1e-2;
		for (int i = 0; i < S.numRows(); i++) {
			if (S.get(i, 0) < minBW) {
				S.set(i, 0, mean);
			}
		}
		SimpleMatrix iF = new SimpleMatrix(0, 0);
		for (int i = 0; i < U.numCols(); i++) {
			SimpleMatrix coliF = U.extractVector(false, i);
			double rowW = Math.pow(S.get(i, 0), 0.5);
			coliF = coliF.scale(rowW).transpose();
			iF = iF.combine(iF.numRows(), 0, coliF);
		}
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

		distribution.setSubspaceInverseCovariance(iF);
		distribution.setmSubspace(subSpace);
		return distribution;
	}

	private static SimpleMatrix projectBandwidthToOriginalSpace(SampleModel distribution, double bandwidthFactor) {
		SimpleMatrix bandwidth = SimpleMatrix.identity(distribution.getGlobalCovariance().numCols());
		SimpleMatrix subSpaceBandwidth = distribution.getSubspaceGlobalCovariance().scale(Math.pow(bandwidthFactor, 2));
		ArrayList<Integer> subspace = distribution.getmSubspace();
		for (int i = 0; i < subSpaceBandwidth.numRows(); i++) {
			for (int j = 0; j < subSpaceBandwidth.numCols(); j++) {
				if (subspace.contains(new Integer(i)) && subspace.contains(new Integer(j)))
					bandwidth.set(i, j, subSpaceBandwidth.get(i, j));
			}
		}
		SimpleMatrix invSubspaceCov = distribution.getSubspaceInverseCovariance();
		bandwidth = invSubspaceCov.transpose().mult(bandwidth).mult(invSubspaceCov);
		return bandwidth;
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
		return hAmise;
	}

	private double getIntSquaredHessian(SimpleMatrix[] means, Double[] weights, SimpleMatrix[] covariance, SimpleMatrix F, SimpleMatrix g) {
		long time = System.currentTimeMillis();
		long d = means[0].numRows();
		long N = means.length;
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
				I = I + f_t * c * w2 * w1 * eta;
			}
		}
		/*time = System.currentTimeMillis()-time;
		if((time) > 100)
			System.out.println("Time for IntSqrdHessian: "+ ((double)time/1000)+"s"+"  loopcount: "+N);*/
		return I;
	}

	public void setSubDistributions(ArrayList<BaseSampleDistribution> subDistributions) {
		this.mSubDistributions = subDistributions;
	}

	public ArrayList<BaseSampleDistribution> getSubDistributions() {
		return mSubDistributions;
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
		for (BaseSampleDistribution d : mSubDistributions)
			covs.add(d.getmGlobalCovarianceSmoothed());
		return covs;
	}

	public ArrayList<SimpleMatrix> getSubMeans() {
		ArrayList<SimpleMatrix> means = new ArrayList<SimpleMatrix>();
		for (BaseSampleDistribution d : mSubDistributions)
			means.add(d.getGlobalMean());
		return means;
	}

	public ArrayList<SimpleMatrix> getSubCovariances() {
		ArrayList<SimpleMatrix> covs = new ArrayList<SimpleMatrix>();
		for (BaseSampleDistribution d : mSubDistributions)
			covs.add(d.getGlobalCovariance());
		return covs;
	}

	public ArrayList<Double> getSubWeights() {
		ArrayList<Double> weights = new ArrayList<Double>();
		for (BaseSampleDistribution d : mSubDistributions)
			weights.add(d.getGlobalWeight());
		return weights;
	}

	public void setSubWeights(ArrayList<Double> weights) {
		for (int i = 0; i < mSubDistributions.size(); i++) {
			mSubDistributions.get(i).setGlobalWeight(weights.get(i));
		}
	}
	
	/**
	 * This method derives the conditional distribution of the sample model kde with distribution p(x).
	 * It marginalizes on the first n dimensions. The number n is given as a parameter.
	 * @param firstDimensions
	 * @return
	 */
	public ConditionalDistribution getMarginalDistribution(int n){
		ArrayList<SimpleMatrix> means = this.getSubMeans();
		ArrayList<SimpleMatrix> marginalMeans = new ArrayList<SimpleMatrix>();

		ArrayList<SimpleMatrix> covs = this.getSubSmoothedCovariances();
		ArrayList<SimpleMatrix> marginalCovs = new ArrayList<SimpleMatrix>();
		
		ArrayList<Double> weights = this.getSubWeights();
		ArrayList<Double> marginalWeights = new ArrayList<Double>();

		ConditionalDistribution result = null;
		
		double a = Math.pow(Math.sqrt(2 * Math.PI), n);

		
		for(int i=0; i<means.size(); i++) {
			SimpleMatrix c = covs.get(i);
			SimpleMatrix m = means.get(i);
			SimpleMatrix m1 = new SimpleMatrix(n,1);
			
			// extract all elements from covariance that correspond only to m1
			// that means extract the block in the left top corner with height=width=n
			SimpleMatrix newC1 = new SimpleMatrix(n,n);
			for(int j=0; j<n; j++) {
				for(int k=0; k<n; k++) {
					newC1.set(j, k, c.get(j,k) );
				}
			}
			//extract last rows from mean to m1
			for(int j=0; j<n; j++) {
				m1.set(j,0,m.get(j,0));
			}
			
			marginalMeans.add(m1);
			marginalCovs.add(newC1);

		}
		result = new ConditionalDistribution(marginalMeans, marginalCovs, weights);
		return result;
	}
	
	/**
	 * This method derives the conditional distribution of the actual sample model kde with distribution p(x).
	 * It takes a condition parameter that is a vector c of dimension m. Using this vector
	 * it finds the conditional distribution p(x*|c) where c=(x_0,...,x_m), x*=(x_m+1,...,x_n).
	 * For detailed description see:
	 * @param condition A vector that defines c in p(x*|c)
	 * @return The conditional distribution of this sample model under the given condition
	 */
	public ConditionalDistribution getConditionalDistribution(SimpleMatrix condition){
		int lenCond = condition.numRows();
		
		ArrayList<SimpleMatrix> means = this.getSubMeans();
		ArrayList<SimpleMatrix> conditionalMeans = new ArrayList<SimpleMatrix>();

		ArrayList<SimpleMatrix> covs = this.getSubSmoothedCovariances();
		ArrayList<SimpleMatrix> conditionalCovs = new ArrayList<SimpleMatrix>();
		
		ArrayList<Double> weights = this.getSubWeights();
		ArrayList<Double> conditionalWeights = new ArrayList<Double>();

		ConditionalDistribution result = null;
		
		double n = condition.numRows();
		double a = Math.pow(Math.sqrt(2 * Math.PI), n);

		
		for(int i=0; i<means.size(); i++) {
			SimpleMatrix c = covs.get(i);
			SimpleMatrix invC = c.invert();
			SimpleMatrix m = means.get(i);
			int lenM1 = m.numRows()-lenCond;
			SimpleMatrix m1 = new SimpleMatrix(lenM1,1);
			SimpleMatrix m2 = new SimpleMatrix(lenCond,1);
			
			// extract all elements from inverse covariance that correspond only to m1
			// that means extract the block in the right bottom corner with height=width=lenM1
			SimpleMatrix newC1 = new SimpleMatrix(lenM1,lenM1);
			for(int j=0; j<lenM1; j++) {
				for(int k=0; k<lenM1; k++) {
					newC1.set(j, k, invC.get(j+lenCond,k+lenCond) );
				}
			}
			// extract all elements from inverse covariance that correspond to m1 and m2
			// from the the block in the left bottom corner with height=width=lenM1
			SimpleMatrix newC2 = new SimpleMatrix(lenM1,lenCond);
			for(int j=0; j<lenM1; j++) {
				for(int k=0; k<lenCond; k++) {
					newC2.set(j, k, invC.get(j+lenCond,k) );
				}
			}		
			
			//extract first rows from mean to m2
			for(int j=0; j<lenCond; j++) {
				m2.set(j,0,m.get(j,0));
			}
			//extract last rows from mean to m1
			for(int j=0; j<lenM1; j++) {
				m1.set(j,0,m.get(j+lenCond,0));
			}
			SimpleMatrix invNewC1 = newC1.invert();
			// calculate new mean and new covariance of conditional distribution
			SimpleMatrix condMean = m1.minus( invNewC1.mult(newC2).mult( condition.minus(m2) ) );
			SimpleMatrix condCovariance = invNewC1;
			conditionalMeans.add(condMean);
			conditionalCovs.add(condCovariance);
			
			// calculate new weights
			
			// extract all elements from inverse covariance that correspond only to m2
			// that means extract the block in the left top corner with height=width=lenCond
			SimpleMatrix newC22 = new SimpleMatrix(lenCond,lenCond);
			for(int j=0; j<lenCond; j++) {
				for(int k=0; k<lenCond; k++) {
					newC22.set(j, k, c.get(j,k) );
				}
			}
			double mahalanobisDistance = condition.minus(m2).transpose().mult(newC22.invert()).mult(condition.minus(m2)).trace();
			double newWeight = ((1 / (a * Math.sqrt(newC22.determinant()))) * Math.exp((-0.5d) * mahalanobisDistance))* weights.get(i);
			conditionalWeights.add(newWeight);
		}
		// normalize weights
		double weightSum = 0;
		for(int i=0; i<conditionalWeights.size(); i++) {
			weightSum += conditionalWeights.get(i);
		}
		for(int i=0; i<conditionalWeights.size(); i++) {
			double weight = conditionalWeights.get(i);
			weight = weight /weightSum;
			conditionalWeights.set(i,weight);
		}
		result = new ConditionalDistribution(conditionalMeans, conditionalCovs, conditionalWeights);
		return result;
	}
	

	/**
	 * Find Maximum by gradient-quadratic search.
	 * First a conditional distribution is derived from the kde.
	 * @param start
	 * @return
	 */
	public SearchResult gradQuadrSearch(SimpleMatrix start){
		
		
		SimpleMatrix condVector = new SimpleMatrix(4,1);
		for(int i=0; i<condVector.numRows(); i++){
			condVector.set(i,0,start.get(i,0));
		}
		ConditionalDistribution conditionalDist = getConditionalDistribution(condVector);
		
		ArrayList<SimpleMatrix> means = conditionalDist.conditionalMeans;
		ArrayList<SimpleMatrix> covs = conditionalDist.conditionalCovs;
		ArrayList<Double> weights = conditionalDist.conditionalWeights;

		SimpleMatrix gradient = new SimpleMatrix(2,1);
		SimpleMatrix hessian = new SimpleMatrix(2,2);
		double n = means.get(0).numRows();
		double a = Math.pow(Math.sqrt(2 * Math.PI), n);
		
		SimpleMatrix x = new SimpleMatrix(2,1);
		x.set(0,0,start.get(start.numRows()-2,0));
		x.set(1,0,start.get(start.numRows()-1,0));
		ArrayList<Double> mahalanobisDistances;
		double step = 1;
		double probability = 0;
		SimpleMatrix gradStep = null;
		do {
			mahalanobisDistances = mahalanobis(x, means, covs);
			//calculate gradient and hessian:
			double prob = 0;
			for (int i = 0; i < means.size(); i++) {
				// check wether the component actually contributes to to the density at given point 
				if(mahalanobisDistances.get(i) < MAX_MAHALANOBIS_DIST) {
					SimpleMatrix m = means.get(i);

					SimpleMatrix dm = m.minus(x);
					SimpleMatrix c = covs.get(i);

					
					SimpleMatrix invC = c.invert();
					double w = weights.get(i);
					//probability p(x,m)
					double p = ((1 / (a * Math.sqrt(c.determinant()))) * Math.exp((-0.5d) * mahalanobisDistances.get(i))) * w;
					prob += p; 
					gradient = gradient.plus( invC.mult(dm).scale(p) );
					hessian = hessian.plus( invC.mult( dm.mult(dm.transpose()).minus(c) ).mult(invC).scale(p) );
				}


			}
			// save x
			SimpleMatrix xOld = new SimpleMatrix(x);
			SimpleEVD hessianEVD = hessian.eig();
			int maxEVIndex = hessianEVD.getIndexMax();
			if(hessianEVD.getEigenvalue(maxEVIndex).getReal() < 0){
				gradStep = hessian.invert().mult(gradient);
				x = xOld.minus(gradStep);
			}
			double prob1 = 	evaluate(x, means, covs, weights);
			if( prob1 <= prob | hessianEVD.getEigenvalue(maxEVIndex).getReal() >= 0) {
				gradStep = gradient.scale(step);
				x = xOld.plus(gradStep);
				while(evaluate(x, means, covs, weights) < prob){
					step = step/2;
					gradStep = gradient.scale(step);
					x = xOld.plus(gradStep);
				}
			}
			probability =	evaluate(x, means, covs, weights); 
		}while(gradStep.elementMaxAbs() > 1E-10);
		
		return new SearchResult(x, probability);
	}

	@Override
	public double evaluate(SimpleMatrix pointVector) {
		ArrayList<SimpleMatrix> means = new ArrayList<SimpleMatrix>();
		ArrayList<SimpleMatrix> covs = new ArrayList<SimpleMatrix>();
		ArrayList<Double> weights = new ArrayList<Double>();
		means = this.getSubMeans();
		covs = this.getSubSmoothedCovariances();
		weights = this.getSubWeights();
		double d = 0d;
		double n = means.get(0).numRows();
		double a = Math.pow(Math.sqrt(2 * Math.PI), n);
		
		ArrayList<Double> mahalanobisDistances = mahalanobis(pointVector, means, covs);

		for (int i = 0; i < means.size(); i++) {
			// check wether the component actually contributes to to the density at given point 
			if(mahalanobisDistances.get(i) < MAX_MAHALANOBIS_DIST) {
				SimpleMatrix m = means.get(i);
				SimpleMatrix c = covs.get(i);
				double w = weights.get(i);
				d += ((1 / (a * Math.sqrt(c.determinant()))) * Math.exp((-0.5d) * mahalanobisDistances.get(i))) * w;
			}
		}
		return d;
	}
	
	public double evaluate(SimpleMatrix pointVector, ArrayList<SimpleMatrix> means, ArrayList<SimpleMatrix> covs, ArrayList<Double> weights) {
		double d = 0d;
		double n = means.get(0).numRows();
		double a = Math.pow(Math.sqrt(2 * Math.PI), n);
		
		ArrayList<Double> mahalanobisDistances = mahalanobis(pointVector, means, covs);

		for (int i = 0; i < means.size(); i++) {
			// check wether the component actually contributes to to the density at given point 
			if(mahalanobisDistances.get(i) < MAX_MAHALANOBIS_DIST) {
				SimpleMatrix m = means.get(i);
				SimpleMatrix c = covs.get(i);
				double w = weights.get(i);
				d += ((1 / (a * Math.sqrt(c.determinant()))) * Math.exp((-0.5d) * mahalanobisDistances.get(i))) * w;
			}
		}
		return d;
	}
	
	public ArrayList<Double> mahalanobis(SimpleMatrix x, ArrayList<SimpleMatrix> means, ArrayList<SimpleMatrix> covs) {
		ArrayList<Double> mahalanobisDistances = new java.util.ArrayList<Double>();
		for (int i = 0; i < means.size(); i++) {
			SimpleMatrix m = means.get(i);
			SimpleMatrix c = covs.get(i);
			double distance = x.minus(m).transpose().mult(c.invert()).mult(x.minus(m)).trace();
			mahalanobisDistances.add(distance);
		}
		return mahalanobisDistances;
	}

	
	public void resetProbabilityCache(){
		mProbabilityCache.clear();
	}
				
	/**
	 * Evaluates the distribution at the given n-dimensional points and returns
	 * the results in a List of double-values.
	 * 
	 * @param points
	 * @return array of double values
	 */
	@Override
	public ArrayList<Double> evaluate(ArrayList<SimpleMatrix> points) {
		ArrayList<Double> resultPoints = new ArrayList<Double>();
		for (SimpleMatrix point : points) {
			resultPoints.add(evaluate(point));
		}
		return resultPoints;
	}

	@Override
	public void setBandwidthMatrix(SimpleMatrix mBandwidthMatrix) {
		this.mBandwidthMatrix = mBandwidthMatrix;
		for (BaseSampleDistribution d : mSubDistributions) {
			d.setBandwidthMatrix(mBandwidthMatrix);
		}
	}

	public void setNoOfCompsThreshold(float threshold) {
		mNoOfCompsThreshold = threshold;
	}

	public float getNoOfCompsThreshold() {
		return this.mNoOfCompsThreshold;
	}
}
