package de.tuhh.luethke.oKDE.model;

import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;
import java.util.List;

import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;
import org.ejml.simple.SimpleMatrix;
import org.ejml.simple.SimpleSVD;

import de.tuhh.luethke.oKDE.Exceptions.EmptyDistributionException;
import de.tuhh.luethke.oKDE.utility.Compressor;
import de.tuhh.luethke.oKDE.utility.MomentMatcher;

public class SampleModel extends BaseSampleDistribution {

	private static final float DEFAULT_NO_OF_COMPS_THRES = 6;
	
	// When mahalanobis distance gets bigger than this the component does not contribute
	// to the density that should be calculated
	// exp(-50) ~ 10E-22
	private static final double MAX_MAHALANOBIS_DIST = 50;

	// compression threshold (maximal hellinger distance)
	public double mCompressionThreshold;
	
	// effective number of observed samples
	protected double mEffectiveNoOfSamples;

	// component distributions
	protected ArrayList<BaseSampleDistribution> mSubDistributions;

	// threshold to determine when compression is necessary
	public float mNoOfCompsThreshold;
	
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
	
	private ArrayList<Double> mahalanobis(SimpleMatrix x, ArrayList<SimpleMatrix> means, ArrayList<SimpleMatrix> covs) {
		ArrayList<Double> mahalanobisDistances = new java.util.ArrayList<Double>();
		for (int i = 0; i < means.size(); i++) {
			SimpleMatrix m = means.get(i);
			SimpleMatrix c = covs.get(i);
			double distance = x.minus(m).transpose().mult(c.invert()).mult(x.minus(m)).trace();
			mahalanobisDistances.add(distance);
		}
		return mahalanobisDistances;
	}
	
	/**
	 * Evaluates the conditional probability using the given point vector. 
	 * 
	 * Example: 
	 * If the point vector is (4,1,3,2) and the conditional dimensions are 0,1 the probability of (4,1,3,2) under the condition 
	 * that x1=4, x2=1 is evaluated.
	 * 
	 * @param pointVector This vector defines where the pdf should be evaluated
	 * @param condDim An array containing the dimensions that shall be taken as conditional variables
	 * @return The conditional probability 
	 */
	public double evaluateConditional(SimpleMatrix pointVector, int[] condDim) {
		long time = System.currentTimeMillis();
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
		double marg = marginal(pointVector,condDim);
		return d/marg;
	}
	
	
	public double marginal(SimpleMatrix pointVector, int[] margDimensions){
		SimpleMatrix tmpMatrix = new SimpleMatrix(margDimensions.length,1);
		for(int i=0; i<margDimensions.length; i++) {
			tmpMatrix.set(i,0,pointVector.get(margDimensions[i],0));
		}
		pointVector = tmpMatrix;
		ArrayList<SimpleMatrix> means = new ArrayList<SimpleMatrix>();
		ArrayList<SimpleMatrix> covs = new ArrayList<SimpleMatrix>();
		ArrayList<Double> weights = new ArrayList<Double>();
		means = this.getSubMeans();
		covs = this.getSubSmoothedCovariances();
		weights = this.getSubWeights();
		for(int i=0; i<means.size(); i++){
			double[][] m = new double[margDimensions.length][1];
			for(int j=0; j<margDimensions.length; j++) {
				m[j][0] = means.get(i).get(margDimensions[j], 0);
			}
			means.set(i, new SimpleMatrix(m));
		}
		for(int i=0; i<covs.size(); i++){
			double[][] m = new double[margDimensions.length][margDimensions.length];
			for(int j=0; j<margDimensions.length; j++) {
				for(int k=0; k<margDimensions.length; k++) {
					m[j][k] = covs.get(i).get(margDimensions[j], margDimensions[k]);
				}
			}
			covs.set(i, new SimpleMatrix(m));
		}
		
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
	
	/**
	 * Finds the volume under the surface described by this KDE. The point that is given as the center
	 * parameter is taken and a square is defined around this point 
	 * Using xSegs number of segments across the x axis and ySegs number of segments across the y axis. 
	 * @param center A column vector 
	 * @param xSegs The number of segments in the x axis.
	 * @param ySegs The number of segments in the y axis.
	 * @return The volume under the pdf.
	 */
	public double cummulativeConditional(SimpleMatrix fixed , SimpleMatrix center, double squareWidth, int xSegs, int ySegs) {
    	int[] condDim = new int[center.numRows()];
    	for(int i=0; i<condDim.length; i++){
    		condDim[i] = i;
    	}
		double a = center.get(0,0)-squareWidth/2;
		double b = center.get(0,0)+squareWidth/2;
		double c = center.get(1,0)-squareWidth/2;
		double d = center.get(1,0)+squareWidth/2;
	    double xSegSize = (b - a) / xSegs; // length of an x segment.
	    double ySegSize = (d - c) / ySegs; // length of a y segment.
	    double volume = 0; // volume under the surface.
    	double[][] tmp = new double[fixed.numRows()+2][1];
    	for(int i=0; i<fixed.numRows(); i++) {
    		tmp[i][0] = fixed.get(i,0);
    	}
    	int dim = tmp.length-1;
	    for (int i = 0; i < xSegs; i++) {
	        for (int j = 0; j < ySegs; j++) {
	        	tmp[dim-1][0] = a + (xSegSize * i);
	        	tmp[dim][0] = c + (ySegSize * j);
	            double height = evaluateConditional(new SimpleMatrix(tmp),condDim);
	            tmp[dim-1][0] = a + (xSegSize * (i + 1));
	        	tmp[dim][0] = c + (ySegSize * j);
	            height += evaluateConditional(new SimpleMatrix(tmp),condDim);
	            tmp[dim-1][0] = a + (xSegSize * (i + 1));
	        	tmp[dim][0] = c + (ySegSize * (j + 1));
	            height += evaluateConditional(new SimpleMatrix(tmp),condDim);
	            tmp[dim-1][0] = a + (xSegSize * i);
	        	tmp[dim][0] = c + (ySegSize * (j + 1));
	            height += evaluateConditional(new SimpleMatrix(tmp),condDim);
	            height /= 4;

	            // height is the average value of the corners of the current segment.
	            // We can use the average value since a box of this height has the same volume as the original segment shape.

	            // Add the volume of the box to the volume.
	            if (height>1)
	            	System.out.println(1);
	            volume += xSegSize * ySegSize * height;
	        }
	    }

	    return volume;
	}
	
	public double cummulativeMarginal(SimpleMatrix fixed , SimpleMatrix center, double squareWidth, int xSegs, int ySegs) {
		int[] condDim = new int[center.numRows()];
    	for(int i=0; i<condDim.length; i++){
    		condDim[i] = i;
    	}
		double a = center.get(0,0)-squareWidth/2;
		double b = center.get(0,0)+squareWidth/2;
		double c = center.get(1,0)-squareWidth/2;
		double d = center.get(1,0)+squareWidth/2;
	    double xSegSize = (b - a) / xSegs; // length of an x segment.
	    double ySegSize = (d - c) / ySegs; // length of a y segment.
	    double volume = 0; // volume under the surface.
    	double[][] tmp = new double[fixed.numRows()+2][1];
    	for(int i=0; i<fixed.numRows(); i++) {
    		tmp[i][0] = fixed.get(i,0);
    	}
    	int dim = tmp.length-1;
	    for (int i = 0; i < xSegs; i++) {
	        for (int j = 0; j < ySegs; j++) {
	        	tmp[dim-1][0] = a + (xSegSize * i);
	        	tmp[dim][0] = c + (ySegSize * j);
	            double height = marginal(new SimpleMatrix(tmp),condDim);
	            tmp[dim-1][0] = a + (xSegSize * (i + 1));
	        	tmp[dim][0] = c + (ySegSize * j);
	            height += marginal(new SimpleMatrix(tmp),condDim);
	            tmp[dim-1][0] = a + (xSegSize * (i + 1));
	        	tmp[dim][0] = c + (ySegSize * (j + 1));
	            height += marginal(new SimpleMatrix(tmp),condDim);
	            tmp[dim-1][0] = a + (xSegSize * i);
	        	tmp[dim][0] = c + (ySegSize * (j + 1));
	            height += marginal(new SimpleMatrix(tmp),condDim);
	            height /= 4;

	            // height is the average value of the corners of the current segment.
	            // We can use the average value since a box of this height has the same volume as the original segment shape.

	            // Add the volume of the box to the volume.
	            if (height>1)
	            	System.out.println(1);
	            volume += xSegSize * ySegSize * height;
	        }
	    }

	    return volume;
	}
	
	public double trapezoidRule1(SimpleMatrix fixed , SimpleMatrix bounds, int xSegs, int ySegs) {
		double a = bounds.get(0,0)-5;
		double b = bounds.get(0,0)+5;
		double c = bounds.get(1,0)-5;
		double d = bounds.get(1,0)+5;
	    double xSegSize = (b - a) / xSegs; // length of an x segment.
	    double ySegSize = (d - c) / ySegs; // length of a y segment.
	    double volume = 0; // volume under the surface.
    	double[][] tmp = new double[fixed.numRows()+2][1];
    	for(int i=0; i<fixed.numRows(); i++) {
    		tmp[i][0] = fixed.get(i,0);
    	}
    	int dim = tmp.length-1;
    	int[] condDim = {0,1};
	    for (int i = 0; i < xSegs; i++) {
	        for (int j = 0; j < ySegs; j++) {
	        	tmp[dim-1][0] = a + (xSegSize * i);
	        	tmp[dim][0] = c + (ySegSize * j);
	            double height = evaluateConditional(new SimpleMatrix(tmp),condDim);
	            tmp[dim-1][0] = a + (xSegSize * (i + 1));
	        	tmp[dim][0] = c + (ySegSize * j);
	            height += evaluateConditional(new SimpleMatrix(tmp),condDim);
	            tmp[dim-1][0] = a + (xSegSize * (i + 1));
	        	tmp[dim][0] = c + (ySegSize * (j + 1));
	            height += evaluateConditional(new SimpleMatrix(tmp),condDim);
	            tmp[dim-1][0] = a + (xSegSize * i);
	        	tmp[dim][0] = c + (ySegSize * (j + 1));
	            height += evaluateConditional(new SimpleMatrix(tmp),condDim);
	            height /= 4;

	            // height is the average value of the corners of the current segment.
	            // We can use the average value since a box of this height has the same volume as the original segment shape.

	            // Add the volume of the box to the volume.
	            volume += xSegSize * ySegSize * height;
	        }
	    }

	    return volume;
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
