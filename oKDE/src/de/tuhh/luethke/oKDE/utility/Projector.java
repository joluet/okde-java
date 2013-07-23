package de.tuhh.luethke.oKDE.utility;

import java.util.Arrays;
import java.util.List;

import org.ejml.simple.SimpleMatrix;
import org.ejml.simple.SimpleSVD;

import de.tuhh.luethke.oKDE.Exceptions.EmptyDistributionException;
import de.tuhh.luethke.oKDE.Exceptions.NoOfArgumentsException;
import de.tuhh.luethke.oKDE.model.OneComponentDistribution;
import de.tuhh.luethke.oKDE.model.BaseSampleDistribution;
import de.tuhh.luethke.oKDE.model.SampleModel;
import de.tuhh.luethke.oKDE.model.TwoComponentDistribution;

/**
 * This class transforms a sample distribution that is lying in a degenerate
 * space into the non-degenerate subspace by projecting
 * 
 * @author Jonas Luethke
 * 
 */
public class Projector {
	private static final double MIN_VALUE = 1E-7d;

	private static final double CONST_SMALL_FACTOR = 1E-10d;

	/**
	 * Forward transform a sample distribution
	 * 
	 * @param distribution
	 * @throws EmptyDistributionException
	 */
	public static ProjectionData projectSampleDistToSubspace(SampleModel distribution) throws EmptyDistributionException {
		ProjectionData projectionData = new ProjectionData();
		MomentMatcher.matchMoments(distribution);
		SimpleMatrix globalSmoothedCov = distribution.getmGlobalCovarianceSmoothed();
		SimpleMatrix globalMean = distribution.getGlobalMean();
		int globalCovSize = globalSmoothedCov.numRows();
		SimpleMatrix identity = SimpleMatrix.identity(globalCovSize);
		globalSmoothedCov = globalSmoothedCov.plus(identity.scale(CONST_SMALL_FACTOR));
		int d = globalSmoothedCov.numCols();
		// calculate eigen directions --> determine the subspace
		projectionData.mSVD = globalSmoothedCov.svd(true);
		SimpleSVD<?> svd = globalSmoothedCov.svd(true);
		SimpleMatrix U = svd.getU();
		SimpleMatrix S = svd.getW();
		SimpleMatrix V = svd.getV();
		
		
		V = U;

		SimpleMatrix s = S.extractDiag();
		s = s.scale(1 / S.elementMaxAbs());

		Double[] validElements = new Double[d];
		int countValidElements = 0;
		SimpleMatrix invS = null;
		// if S is almost zero --> singular
		if (s.elementMaxAbs() < MIN_VALUE || Double.isNaN(s.elementSum())) {
			S = SimpleMatrix.identity(S.numCols()).transpose();
			invS = SimpleMatrix.identity(d).scale(2d / MIN_VALUE);
			for (int i = 0; i < validElements.length; i++)
				validElements[i] = 1d;
			countValidElements = validElements.length;
		} else {
			// if S ins not completely singular: check what elements are >0
			// store them in validElements
			for (int i = 0; i < validElements.length; i++) {
				if (s.get(i, 0) > MIN_VALUE) {
					validElements[i] = 1d;
					countValidElements++;
				} else
					validElements[i] = s.get(i, 0);
			}
			// create inverse matrix of S
			S = MatrixOps.elemPow(S, -1);
			invS = SimpleMatrix.identity(d).scale(0);
			for (int i = 0; i < validElements.length; i++) {
				if (validElements[i] == 1)
					invS.set(i, i, S.get(i, i));
				else
					invS.set(i, i, 1d / validElements[i]);
			}
		}
		//TODO: remove
		validElements = new Double[2];
		validElements[0] = 1d;
		validElements[1] = 0d;
		countValidElements = 1;
		
		SimpleMatrix trnsF = MatrixOps.elemSqrt(MatrixOps.abs(invS)).mult(V.invert());

		// forward transform the pdf and remove non-valid eigendirections
		SimpleMatrix trnsBandwidthMatrix = transformMatrix(trnsF, distribution.getBandwidthMatrix(), validElements, countValidElements);
		List<BaseSampleDistribution> subDistributions = distribution.getSubDistributions();
		for (int i = 0; i < subDistributions.size(); i++) {
			SimpleMatrix[] subSubMeans = new SimpleMatrix[1];
			SimpleMatrix[] subSubCovs = new SimpleMatrix[1];
			if (subDistributions.get(i).getClass() == TwoComponentDistribution.class) {
				TwoComponentDistribution subDist = (TwoComponentDistribution) subDistributions.get(i);
				subSubMeans = subDist.getSubMeans();
				subSubCovs = subDist.getSubCovariances();
			} else {
				OneComponentDistribution subDist = (OneComponentDistribution) subDistributions.get(i);
				subSubMeans[0] = subDist.getGlobalMean();
				subSubCovs[0] = subDist.getmGlobalCovarianceSmoothed();
			}
			if (subSubMeans.length > 1) {
				for (int j = 0; j < subSubMeans.length; j++) {
					subSubCovs[j] = transformMatrix(trnsF, subSubCovs[j], validElements, countValidElements);
					SimpleMatrix tmp = trnsF.mult(subSubMeans[i].minus(globalMean));
					tmp = MatrixOps.deleteElementsFromVector(tmp, Arrays.asList(validElements), countValidElements);
					subSubMeans[i] = tmp;
				}
				try {
					((TwoComponentDistribution) subDistributions.get(i)).setSubMeans(subSubMeans);
					((TwoComponentDistribution) subDistributions.get(i)).setSubCovariances(subSubCovs);
				} catch (NoOfArgumentsException e) {
					// This exception can't be thrown because sample models can
					// only contain subcomponents with at most 2 components
					e.printStackTrace();
				}
				MomentMatcher.matchMoments(((TwoComponentDistribution) subDistributions.get(i)));
			} else {
				SimpleMatrix subMean = subDistributions.get(i).getGlobalMean();
				SimpleMatrix tmp = trnsF.mult(subMean.minus(globalMean));
				tmp = MatrixOps.deleteElementsFromVector(tmp, Arrays.asList(validElements), countValidElements);
				subDistributions.get(i).setGlobalMean(tmp);
			}
			subDistributions.get(i).setBandwidthMatrix(trnsBandwidthMatrix);
		}
		// transform also the global covariance
		globalSmoothedCov = transformMatrix(trnsF, globalSmoothedCov, validElements, countValidElements);
		distribution.setBandwidthMatrix(trnsBandwidthMatrix);
		projectionData.mCountValidElements = countValidElements;
		projectionData.mValidElements = validElements;
		return projectionData;
	}

	private static SimpleMatrix transformMatrix(SimpleMatrix trnsF, SimpleMatrix matrix, Double[] validElements, int countValidElements) {
		// forward transform the pdf and remove non-valid eigendirections
		SimpleMatrix tmp = trnsF.mult(matrix).mult(trnsF.transpose());
		SimpleMatrix trnsMatrix = new SimpleMatrix(countValidElements, countValidElements);
		for (int row = 0; row < trnsMatrix.numRows(); row++) {
			for (int column = 0; column < trnsMatrix.numCols(); column++) {
				if (validElements[row] == 1)
					trnsMatrix.set(row, column, tmp.get(row, column));
				else if (validElements[column] == 1)
					trnsMatrix.set(row, column, tmp.get(row, column));
				else
					trnsMatrix.set(0);
			}
		}
		return trnsMatrix;
	}
	
	private static SimpleMatrix backTransformMatrix(SimpleMatrix matrix, SimpleMatrix bandwidth, Double[] validElements) {
		// add removed eigendirections and backwards transform
		for (int row = 0; row < matrix.numRows(); row++) {
			for (int column = 0; column < matrix.numCols(); column++) {
				if (validElements[row] == 1)
					matrix.set(row, column, bandwidth.get(row, column));
				else if (validElements[column] == 1)
					matrix.set(row, column, bandwidth.get(row, column));
				else
					matrix.set(0);
			}
		}
		
		return matrix;
	}
	
	
	public static void projectSampleDistToOriginalSpace(SampleModel distribution, ProjectionData projectionData) throws EmptyDistributionException {
		// add missing nullspace to pdf
		/*num_nullDir = size(svdRes.S,1) - length(svdRes.id_valid) ;
		F_trns = svdRes.V*sqrt(svdRes.S) ;
		C_prot = zeros(size(svdRes.S)) ;*/
		SimpleMatrix bandwidth = distribution.getBandwidthMatrix();
		SimpleSVD<?> svd = projectionData.mSVD;
		SimpleMatrix U = svd.getU();
		SimpleMatrix S = svd.getW();
		SimpleMatrix V = svd.getV();
		Double[] validElements = projectionData.mValidElements;
		int countValidElements = projectionData.mCountValidElements;
		int noOfNullDirs = S.numCols() - countValidElements;
		SimpleMatrix trnsF = V.mult(MatrixOps.elemSqrt(S));
		SimpleMatrix protC = S.scale(0);
		protC = backTransformMatrix(protC, bandwidth, validElements);
		bandwidth = trnsF.mult(protC).mult(trnsF.transpose());
		
		List<BaseSampleDistribution> subDistributions = distribution.getSubDistributions();
		for (int i = 0; i < subDistributions.size(); i++) {
			SimpleMatrix[] subSubMeans = new SimpleMatrix[1];
			SimpleMatrix[] subSubCovs = new SimpleMatrix[1];
			if (subDistributions.get(i).getClass() == TwoComponentDistribution.class) {
				TwoComponentDistribution subDist = (TwoComponentDistribution) subDistributions.get(i);
				subSubMeans = subDist.getSubMeans();
				subSubCovs = subDist.getSubCovariances();
			} else {
				OneComponentDistribution subDist = (OneComponentDistribution) subDistributions.get(i);
				subSubMeans[0] = subDist.getGlobalMean();
				subSubCovs[0] = subDist.getmGlobalCovarianceSmoothed();
			}
			if (subSubMeans.length > 1) {
				for (int j = 0; j < subSubMeans.length; j++) {
					subSubCovs[j] = transformMatrix(trnsF, subSubCovs[j], validElements, countValidElements);
					SimpleMatrix tmp = trnsF.mult(subSubMeans[i].minus(globalMean));
					tmp = MatrixOps.deleteElementsFromVector(tmp, Arrays.asList(validElements));
					subSubMeans[i] = tmp;
				}
				try {
					((TwoComponentDistribution) subDistributions.get(i)).setSubMeans(subSubMeans);
					((TwoComponentDistribution) subDistributions.get(i)).setSubCovariances(subSubCovs);
				} catch (NoOfArgumentsException e) {
					// This exception can't be thrown because sample models can
					// only contain subcomponents with at most 2 components
					e.printStackTrace();
				}
				MomentMatcher.matchMoments(((TwoComponentDistribution) subDistributions.get(i)));
			} else {
				SimpleMatrix subMean = subDistributions.get(i).getGlobalMean();
				SimpleMatrix tmp = trnsF.mult(subMean.minus(globalMean));
				tmp = MatrixOps.deleteElementsFromVector(tmp, Arrays.asList(validElements));
				subDistributions.get(i).setGlobalMean(tmp);
			}
			subDistributions.get(i).setBandwidthMatrix(trnsBandwidthMatrix);
		}
	}
}
