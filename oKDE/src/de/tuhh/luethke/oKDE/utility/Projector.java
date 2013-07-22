package de.tuhh.luethke.oKDE.utility;

import java.util.Arrays;
import java.util.List;

import org.ejml.simple.SimpleMatrix;
import org.ejml.simple.SimpleSVD;

import de.tuhh.luethke.oKDE.Exceptions.EmptyDistributionException;
import de.tuhh.luethke.oKDE.model.OneComponentDistribution;
import de.tuhh.luethke.oKDE.model.SampleDist;
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
	public static void projectSampleDistToSubspace(SampleModel distribution) throws EmptyDistributionException {
		MomentMatcher.matchMoments(distribution, true);
		SimpleMatrix globalCov = distribution.getGlobalCovariance();
		SimpleMatrix globalMean = distribution.getGlobalMean();
		int globalCovSize = globalCov.numRows();
		SimpleMatrix identity = SimpleMatrix.identity(globalCovSize);
		globalCov = globalCov.plus(identity.scale(CONST_SMALL_FACTOR));
		int d = globalCov.numCols();
		// calculate eigen directions --> determine the subspace
		SimpleSVD svd = globalCov.svd(true);
		SimpleMatrix U = svd.getU();
		SimpleMatrix S = svd.getW();
		System.out.println("S" + S);
		SimpleMatrix V = svd.getV();

		V = U;

		SimpleMatrix s = S.extractDiag();
		s = s.scale(1 / S.elementMaxAbs());

		Double[] validElements = new Double[d];
		int countValidElements = 0;
		boolean isCompletelySingular = false;
		SimpleMatrix invS = null;
		// if S is almost zero --> singular
		if (s.elementMaxAbs() < MIN_VALUE || Double.isNaN(s.elementSum())) {
			S = SimpleMatrix.identity(S.numCols()).transpose();
			invS = SimpleMatrix.identity(d).scale(2d / MIN_VALUE);
			for (int i = 0; i < validElements.length; i++)
				validElements[i] = 1d;
			countValidElements = validElements.length;
			isCompletelySingular = true;
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
		SimpleMatrix trnsF = MatrixOps.elemSqrt(MatrixOps.abs(invS)).mult(V.invert());

		// forward transform the pdf and remove non-valid eigendirections
		SimpleMatrix trnsBandwidthMatrix = transformMatrix(trnsF, distribution.getmBandwidthMatrix(), validElements, countValidElements);
		System.out.println(trnsBandwidthMatrix);

		List<SampleDist> subDistributions = distribution.getSubDistributions();
		for (int i = 0; i < subDistributions.size(); i++) {
			SimpleMatrix[] subSubMeans = new SimpleMatrix[1];
			SimpleMatrix[] subSubCovs = new SimpleMatrix[1];
			if (subDistributions.get(i).getClass() == TwoComponentDistribution.class) {
				TwoComponentDistribution subDist = (TwoComponentDistribution) subDistributions.get(i);
				subSubMeans = subDist.getSubMeans();
				subSubCovs = subDist.getSubSmoothedCovariances();
			} else {
				OneComponentDistribution subDist = (OneComponentDistribution) subDistributions.get(i);
				subSubMeans[1] = subDist.getGlobalMean();
				subSubCovs[1] = subDist.getmGlobalCovarianceSmoothed();
			}
			if (subSubMeans.length > 1) {
				for (int j = 0; j < subSubMeans.length; j++) {
					subSubCovs[j] = transformMatrix(trnsF, subSubCovs[j], validElements, countValidElements);
					SimpleMatrix tmp = trnsF.mult(subSubMeans[i].minus(globalMean));
					tmp = MatrixOps.deleteElementsFromVector(tmp, Arrays.asList(validElements));
					subSubMeans[i] = tmp;
				}
				((TwoComponentDistribution) subDistributions.get(i)).setSubMeans(subSubMeans[0], subSubMeans[1]);
				((TwoComponentDistribution) subDistributions.get(i)).setSubCovariances(subSubCovs[0], subSubCovs[1]);
				MomentMatcher.matchMoments(((TwoComponentDistribution) subDistributions.get(i)), true);
			} else {
				SimpleMatrix subMean = subDistributions.get(i).getGlobalMean();
				SimpleMatrix subCov = subDistributions.get(i).getmGlobalCovarianceSmoothed();
				SimpleMatrix tmp = trnsF.mult(subMean.minus(globalMean));
				tmp = MatrixOps.deleteElementsFromVector(tmp, Arrays.asList(validElements));
				subDistributions.get(i).setGlobalMean(tmp);
				subDistributions.get(i).setmGlobalCovarianceSmoothed(transformMatrix(trnsF, subCov, validElements, countValidElements));
			}
			SimpleMatrix subCov = subDistributions.get(i).getGlobalCovariance();
			subDistributions.get(i).setmGlobalCovarianceSmoothed(subCov.plus(trnsBandwidthMatrix));
		}
		// transform also the global covariance
		globalCov = transformMatrix(trnsF, globalCov, validElements, countValidElements);
		System.out.println(globalCov);

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
}
