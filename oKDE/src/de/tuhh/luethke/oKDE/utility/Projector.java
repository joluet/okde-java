package de.tuhh.luethke.oKDE.utility;

import java.util.Arrays;
import java.util.List;

import org.ejml.simple.SimpleMatrix;
import org.ejml.simple.SimpleSVD;

import de.tuhh.luethke.oKDE.Exceptions.EmptyDistributionException;
import de.tuhh.luethke.oKDE.model.SampleDist;

public class Projector {
	private static final double MIN_VALUE = 1E-7d;

	private static final double CONST_SMALL_FACTOR = 1E-10d;

	public static void projectSampleDistToSubspace(SampleDist distribution) throws EmptyDistributionException {
		MomentMatcher.matchMoments(distribution);
		SimpleMatrix globalCov = distribution.getGlobalCovariance();
		SimpleMatrix globalMean = distribution.getGlobalMean();
		int globalCovSize = globalCov.numRows();
		SimpleMatrix identity = SimpleMatrix.identity(globalCovSize);
		globalCov = globalCov.plus(identity.scale(CONST_SMALL_FACTOR));
		int d = globalCov.numCols();
		// calculate eigen directions --> determine the subspace
		SimpleSVD svd = globalCov.svd(true);
		double[][] data = { { -0.696148448325867, -0.717897860348871 }, { -0.717897860348871, 0.696148448325867 } };
		SimpleMatrix U = new SimpleMatrix(data);// svd.getU();
		double[][] data1 = { { 10.092428727058593, 0 }, { 0, 0.018741000612343 } };
		SimpleMatrix S = new SimpleMatrix(data1);// svd.getW();
		System.out.println("S" + S);
		SimpleMatrix V = new SimpleMatrix(data);// svd.getV();

		V = U;

		SimpleMatrix s = S.extractDiag();
		s = s.scale(1 / S.elementMaxAbs());

		Double[] validElements = new Double[d];
		int countValidElements = 0;
		boolean isCompletelySingular = false;
		SimpleMatrix invS = null;
		if (s.elementMaxAbs() < MIN_VALUE || Double.isNaN(s.elementSum())) {
			S = SimpleMatrix.identity(S.numCols()).transpose();
			invS = SimpleMatrix.identity(d).scale(2d / MIN_VALUE);
			for (int i = 0; i < validElements.length; i++)
				validElements[i] = 1d;
			countValidElements = validElements.length;
			isCompletelySingular = true;
		} else {

			for (int i = 0; i < validElements.length; i++) {
				if (s.get(i, 0) > MIN_VALUE) {
					validElements[i] = 1d;
					countValidElements++;
				} else
					validElements[i] = s.get(i, 0);
			}

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
			List<SimpleMatrix> subSubMeans = subDistributions.get(i).getSubMeans();
			List<SimpleMatrix> subSubCovs = subDistributions.get(i).getSubCovariances();
			if(subSubMeans.size() > 0)
				for (int j = 0; j < subSubMeans.size(); j++) {
					subSubCovs.set(j, transformMatrix(trnsF, subSubCovs.get(j), validElements, countValidElements));
					SimpleMatrix tmp = trnsF.mult(subSubMeans.get(i).minus(globalMean));
					tmp = MatrixOps.deleteElementsFromVector(tmp, Arrays.asList(validElements));
					subSubMeans.set(i, tmp);
					
				}
			else {
				SimpleMatrix subMean = subDistributions.get(i).getGlobalMean();
				SimpleMatrix subCov = subDistributions.get(i).getGlobalCovariance();
				SimpleMatrix tmp = trnsF.mult(subMean.minus(globalMean));
				tmp = MatrixOps.deleteElementsFromVector(tmp, Arrays.asList(validElements));
				subDistributions.get(i).setGlobalMean(tmp);
				subDistributions.get(i).setGlobalCovariance(transformMatrix(trnsF, subCov, validElements, countValidElements));
			}
			MomentMatcher.matchMoments(subDistributions.get(i));
			SimpleMatrix subCov = subDistributions.get(i).getGlobalCovariance();
			subDistributions.get(i).setGlobalCovariance(subCov.plus(trnsBandwidthMatrix));
		}
		/*% transform also the global covariance
		globalCov = F_trns*globalCov*F_trns' ;
		globalCov = globalCov(id_valid, id_valid) ;*/
		globalCov = trnsF.mult(globalCov).mult(globalCov);
		globalCov = transformMatrix(trnsF, globalCov, validElements, countValidElements);
		globalCov = globalCov;

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
