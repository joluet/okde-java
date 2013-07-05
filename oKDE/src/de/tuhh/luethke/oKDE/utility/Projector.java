package de.tuhh.luethke.oKDE.utility;

import java.util.List;

import org.ejml.simple.SimpleMatrix;
import org.ejml.simple.SimpleSVD;

import de.tuhh.luethke.oKDE.Exceptions.EmptyDistributionException;
import de.tuhh.luethke.oKDE.model.SampleDist;

public class Projector {
    private static final double MIN_VALUE = 1E-7d; 
    
    private static final double CONST_SMALL_FACTOR = 1E-10d; 
    
    public static void projectSampleDistToSubspace(SampleDist distribution) throws EmptyDistributionException {
	/*MomentMatcher.matchMoments(distribution);*/
	SimpleMatrix globalCov = distribution.getGlobalCovariance();
	int globalCovSize = globalCov.numRows();
	SimpleMatrix identity = SimpleMatrix.identity(globalCovSize);
	globalCov = globalCov.plus(identity.scale(CONST_SMALL_FACTOR));
	int d = globalCov.numCols();
	// calculate eigen directions --> determine the subspace
	SimpleSVD svd = globalCov.svd(true);
	SimpleMatrix U = svd.getU();
	SimpleMatrix S = svd.getW();
	SimpleMatrix V = svd.getV();
	V = U;  
	
	SimpleMatrix s = S.extractDiag();
	s = s.scale(1/S.elementMaxAbs());
	
	double[] validElements = new double[d];
	int countValidElements = 0;
	boolean isCompletelySingular = false;
	SimpleMatrix invS = null;
	if(s.elementMaxAbs() < MIN_VALUE || Double.isNaN(s.elementSum())){
	    S = SimpleMatrix.identity(S.numCols()).transpose();
	    invS = SimpleMatrix.identity(d).scale(2d/MIN_VALUE);
	    for(int i=0; i<validElements.length; i++)
		validElements[i] = 1;	
	    countValidElements = validElements.length;
	    isCompletelySingular = true;
	}else{
	    
	    for(int i=0; i<validElements.length; i++){
		if(s.get(i, 0) > MIN_VALUE){
		    validElements[i] = 1;
		    countValidElements++;
		}else
		    validElements[i] = s.get(i,0);
	    }

	    
	    S = MatrixOps.elemPow(S, -1);
	    invS = SimpleMatrix.identity(d).scale(0);
	    
	    for(int i=0; i<validElements.length; i++){
		if(validElements[i] == 1)
		    invS.set(i,i,S.get(i,i));
		else
		    invS.set(i,i,1d/validElements[i]);
	    }
	}
	SimpleMatrix trnsF = MatrixOps.elemSqrt(MatrixOps.abs(invS)).mult(V.invert());	    
	
	// forward transform the pdf and remove non-valid eigendirections
	SimpleMatrix trnsBandwidthMatrix = transformMatrix(trnsF, distribution.getmBandwidthMatrix(),
		validElements, countValidElements);
	System.out.println(trnsBandwidthMatrix);
	
	List<SimpleMatrix> oldMeans = distribution.getMeans();
	
	SimpleMatrix globalMean = distribution.getGlobalMean();
	
	for(int i=0; i<oldMeans.size(); i++){
	    oldMeans.set(i, oldMeans.get(i).minus(globalMean));
	}
	
	List<SimpleMatrix> oldCovs = distribution.getCovariances();
	for(int i=0; i<oldCovs.size(); i++){
	    oldCovs.set(i, transformMatrix(trnsF, oldCovs.get(i),
			validElements, countValidElements));
	}
	
    }

    private static SimpleMatrix transformMatrix(SimpleMatrix trnsF, SimpleMatrix matrix, double[] validElements, 
	    int countValidElements) {
	// forward transform the pdf and remove non-valid eigendirections
	SimpleMatrix tmp = trnsF.mult(matrix).mult(
		trnsF.transpose());
	SimpleMatrix trnsMatrix = new SimpleMatrix(countValidElements,
		countValidElements);
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
