package de.tuhh.luethke.oKDE.model;

import java.util.ArrayList;

import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;
import org.ejml.simple.SimpleMatrix;

public class SampleDist {
	
	private ArrayList<SimpleMatrix> mSampleData; 
	
	//kernel weights
	private double[] mMixWeights;
	
	private  ArrayList<Double> mWeights;

	
	private double mWeightSum;
	
	private double mForgettingFactor;
	
	private double mN_eff;
	
	//kernel covariances
	private SimpleMatrix[] mCovariance;
	
	public SampleDist(ArrayList<SimpleMatrix> sampleData, ArrayList<Double> weights,
		SimpleMatrix[] covariances) {
	    super();
	    mSampleData = sampleData;
	    mMixWeights = new double[2];
	    mCovariance = covariances;
	    mWeightSum = 0;
	    for(double w : weights){
		mWeightSum += w;
	    }
	    mWeights = new ArrayList<Double>();
	    mN_eff = weights.size();
	    mForgettingFactor = 1;
	}

	public double[] getWeights() {
	    return mMixWeights;
	}
	
	public void addWeights(double[] weights) {
	    double sumOfNewWeights = 0;
	    for(double w : weights) {
		sumOfNewWeights += w;
		mWeights.add(w);
	    }
	    
	    
	    mN_eff = mN_eff*mForgettingFactor + weights.length;
	    
	    mMixWeights[0] = mWeightSum/(mWeightSum*mForgettingFactor + sumOfNewWeights);
	    mMixWeights[1] = sumOfNewWeights/(mWeightSum*mForgettingFactor + sumOfNewWeights);
	    
	    mWeightSum = mWeightSum*mForgettingFactor + sumOfNewWeights;
	    
	    for(int i=0; i<mWeights.size()-weights.length; i++){
		mWeights.set(i, mWeights.get(i)*mMixWeights[0]);
	    }
	    for(int i=mWeights.size()-weights.length; i<mWeights.size(); i++){
		mWeights.set(i, mWeights.get(i)*mMixWeights[1]);
	    }
	    
	    System.out.println(mMixWeights[0]+"-"+mMixWeights[1]+" "+mWeightSum+" "+mWeights);
	}

	public ArrayList<SimpleMatrix> getSampleData() {
		return mSampleData;
	}
	public SimpleMatrix[] getCovariance() {
		return mCovariance;
	}

	public double getWeightSum() {
	    return mWeightSum;
	}
	
	public static void estimateBandwidth(SimpleMatrix mu, SimpleMatrix[] covariance, double[] weights, SimpleMatrix Cov_smp, int N_eff) {
	    

	  
	    double d = mu.numRows() ;
	    
	    
	    
	    SimpleMatrix G = Cov_smp.scale( Math.pow( (4/((d+2)*N_eff)),(2/(d+4)) ) ) ;

	    System.out.println("G: "+G);
	    System.out.println(G.get(0, 0));
	    float alphaScale = 1 ;
	    SimpleMatrix F = Cov_smp.scale(alphaScale);
	    
	    double Rf2 = getIntSquaredHessian( mu, weights, covariance, G ) ;
	    System.out.println("Rf2: "+Rf2);
	    double hAmise = Math.pow(( Math.pow(N_eff,(-1)) * Math.pow(F.determinant(),(-1/2))  /  ( Math.pow(Math.sqrt(4*Math.PI),d) * Rf2 * d )),(1/(d+4))) ;
	    
	    System.out.println("hAmise: "+hAmise);
	}
	
	public static double getIntSquaredHessian(SimpleMatrix mu, double[] weights, SimpleMatrix[] covariance, SimpleMatrix g) {
	    
	    long d = mu.numRows();
	    long N = mu.numCols();
	    System.out.println("d:" +d );
	    //normalizer
	    double constNorm = Math.pow((1/(2*Math.PI)),(d/2));
	    
	    double w1,w2,m,I=0,eta;
	    SimpleMatrix s1,s2,mu1, mu2,dm,ds;
	    for(int i1=0; i1<N; i1++){
		s1 = covariance[i1].plus(g) ;
		mu1 = mu.extractVector(false, i1);
		w1 = weights[i1];
		for(int i2=i1; i2<N; i2++){
		    s2 = covariance[i2];
		    mu2 = mu.extractVector(false, i2);
		    w2 = weights[i2];
		    SimpleMatrix A = s1.plus(s2).invert();
		    dm = mu1.minus(mu2);
		    /*if(true F is not identity){
			ds = dm.transpose().mult(A);
			b = ds'*ds ;
			B = A - 2*b ;
			C = A - b ;
			            
			            f_t = constNorm*sqrt(det(A))*exp(-0.5*ds*dm) ;
			            c = 2*trace(F*A*F*B) + trace(F*C)^2 ;
		    }else {*/
        		m = dm.transpose().mult(A).mult(dm).get(0);
        		System.out.println("m: "+m);
        		System.out.println("A:"+A);
        		double f_t = constNorm * Math.sqrt(A.determinant())
        			* Math.exp(-0.5 * m);
        		
        		DenseMatrix64F A_sqr = new DenseMatrix64F(A.numRows(),
        			A.numCols());
        		CommonOps.elementMult(A.getMatrix(), A.transpose().getMatrix(),
        			A_sqr);
        		double sum = CommonOps.elementSum(A_sqr);
        		double c = 2 * sum * (1 - 2 * m) + Math.pow((1 - m), 2)
        			* Math.pow(A.trace(), 2);
        		// }
        		// determine the weight of the current term
        		if (i1 == i2)
        		    eta = 1 ;
        		else
        		    eta = 2 ;
        		System.out.println("I: "+I+" f_t: "+f_t+" c: "+c+" w2*w1:"+(w2*w1));
        		I = I + f_t*c*w2*w1*eta ;
        		System.out.println("I: "+I);
		}
	    }
	    
	    
		        
	    
	    return I;
	}
	
	
	
	/*FAST VERSION? 
	private double getIntSquaredHessian(double[] muValues, double[] wValues, double[] cValues, double[] gValue) {
	    //double *muValues, *cValues, *wValues, *gValue, tmp, *ptr_out ;
	    
	    double A, dm, m, trFA_in, trFA_out, detA ;
	    int i, j, rowLen, colLen, k_i, idx_i, idx_j ;
	    double eta, I, const_norm ;
	    
	    if ( nrhs != 4 )
	        mexErrMsgTxt("To few input parameters! (Mu, w, Cov)") ;
	    	    
	    
	    //Get matrix mu
	    colLen = muValues.numCols();
	    rowLen = muValues.numRows();;
	     
	    if ( mxGetM(prhs[0]) != mxGetM(prhs[2]) )
	         mexErrMsgTxt("rows of Mu and Cov not equal!") ; 
	    if ( mxGetN(prhs[0]) != mxGetN(prhs[1]) )
	         mexErrMsgTxt("columns of Mu and Cov not equal!") ;
	    if ( mxGetN(prhs[0]) != mxGetN(prhs[2]) )
	         mexErrMsgTxt("columns of Mu and w not equal!") ;
	    
	    
	    

	    const_norm = Math.pow(1.0/(2.0*Math.PI), (((double)rowLen)/2.0)) ;
	    I = 0.0 ;
	    for( i = 0 ; i < colLen ; i++ ) {   
	        idx_i = i*rowLen ;
	        
	        for ( j = i ; j < colLen ; j++ ) {
	            idx_j = j*rowLen ;
	            if ( i == j ) {
	                eta = 1.0 ;
	            } else {
	                eta = 2.0 ;
	            }
	            
	            m = 0.0 ;
	            trFA_in = 0.0 ;
	            trFA_out = 0.0 ;
	            detA = 1.0 ;
	            for ( k_i = 0 ; k_i < rowLen ; ++k_i ) {
	                A = 1.0 / (gValue[0] + cValues[idx_i + k_i] + cValues[idx_j + k_i]) ;
	        
	                dm = muValues[idx_i + k_i] - muValues[idx_j + k_i] ;
	                m += dm*dm*A ;
	                
	                trFA_in += A*A ;
	                trFA_out += A ;
	                detA *= A ;
	            }
	            trFA_out *= trFA_out ;
	             

	            I += wValues[i] * wValues[j] * eta * ( const_norm*Math.sqrt(detA)*Math.exp(-0.5*m) ) * ( 2*trFA_in*(1.0-2.0*m) + (1.0-m)*(1.0-m)*trFA_out ) ;

	                
	        }
	    }
	}
	*/
}
