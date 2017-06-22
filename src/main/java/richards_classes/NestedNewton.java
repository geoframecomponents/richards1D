/**
 *  @author Niccolo' Tubini
 */

package richards_classes;
/**
 * This class carries out the Nested-Newton algorithm
 * (A NESTED NEWTON-TYPE ALGORITHM FOR FINITE VOLUME METHODS SOLVING RICHARDS' EQUATION IN MIXED FORM, Casulli V., Zanolli P., Journal Scientific Computing, 2010)
 */

public class NestedNewton {
	private double outerResidual;
	private double innerResidual;
	
	int nestedNewton;
	int MAXITER_NEWT;
	int NUM_CONTROL_VOLUMES;
	
	double newtonTolerance;
	
	double[] psis;
	double[] mainDiagonal;
	double[] upperDiagonal;
	double[] lowerDiagonal;
	double[] rhss;
		
	double[] fs;
	double[] fks;
	double[] bb;
	double[] cc;
	double[] dis;
	double[] dpsis;
	double[] psis_outer;
	
	SoilParametrization soilPar;
	
	Thomas thomasAlg = new Thomas();
	JordanDecomposition jordanDecomposition;
	
	/**
	 * @param nestedNewton control parameter to choose between simple Newton method (0), or the nested Newton one (1)
	 * @param newtonTolerance prefixed tolerance representing the maximum mass balance error allowed  
	 * @param MAXITER_NEWT prefixed maximum number of iteration
	 * @param NUM_CONTROL_VOLUMES number of control volumes
	 * @param soilPar is the class to compute the soil hydraulic properties
	 */
	public NestedNewton(int nestedNewton, double newtonTolerance, int MAXITER_NEWT, int NUM_CONTROL_VOLUMES, SoilParametrization soilPar){
		this.nestedNewton = nestedNewton;
		this.newtonTolerance = newtonTolerance;
		this.MAXITER_NEWT = MAXITER_NEWT;
		this.NUM_CONTROL_VOLUMES = NUM_CONTROL_VOLUMES;
		this.soilPar = soilPar;
		
		jordanDecomposition = new JordanDecomposition(this.soilPar);
		fs			  = new double[this.NUM_CONTROL_VOLUMES];
		fks			  = new double[this.NUM_CONTROL_VOLUMES];
		bb            = new double[this.NUM_CONTROL_VOLUMES]; 
		cc 			  = new double[this.NUM_CONTROL_VOLUMES];
		dis			  = new double[this.NUM_CONTROL_VOLUMES];
		dpsis		  = new double[this.NUM_CONTROL_VOLUMES];
		psis_outer	  = new double[this.NUM_CONTROL_VOLUMES];
	}
	
	/**
	 * @param psis vector contains the suction values, it is a vector of length NUM_CONTROL_VOLUMES
	 * @param upperDiagonal upper diagonal of the coefficient matrix A of the linear system, it is a vector of length NUM_CONTROL_VOLUMES
	 * @param mainDiagonal main diagonal of the coefficient matrix A of the linear system, it is a vector of length NUM_CONTROL_VOLUMES
	 * @param lowerDiagonal lower diagonal of the coefficient matrix A of the linear system, it is a vector of length NUM_CONTROL_VOLUMES
	 * @param rhss right hand side term of the linear system, it is a vector of length NUM_CONTROL_VOLUMES
	 */
	public void set(double[] psis, double[] mainDiagonal, double[] upperDiagonal, double[] lowerDiagonal, double[] rhss){
		this.psis = psis;
		this.mainDiagonal = mainDiagonal;
		this.upperDiagonal = upperDiagonal;
		this.lowerDiagonal = lowerDiagonal;
		this.rhss = rhss;
	}
	
	public double[] solver(){

		// Initial guess of psis
		for(int i = 0; i < NUM_CONTROL_VOLUMES; i++) {
			psis[i] = Math.min(psis[i], soilPar.getPsiStar());
		}

		//// OUTER CYCLE ////
		for(int i = 0; i < MAXITER_NEWT; i++) {
			// I have to assign 0 to outerResidual otherwise I will take into account of the previous error
			outerResidual = 0.0;
			for(int j = 0; j < NUM_CONTROL_VOLUMES; j++) {
				fs[j] = soilPar.waterContent(psis[j]) - rhss[j];
				if(j == 0) {
					fs[j] = fs[j] + mainDiagonal[j]*psis[j] + upperDiagonal[j]*psis[j+1];
				}else if(j == NUM_CONTROL_VOLUMES -1) {
					fs[j] = fs[j] + lowerDiagonal[j]*psis[j-1] + mainDiagonal[j]*psis[j];
				}else {
					fs[j] = fs[j] + lowerDiagonal[j]*psis[j-1] + mainDiagonal[j]*psis[j] + upperDiagonal[j]*psis[j+1];
				}
				dis[j] = soilPar.dWaterContent(psis[j]);
				outerResidual += fs[j]*fs[j];
			}
			outerResidual = Math.pow(outerResidual,0.5);  
			System.out.println("   Outer iteration " + i + " with residual " +  outerResidual);
			if(outerResidual < newtonTolerance) {
				break;
			}
			if(nestedNewton == 0){
				bb = mainDiagonal.clone();
				cc = upperDiagonal.clone();
				for(int y = 0; y < NUM_CONTROL_VOLUMES; y++) {
					bb[y] += dis[y];
				}
				thomasAlg.set(cc,bb,lowerDiagonal,fs);
				dpsis = thomasAlg.solver();

				//// PSIS UPDATE ////
				for(int s = 0; s < NUM_CONTROL_VOLUMES; s++) {
					psis[s] = psis[s] - dpsis[s];
				}
			}else{
				psis_outer = psis.clone();

				// Initial guess of psis
				for(int j = 0; j < NUM_CONTROL_VOLUMES; j++) {
					psis[j] = Math.max(psis[j], soilPar.getPsiStar());
				}

				//// INNER CYCLE ////
				for(int j = 0; j < MAXITER_NEWT; j++) {
					// I have to assign 0 to innerResidual otherwise I will take into account of the previous error
					innerResidual = 0.0; 
					for(int l=0; l < NUM_CONTROL_VOLUMES; l++) {
						fks[l] = jordanDecomposition.waterContent1(psis[l]) - (jordanDecomposition.waterContent2(psis_outer[l]) + jordanDecomposition.dWaterContent2(psis_outer[l])*(psis[l] - psis_outer[l])) - this.rhss[l];
						if(l == 0) {
							fks[l] = fks[l] + mainDiagonal[l]*psis[l] + upperDiagonal[l]*psis[l+1];
						}else if(l == NUM_CONTROL_VOLUMES -1) {
							fks[l] = fks[l] + lowerDiagonal[l]*psis[l-1] + mainDiagonal[l]*psis[l];
						}else {
							fks[l] = fks[l] + lowerDiagonal[l]*psis[l-1] + mainDiagonal[l]*psis[l] + upperDiagonal[l]*psis[l+1];
						}
						dis[l] = jordanDecomposition.dWaterContent1(psis[l]) - jordanDecomposition.dWaterContent2(psis_outer[l]);
						innerResidual += fks[l]*fks[l];
					}
					innerResidual = Math.pow(innerResidual,0.5);

					System.out.println("     -Inner iteration " + j + " with residual " +  innerResidual);    

					if(innerResidual < newtonTolerance) {
						break;
					}

					//// THOMAS ALGORITHM////
					// Attention: the main diagonal of the coefficient matrix must not change!! The same for the upper diagonal

					bb = mainDiagonal.clone();
					cc = upperDiagonal.clone();
					for(int y = 0; y < NUM_CONTROL_VOLUMES; y++) {
						bb[y] += dis[y];
					}
					thomasAlg.set(cc,bb,lowerDiagonal,fks);
					dpsis = thomasAlg.solver();

					//// PSIS UPDATE ////
					for(int s = 0; s < NUM_CONTROL_VOLUMES; s++) {
						psis[s] = psis[s] - dpsis[s];
					}
				}
			} //// INNER CYCLE END ////
		}
		return psis;
	}
	
}
