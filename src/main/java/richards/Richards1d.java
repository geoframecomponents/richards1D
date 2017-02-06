package richards;

/*
 * GNU GPL v3 License
 *
 * Copyright 2015 Marialaura Bancheri
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


// I want to try it with as little dependencies as possible, first
//import java.io.*;
import java.lang.Math;
import richards_classes.*;
public class Richards1d {
	// "Static" keyword defines automatically a globally accessible variable
	// Useful for soil parameters that does not change in time: we can assume them
	// to be the same for the entire program cycle, thus making them accessible
	// to functions without having to pass them as arguments (improves readability)
	
	private  int 	days;
	private	double 	ks;         	// [meter/second]
	private	double 	theta_s;        // Saturated water content
	private	double 	theta_r;        // Residual water content
	private	double 	n;              // For Van Genuchten
	private	double 	m;              // For Van Genuchten
	private	double 	alpha;          // For Van Genuchten

	// Space
	private 	double 	space_bottom;
	private 	double 	space_top;
	private 	int 	NUM_CONTROL_VOLUMES;
	private 	double 	space_delta;
	private 	double[] space_cv_centres;

	// Time
	private 	double 	time_end;
	private 	double 	time_initial;
	private 	double 	time_delta;

	// Time and space
	private	double 	gridvar;
	private	double 	gridvarsq;

	// Cycle variables
	private 	int 	MAXITER;
	private  int 	MAXITER_NEWT;
	private	double 	newton_tolerance;
	
	public Richards1d(int days, double ks, double theta_s, double theta_r, double n, double alpha, double space_bottom, double space_top, int NUM_CONTROL_VOLUMES, double space_delta, double[] space_cv_centres, double time_end, double time_initial, double time_delta, double gridvar, double gridvarsq, int MAXITER, int MAXITER_NEWT, double newton_tolerance){
		this.days = days;
		this.ks = ks;         	
		this.theta_s = theta_s;       
		this.theta_r = theta_r;       
		this.n = n;             
		this.m = 1-1/this.n;             
		this.alpha = alpha;          

		this.space_bottom = space_bottom;
		this.space_top = space_top;
		this.NUM_CONTROL_VOLUMES = NUM_CONTROL_VOLUMES; 
		this.space_delta = space_delta; 			
		this.space_cv_centres = space_cv_centres;

		this.time_end = time_end;            
		this.time_initial = time_initial;
		this.time_delta = time_delta;

		this.gridvar = gridvar;
		this.gridvarsq = gridvarsq;		

		this.MAXITER = MAXITER;
		this.MAXITER_NEWT = MAXITER_NEWT;
		this.newton_tolerance = newton_tolerance;
	}

	public void solve() {

		// Model parameters - SI UNITS
		double 			psi_r				= 0.0;			// Right boundary condition for pressure	
		double 			psi_l				= 0.0;			// Left boundary condition for pressure

		// Time
		int 			time 				= 0;

		// Working variables
		double[] 		psis 				= new double[NUM_CONTROL_VOLUMES];	
		double[] 		dpsis 				= new double[NUM_CONTROL_VOLUMES];	
		double[] 		psis_outer			= new double[NUM_CONTROL_VOLUMES];	// Used for saving \psi of the outer iteration to be plugged inside the inner iteration
		double 			k_r					= 0.0;								
		double 			k_l					= 0.0;
		double[] 		kappas				= new double[NUM_CONTROL_VOLUMES];
		double[] 		thetas				= new double[NUM_CONTROL_VOLUMES];
		double[] 		lowerDiagonal		= new double[NUM_CONTROL_VOLUMES];
		double[] 		mainDiagonal		= new double[NUM_CONTROL_VOLUMES];
		double[]		bb					= new double[NUM_CONTROL_VOLUMES]; // This is the vector which element are those of the principal diagonal and it is used in the Thomas algorithm
		double[] 		upperDiagonal		= new double[NUM_CONTROL_VOLUMES];
		double[]		cc 					= new double[NUM_CONTROL_VOLUMES];
		double[] 		fs					= new double[NUM_CONTROL_VOLUMES];
		double[] 		fks					= new double[NUM_CONTROL_VOLUMES];
		double[] 		dis					= new double[NUM_CONTROL_VOLUMES];
		double[] 		rhss				= new double[NUM_CONTROL_VOLUMES];
		double 			k_p					= 0.0;
		double 			k_m					= 0.0;
		double 			outer_residual		= 0.0;
		double 			inner_residual		= 0.0;

		Thomas thomasAlg = new Thomas();
		PrintTXT print = new PrintTXT();
		SoilParametrization soilPar = new VanGenuchten(n,alpha,theta_r,theta_s,ks);
		JordanDecomposition jordanDecomposition = new JordanDecomposition(soilPar);
		
		// Initial domain conditions
		for(int i = 0; i < NUM_CONTROL_VOLUMES; i++) {
			psis[i] = -space_cv_centres[i];
		}

		////////////////////
		//// MAIN CYCLE ////
		////////////////////
		for(int niter = 0; niter < MAXITER; niter++) {

		    System.out.println("Main cycle iteration:" + niter);		    

		    // Not sure about what that means... Sort of Boundary conditions, but why that form? 
		    if(time <= Math.pow(10,5)) {
		        psi_r = -0.05 + 0.03 * Math.sin(2 * Math.PI * time/Math.pow(10,5));
		    } else if(time > Math.pow(10,5) && time <= 1.8 * Math.pow(10,5)) {
		        psi_r = +0.1;
		    }
		    else {
		        psi_r = -0.05 + 2952.45 * Math.exp(-time / 18204.8);
		    }

			//// KAPPAS ////
		    k_r = soilPar.hydraulicConductivity(psi_r);
		    k_l = soilPar.hydraulicConductivity(psi_l);
		    for(int i = 0; i < NUM_CONTROL_VOLUMES; i++) {      
		    	thetas[i] = soilPar.waterContent(psis[i]);
		    	kappas[i] = soilPar.hydraulicConductivity(psis[i]);
		    }
		   	// "Your time has come"
		   	if(time > time_end) {
		   		break;
		   	}

			/* COEFFICIENT MATRIX IS BUILD BY THREE VECTORS COLLECTING ELEMENTS OF THE THREE DIAGONAL:
		   	 a lower diagonal
		   	 b main diagonal
		   	 c upper diagonal
		   	 RIGHT HAND SIDE */
			for(int i = 0; i < NUM_CONTROL_VOLUMES; i++) {
				if( i == 0 ) {

		            k_p = 0.5*(kappas[i] + kappas[i+1]);
		            k_m = 0.5*(kappas[i] + k_l);
		            lowerDiagonal[i] = 0;
		            mainDiagonal[i] = gridvarsq * (2 * k_m + k_p);
		            upperDiagonal[i] = -k_p * gridvarsq;
		            rhss[i] = thetas[i] + gridvar * (k_p-k_m) + 2 * k_m * gridvarsq * psi_l; 

				} else if(i == NUM_CONTROL_VOLUMES -1) {

		            k_p = 0.5*(kappas[i] + k_r);
		            k_m = 0.5*(kappas[i] + kappas[i-1]);
		            lowerDiagonal[i] = -k_m * gridvarsq;
		            mainDiagonal[i] = gridvarsq * (k_m + 2*k_p);
		            upperDiagonal[i] = 0;
		            rhss[i] = thetas[i] + gridvar * (k_p-k_m) + 2 * k_p * gridvarsq * psi_r; 

				} else {

		            k_p = 0.5*(kappas[i] + kappas[i+1]);
		            k_m = 0.5*(kappas[i] + kappas[i-1]);
		            lowerDiagonal[i] = -k_m * gridvarsq;
		            mainDiagonal[i] = gridvarsq * (k_m + k_p);
		            upperDiagonal[i] = -k_p * gridvarsq;
		            rhss[i] = thetas[i] + gridvar * (k_p-k_m); 
     
				}
			}

			// Initial guess of psis
			for(int i = 0; i < NUM_CONTROL_VOLUMES; i++) {
				psis[i] = Math.min(psis[i], soilPar.getPsiStar());
			}

			//// NESTED NEWTON ////
		    //// OUTER CYCLE, linearizes one of q_{1}, q_{2} ////
		    for(int i = 0; i < MAXITER_NEWT; i++) {
		    	// I have to assign 0 to outer_residual otherwise I will take into account of the previous error
		    	outer_residual = 0.0;
		    	for(int j = 0; j < NUM_CONTROL_VOLUMES; j++) {

		    
		    		fs[j] = soilPar.waterContent(psis[j]) - rhss[j];
			        if(j == 0) {
			            fs[j] = fs[j] + mainDiagonal[j]*psis[j] + upperDiagonal[j]*psis[j+1];
			        }
			        else if(j == NUM_CONTROL_VOLUMES -1) {
			            fs[j] = fs[j] + lowerDiagonal[j]*psis[j-1] + mainDiagonal[j]*psis[j];
			        }
			        else {
			            fs[j] = fs[j] + lowerDiagonal[j]*psis[j-1] + mainDiagonal[j]*psis[j] + upperDiagonal[j]*psis[j+1];
			        }
			        outer_residual += fs[j]*fs[j];
		    	}
		    	outer_residual = Math.pow(outer_residual,0.5);

			    System.out.println("Outer iteration " + i + " with residual " +  outer_residual);    

    			if(outer_residual < newton_tolerance) {
    				break;
    			}

    			psis_outer = psis.clone();

				// Initial guess of psis
				for(int j = 0; j < NUM_CONTROL_VOLUMES; j++) {
					psis[j] = Math.max(psis[j], soilPar.getPsiStar());
				}

			    //// INNER CYCLE ////
				for(int j = 0; j < MAXITER_NEWT; j++) {
					// I have to assign 0 to inner_residual otherwise I will take into account of the previous error
					inner_residual = 0.0; 
					for(int l=0; l < NUM_CONTROL_VOLUMES; l++) {

				        //fks[l] = theta1(psis[l]) - (theta2(psis_outer[l]) + dtheta2(psis_outer[l])*(psis[l] - psis_outer[l])) - rhss[l];
						fks[l] = jordanDecomposition.waterContent1(psis[l]) - (jordanDecomposition.waterContent2(psis_outer[l]) + jordanDecomposition.dWaterContent2(psis_outer[l])*(psis[l] - psis_outer[l])) - rhss[l];
			            if(l == 0) {
			                fks[l] = fks[l] + mainDiagonal[l]*psis[l] + upperDiagonal[l]*psis[l+1];
			            }
			            else if(l == NUM_CONTROL_VOLUMES -1) {
			                fks[l] = fks[l] + lowerDiagonal[l]*psis[l-1] + mainDiagonal[l]*psis[l];
			            }
			            else {
			                fks[l] = fks[l] + lowerDiagonal[l]*psis[l-1] + mainDiagonal[l]*psis[l] + upperDiagonal[l]*psis[l+1];
			            }
			            dis[l] = jordanDecomposition.dWaterContent1(psis[l]) - jordanDecomposition.dWaterContent2(psis_outer[l]);
				        inner_residual += fks[l]*fks[l];

				    }
				    inner_residual = Math.pow(inner_residual,0.5);

			    	System.out.println("Inner iteration " + j + " with residual " +  inner_residual);    

			        if(inner_residual < newton_tolerance) {
			            break;
			        }

				    //// CONJGRAD ////
			        // Attention: the main diagonal of the coefficient matrix must not change!! The same for the upper diagonal
			        
			        bb = mainDiagonal.clone();
			        cc = upperDiagonal.clone();
				    for(int y = 0; y < NUM_CONTROL_VOLUMES; y++) {
				    	bb[y] += dis[y];
				    }
				    thomasAlg.set(cc,bb,lowerDiagonal,fks);
				    dpsis = thomasAlg.solver();
				    //dpsis = thomas(lowerDiagonal,bb,cc,fks);
				    
				    //// PSIS UPDATE ////
			        for(int s = 0; s < NUM_CONTROL_VOLUMES; s++) {
			        	psis[s] = psis[s] - dpsis[s];
			        }

				} //// INNER CYCLE END ////
			} //// OUTER CYCLE END ////
		    //// Print value of psis with PrintTXT class  ////
		    
			print.setValueFirstVector(psis);
			print.setValueSecondVector(space_cv_centres);
			
			print.PrintTwoVectors("PsiTestJoradan_"+time+"_s.txt", "Psi values at time: "+time+" seconds", "Psi[m] Depth[m]");
			print.setValueFirstVector(dpsis);
			print.PrintTwoVectors("dPsiTestJoradan_"+time+"_s.txt", "Psi values at time: "+time+" seconds", "Psi[m] Depth[m]");
			// "The show must go on"
	   	time += time_delta;

		} //// MAIN CYCLE END ////
		//print(psis);
	} //// MAIN END ////

	////////////////////
	//// METHODS ////
	////////////////////



	// UTILITY METHODS

	private static void print(double[] arr) {
		System.out.println(java.util.Arrays.toString(arr));
	}
	private static void print(double num) {
		System.out.println(num);
	}
}