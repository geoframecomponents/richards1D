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
import java.lang.Math;

public class Richards1d {
	// "Static" keyword defines automatically a globally accessible variable
	// Useful for soil parameters that does not change in time: we can assume them
	// to be the same for the entire program cycle, thus making them accessible
	// to functions without having to pass them as arguments (improves readability)
	static  int 	days				= 24*3600;
	static	double 	ks 					= 0.062/days;  	// [meter/second]
	static	double 	theta_s				= 0.41;         // Saturated water content
	static	double 	theta_r				= 0.095;        // Residual water content
	static	double 	n					= 1.31;         // For Van Genuchten
	static	double 	m					= 1-1/n;        // For Van Genuchten
	static	double 	alpha				= 1.9;          // For Van Genuchten
	static	double 	psi_crit			= (-1.0/alpha)*Math.pow((n-1.0)/n,1.0/n);  // Where \frac{\partial\psi(\theta)}{\theta}=0


	// Space
	static 	double 	space_bottom		= 0.0;
	static 	double 	space_top			= 2.0;
	static 	int 	NUM_CONTROL_VOLUMES	= 100; 
	static 	double 	space_delta			= (space_top - space_bottom) / NUM_CONTROL_VOLUMES; 			// delta
	static 	double[] space_cv_centres	= seq(space_bottom + space_delta / 2,space_top - space_delta / 2,NUM_CONTROL_VOLUMES); // Centres of the "control volumes"

	// Time
	static 	double 	time_end 			= 3*Math.pow(10,8);            
	static 	double 	time_initial 		= 0.0;
	static 	double 	time_delta 			= 1000000.0;

	// Time and space
	static	double 	gridvar				= time_delta / space_delta;
	static	double 	gridvarsq			= time_delta / Math.pow(space_delta,2);		

	// Cycle variables
	static 	int 	MAXITER 			= 1000;
	static  int 	MAXITER_NEWT 		= 100000;
	static	double 	newton_tolerance	= Math.pow(10,-12);


	public static void main(String[] args) {

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
		double[] 		a 					= new double[NUM_CONTROL_VOLUMES];
		double[] 		b 					= new double[NUM_CONTROL_VOLUMES];
		double[] 		c 					= new double[NUM_CONTROL_VOLUMES];
		double[] 		fs					= new double[NUM_CONTROL_VOLUMES];
		double[] 		fks					= new double[NUM_CONTROL_VOLUMES];
		double[] 		dis					= new double[NUM_CONTROL_VOLUMES];
		double[] 		rhss				= new double[NUM_CONTROL_VOLUMES];
		double 			k_p					= 0.0;
		double 			k_m					= 0.0;
		double 			outer_residual		= 0.0;
		double 			inner_residual		= 0.0;



		// Initial domain conditions
		for(int i = 0; i < NUM_CONTROL_VOLUMES; i++) {
			psis[i] = -space_cv_centres[i];
		}

		////////////////////
		//// MAIN CYCLE ////
		////////////////////
		System.out.println("modifica");
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
		    k_r = kappa(psi_r); 
		    k_l = kappa(psi_l);
		    for(int i = 0; i < NUM_CONTROL_VOLUMES; i++) {
			    thetas[i] = thetaf(psis[i]);
			    kappas[i] = kappa(psis[i]);           
			}

		   	// "Your time has come"
		   	if(time > time_end) {
		   		break;
		   	}

			//// RIGHT HAND SIDE ////
			for(int i = 0; i < NUM_CONTROL_VOLUMES; i++) {
				if( i == 0 ) {

		            k_p = 0.5*(kappas[i] + kappas[i+1]);
		            k_m = 0.5*(kappas[i] + k_l);
		            a[i] = 0;
		            b[i] = gridvarsq * (2 * k_m + k_p);
		            c[i] = -k_p * gridvarsq;
		            rhss[i] = thetas[i] + gridvar * (k_p-k_m) + 2 * k_m * gridvarsq * psi_l; 

				} else if(i == NUM_CONTROL_VOLUMES -1) {

		            k_p = 0.5*(kappas[i] + k_r);
		            k_m = 0.5*(kappas[i] + kappas[i-1]);
		            a[i] = -k_m * gridvar;
		            b[i] = gridvarsq * (k_m + 2*k_p);
		            c[i] = 0;
		            rhss[i] = thetas[i] + gridvar * (k_p-k_m) + 2 * k_p * gridvarsq * psi_r; 

				} else {

		            k_p = 0.5*(kappas[i] + kappas[i+1]);
		            k_m = 0.5*(kappas[i] + kappas[i-1]);
		            a[i] = -k_m * gridvarsq;
		            b[i] = gridvarsq * (k_m + k_p);
		            c[i] = -k_p * gridvarsq;
		            rhss[i] = thetas[i] + gridvar * (k_p-k_m); 
     
				}
			}


			// Initial guess of psis
			for(int i = 0; i < NUM_CONTROL_VOLUMES; i++) {
				psis[i] = Math.min(psis[i], psi_crit);
			}

			//// NESTED NEWTON ////
		    //// OUTER CYCLE, linearizes one of q_{1}, q_{2} ////
		    for(int i = 0; i < MAXITER_NEWT; i++) {

		    	for(int j = 0; j < NUM_CONTROL_VOLUMES; j++) {

		    		fs[j] = thetaf(psis[j]) - rhss[j];        
			        if(j == 0) {
			            fs[j] = fs[j] + b[j]*psis[j] + c[j]*psis[j+1];
			        }
			        else if(j == NUM_CONTROL_VOLUMES -1) {
			            fs[j] = fs[j] + a[j]*psis[j-1] + b[j]*psis[j];
			        }
			        else {
			            fs[j] = fs[j] + a[j]*psis[j-1] + b[j]*psis[j] + c[j]*psis[j+1];
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
					psis[j] = Math.max(psis[j], psi_crit);
				}

			    //// INNER CYCLE ////
				for(int j = 0; j < MAXITER_NEWT; j++) {

					for(int l=0; l < NUM_CONTROL_VOLUMES; l++) {

				        fks[l] = theta1(psis[l]) - (theta2(psis_outer[l]) + dtheta2(psis_outer[l])*(psis[l] - psis_outer[l])) - rhss[l];

			            if(l == 0) {
			                fks[l] = fks[l] + b[l]*psis[l] + c[l]*psis[l+1];
			            }
			            else if(l == NUM_CONTROL_VOLUMES -1) {
			                fks[l] = fks[l] + a[l]*psis[l-1] + b[l]*psis[l];
			            }
			            else {
			                fks[l] = fks[l] + a[l]*psis[l-1] + b[l]*psis[l] + c[l]*psis[l+1];
			            }
			            dis[l] = dtheta1(psis[l]) - dtheta2(psis_outer[l]);
				        inner_residual += fks[l]*fks[l];

				    }
				    inner_residual = Math.pow(inner_residual,0.5);

			    	System.out.println("Inner iteration " + j + " with residual " +  inner_residual);    

			        if(inner_residual < newton_tolerance) {
			            break;
			        }

				    //// CONJGRAD ////
				    for(int y = 0; y < NUM_CONTROL_VOLUMES; y++) {
				    	b[y] += dis[y];
				    }
			        dpsis = thomas(a,b,c,fks);

				    //// PSIS UPDATE ////
			        for(int s = 0; s < NUM_CONTROL_VOLUMES; s++) {
			        	psis[s] = psis[s] - dpsis[s];
			        }

				} //// INNER CYCLE END ////
			} //// OUTER CYCLE END ////

		// "The show must go on"
	   	time += time_delta;

		} //// MAIN CYCLE END ////
		print(psis);
	} //// MAIN END ////

	////////////////////
	//// METHODS ////
	////////////////////

	/**
	 * Returns a K, given the soil parameters
	 *
	 * @param  psi  	 
	 * @return 			a real number representing the control volume's K
	 */	
	public static double kappa(double psi) {

		double kappa;
		double saturation;

		saturation = (thetaf(psi) - theta_r) / (theta_s - theta_r); 
		kappa = ks * Math.pow(saturation, 0.5 ) * Math.pow(1.0 - Math.pow(1.0 - Math.pow(saturation, 1.0/m), m), 2.0);
		return kappa;
	}

	/**
	 * Returns the thetaf term
	 *
	 * @param  psi  	 
	 * @return 			a real number representing the \theta_{f}
	 */	
	public static double thetaf(double psi) {

		double theta_f;

		if(psi <= 0) {
		    theta_f = theta_r + (theta_s - theta_r) / Math.pow(1.0 + Math.pow(Math.abs(alpha*psi), n), m);
		} else {
		    theta_f = theta_s;
		}
		return theta_f;
	}

	/**
	 * Returns the dtheta term
	 *
	 * @param  psi   
	 * @return 			a real number representing the $\delta(\theta)$
	 */	
	public static double dtheta(double psi) {

		double dtheta;

		if (psi <= 0) {
			// Unsaturated case
		    dtheta = alpha*n*m*(theta_s - theta_r) / Math.pow(1.0 + Math.pow(Math.abs(alpha*psi), n), m + 1.0)*Math.pow(Math.abs(alpha*psi), n - 1.0);
		} else {
		    dtheta = 0;
		}
		return dtheta;
	}

	/**
	 * Returns the dtheta1 term
	 *
	 * @param  psi   
	 * @return 			a real number representing the $\delta(\theta)1$
	 */	
	public static double dtheta1(double psi) {

		double dtheta1;
		
		if (psi <= psi_crit) {
		    // left of critical value, take the original derivative
		    dtheta1 = dtheta(psi);
		}
		else {
		    // on the right of the critical value, keep the maximum derivative
		    dtheta1 = dtheta(psi_crit);
		}

		return dtheta1;

	}

	/**
	 * Returns the dtheta2 term
	 *
	 * @param  psi   
	 * @return 			a real number representing the $\delta(\theta)2$
	 */	
	public static double dtheta2(double psi) {

		double dtheta2;
		
		dtheta2 = dtheta1(psi) - dtheta(psi);

		return dtheta2;
	}

	/**
	 * Returns the $\theta_{1}$ term
	 *
	 * @param  psi   
	 * @return 			a real number representing the $\theta_{1}$
	 */	
	public static double theta1(double psi) {

		double theta1;
		
		if(psi <= psi_crit) {
			theta1 = thetaf(psi);
		} else {
			theta1 = thetaf(psi_crit) + dtheta(psi_crit)*(psi - psi_crit);
		}

		return theta1;
	}


	/**
	 * Returns the $\theta_{2}$ term
	 *
	 * @param  psi   
	 * @return 			a real number representing the $\theta_{2}$
	 */	
	public static double theta2(double psi) {

		double theta2;
		
		theta2 = theta1(psi) - thetaf(psi);

		return theta2;	
	}


	/**
	 * Linear system solver with Thomas method
	 *
	 * @param  a   
	 * @param  b   
	 * @param  c   
	 * @param  d   
	 * @return solution the vector of solutions
	 */	

	public static double[] thomas(double[] a, double[] b, double[] c, double[] d) {

		int DIM = d.length;
		double gamma = 0.0;
		double[] solution = new double[DIM];

		c[0] = c[0] / b[0];
		d[0] = d[0] / b[0];

		for(int i = 1; i < DIM; i++) {
			gamma = 1 / (b[i] - c[i-1]*a[i]);
			c[i] = c[i]*gamma;
			d[i] = (d[i] - a[i]*d[i-1])*gamma;
		}

		solution[DIM-1] = d[DIM-1];

		for(int i = DIM-2; i > -1; i--) {
			solution[i] = d[i] - c[i]*solution[i+1];
		}

		return solution;

	}

	// UTILITY METHODS

	/**
	 * Returns a vector (array of doubles) of equally spaced numbers between a given 
	 * lower and upper bound
	 * <p>
	 * The lower and upper bounds need not be in growing order: in 
	 * case min>max, the method generates a sequence of decreasing 
	 * real numbers
	 *
	 * @param  max  	a real number, the value of the lower bound of the sequence
	 * @param  max  	a real number, the value of the upper bound of the sequence
	 * @param  points  	an integer, the number of equally spaced points in the sequence 
	 * @return 			a vector of equally spaced real numbers
	 */	
	public static double[] seq(double min, double max, int points) {  
	    double[] sequence = new double[points];  
	    for (int i = 0; i < points; i++){  
	        sequence[i] = min + i * (max - min) / (points - 1);  
	    }  
	    return sequence;  
	}	

	private static void print(double[] arr) {
		System.out.println(java.util.Arrays.toString(arr));
	}
	private static void print(double num) {
		System.out.println(num);
	}
}