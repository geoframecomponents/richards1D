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
import java.lang.Math

public class Richards1d {
	public static void main(String[] args) {
		// Model parameters - SI UNITS

		// "Static" keyword defines automatically a globally accessible variable
		// Useful for soil parameters that does not change in time: we can assume them
		// to be the same for the entire program cycle, thus making them accessible
		// to functions without having to pass them as arguments (improves readability)
		int 	static 	days				= 24*3600;
		double 	static	ks 					= 0.062/day;  	//[meter/second]
		double 	static	theta_s				= 0.41;         //[-] saturated water content
		double 	static	theta_r				= 0.095;        //[-] residuel water content
		double 	static	n					= 1.31;         // For Van Genuchten
		double 	static	m					= 1-1/n;        // For Van Genuchten
		double 	static	alpha				= 1.9;          // For Van Genuchten
		double 	static	psi_crit			= Math.pow(-1/alpha * (n-1)/n , 1/n);  // Where \frac{\partial\psi(\theta)}{\theta}=0
		double 			psi_r				= 0;			// Right boundary condition for pressure	
		double 			psi_l				= 0;			// Left boundary condition for pressure

		// Space
		double 	static 	space_bottom		= 0;
		double 	static 	space_top			= 2;
		double 	static 	NUM_CONTROL_VOLUMES	= 100; 
		double 	static 	space_delta			= (space_top - space_bottom) / NUM_CONTROL_VOLUMES; 			// delta
		double 	static 	space_cv_centres[]	= seq(space_bottom + space_delta / 2,space_top - space_delta / 2,NUM_CONTROL_VOLUMES); // Centres of the "control volumes"

		// Time
		double 	static 	time_end 			= Math.exp(3,5);            
		double 	static 	time_initial 		= 0.0;
		int 	static 	time_delta 			= 1000;
		int 			time 				= 0;

		// Time and space
		double 	static	gridvar				= time_delta / space_delta;
		double 	static	gridvarsq			= Math.pow((time_delta / space_delta),2);		

		// Working variables
		double[] 		psis 				= new double[NUM_CONTROL_VOLUMES];	
		double[] 		dpsis 				= new double[NUM_CONTROL_VOLUMES];	
		double[] 		psis_outer			= new double[NUM_CONTROL_VOLUMES];	// Used for saving \psi of the outer iteration to be plugged inside the inner iteration
		double 			k_r					= 0;								
		double 			k_l					= 0;
		double[] 		kappas				= new double[NUM_CONTROL_VOLUMES];
		double[] 		a 					= new double[NUM_CONTROL_VOLUMES];
		double[] 		b 					= new double[NUM_CONTROL_VOLUMES];
		double[] 		c 					= new double[NUM_CONTROL_VOLUMES];
		double[] 		fs					= new double[NUM_CONTROL_VOLUMES];
		double[] 		fks					= new double[NUM_CONTROL_VOLUMES];
		double[] 		deltas				= new double[NUM_CONTROL_VOLUMES];
		double[] 		rhss				= new double[NUM_CONTROL_VOLUMES];
		double 			k_p					= 0;
		double 			k_m					= 0;
		double 			outer_residual		= 0;
		double 			inner_residual		= 0;

		// Cycle variables
		int    	static 	MAXITER 		= 100000;
		int 	static  MAXITER_NEWT = 100;
		int 	static	newton_tolerance	= Math.pow(10,-12)

		// Initial domain conditions
		for(int i = 0; i <= NUM_CONTROL_VOLUMES; i++) {
			psi[i] = -y[i];
		}

		////////////////////
		//// MAIN CYCLE ////
		////////////////////
		for(int niter = 0; niter < MAXITER; niter++) {
		    if(time_initial + time_delta > time_end) {
		        time_delta = time_end-time_initial;
		    }

		    // TODO
		    // Not sure about what that means... Sort of Boundary conditions, but why that form? 
		    if(time <= Math.pow(10,5)) {
		        psi_r = -0.05 + 0.03 * Math.sin(2*Math.PI * time/Math.pow(10,5));
		    } else if(time > Math.pow(10,5) && time <= 1.8 * Math.pow(10,5) {
		        psi_r = +0.1;
		    }
		    else {
		        psi_r = -0.05 + 2952.45 * Math.exp(-time / 18204.8);
		    }

			//// KAPPAS ////
		    k_r = kappa(psi_r); 
		    k_l = kappa(psi_l);
		    for(int i = 0; i < MAXITER; i++) {
			    theta[i] = thetaf(theta_s, theta_r, alpha, psi, n, m, psis[i]);
			    kappas[i] = kappa(alpha, theta_2, theta_r, psi[i]);           
			}

			//// RIGHT HAND SIDE ////
			for(int i = 0; i < MAXITER; i++) {
				if( i==1 ) {

		            k_p = 0.5*(kappas[i] + kappas[i+1]);
		            k_m = 0.5*(kappas[i] + k_l);
		            a[i] = 0;
		            b[i] = gridvarsq * (2 * k_m + k_p);
		            c[i] = -k_p * gridvarsq;
		            rhs[i] = thetas[i] + gridvar * (k_p-k_m) + 2 * k_m * gridvarsq * psi_l; 

				} else if(i == MAXITER) {

		            k_p = 0.5*(kappas[i] + k_r);
		            k_m = 0.5*(kappas[i] + k_l);
		            a[i] = -k_m * gridvar;
		            b[i] = gridvarsq * (k_m + 2*k_p);
		            c[i] = 0;
		            rhs[i] = thetas[i] + gridvar * (k_p-k_m) + 2 * k_p * gridvarsq * psi_r; 

				} else {

		            k_p = 0.5*(kappas[i] + kappas[i+1]);
		            k_m = 0.5*(kappas[i] + kappas[i-1]);
		            a[i] = -k_m * gridvarsq;
		            b[i] = gridvarsq * (k_m + k_p);
		            c[i] = 0;
		            rhs[i] = thetas[i] + gridvar * (k_p-k_m); 
     
				}
			}

			// Initial guess of psis
			for(int i = 0; i < NUM_CONTROL_VOLUMES; i++) {
				psis[i] = Math.min(psis[i], psi_crit);
			}

			//// NESTED NEWTON ////
		    //// OUTER CYCLE, linearizes one of q_{1}, q_{2} ////
		    for(int i = 0; i < MAXITER_NEWT; i++) {

		    	for(int ii=0; ii < MAXITER; ii++) {

		    		fs[i] = thetaf(psi[i]) - rhs[i];        

			        if(i == 1) {
			            f[i] = f[i] + b[i]*psi[i] + c[i]*psi[i+1];
			        }
			        else if(i == MAXITER) {
			            f[i] = f[i] + a[i]*psi[i-1] + b[i]*psi[i];
			        }
			        else {
			            f[i] = f[i] + a[i]*psi[i-1] + b[i]*psi[i] + c[i]*psi[i+1];
			        }
			        outer_residual += f[i];
		    	}
		    
			    System.out.println('Outer iteration ' + i ' with residual ' +  outer_residual));    
    			if(outer_residual < newton_tolerance) {
    				break;
    			}

    			System.arraycopy(psis, 0, psis_outer, 0);

				// Initial guess of psis
				for(int ii = 0; ii < NUM_CONTROL_VOLUMES; ii++) {
					psis[ii] = Math.max(psis[ii], psi_crit);
				}

			    //// INNER CYCLE ////
				for(int j = 0; j < NUM_CONTROL_VOLUMES; j++) {

					for(int l=0; l < MAXITER; l++) {

				        fks(l)=theta1(psis[l]) - (theta2(psis_outer[l]) + dtheta2(psis_outer[l])*(psis[l] - psis_outer[l])) - rhs[l];
			            if(l == 1) {
			                fks[i] = fks[i] + b[l]*psis[l] + c[l]*psis[l+1];
			            }
			            else if(l == MAXITER) {
			                fks[l]=fks[l] + a[l]*psi[l-1] + b[l]*psis[l];
			            }
			            else {
			                fks[l]=fks[i] + a[l]*psis[l-1] + b[i]*psis[l] + c[l]*psis[l+1];
			            }
			            deltas[l] = dtheta1(psis[l]) - dtheta2(psis_outer[l]);
				        inner_residual += fks[l];

			    		System.out.println('Inner iteration ' + l + ' with residual ' +  inner_residual));    
				        if(inner_residual < newton_tolerance) {
				            break;
				        }
					    // TODO:_ conjgrad (MUST RETURN AN ARRAY)
					    //// CONJGRAD ////
				        dpsis = Thomas(a,b+deltas,c,fks);
					    //// PSIS UPDATE ////
				        for(int s = 0; s < NUM_CONTROL_VOLUMES; s++) {
				        	psis[s] = psis[s] - dpsis[s];
				        }

					}

			// "The show must go on"
		   	time += time_delta;

		   	// "Your time has come"
		   	if(time > time_end) {
		   		break;
		   	}
		}



		// TODO: dtheta1, dtheta2, theta1, theta2
		////////////////////
		//// METHODS ////
		////////////////////

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

		/**
		 * Returns a K, given the soil parameters
		 *
		 * @param  psi  	 
		 * @return 			a real number representing the control volume's K
		 */	
		public static double k(double psi) {

			double kappa;
			double saturation;

			saturation = (thetaf(psi) - theta_r) / (theta_s - theta_r); 
			kappa = ks * sqrt(saturation) * Math.pow(1 - (1 - Math.pow(Math.pow(sat,1 / m)), m), 2);

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
			    theta_f = theta_r + (theta_s - theta_r) / ( Math.pow(1+Math.pow(Math.abs(alpha*psi), n), m);
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
			    dtheta = alpha*n*m*(thetas - thetar) / ((1 + Math.pow(Math.pow(Math.abs(alpha*psi), n), (m + 1))*Math.pow(Math.abs(alpha*psi), (n - 1));
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

			theta1 = dtheta1(psi) - dtheta(psi);

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


}