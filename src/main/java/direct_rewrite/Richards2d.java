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

public class Richards2d {
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
	static	double 	psi_crit			= (-1/alpha)*Math.pow((n-1)/n,1/n);  // Where \frac{\partial\psi(\theta)}{\theta}=0


	// Space
	static 	double 	space_bottom		= 0.0; // x_{l}
	static 	double 	space_top			= 2.0; // x_{u}
	static 	double 	space_left			= -1.0;// y_{l} 
	static 	double 	space_right			= 1.0; // y_{r}
	static 	int 	NUM_CV_X			= 100; 
	static 	int 	NUM_CV_Y			= 50; 
	static 	double 	space_mesh_delta_x	= (space_top - space_bottom) / NUM_CV_X;
	static 	double 	space_mesh_delta_y	= (space_right - space_left) / NUM_CV_Y;	
	static 	double[] space_cv_centres_x	= seq(space_bottom + space_mesh_delta_x / 2,space_top - space_mesh_delta_x / 2,NUM_CV_X);
	static 	double[] space_cv_centres_y	= seq(space_left + space_mesh_delta_y / 2,space_right - space_mesh_delta_y / 2,NUM_CV_Y);

	// Time
	static 	double 	time_end 			= 3*Math.pow(10,5);            
	static 	double 	time_initial 		= 0.0;
	static 	double 	time_delta 			= 1000000.0;

	// Time and space
	static	double 	gridvar				= time_delta / space_mesh_delta_x;
	static	double 	gridvarsq			= time_delta / Math.pow(space_mesh_delta_x,2);		

	// Cycle variables
	static 	int 	MAXITER 			= 100000;
	static  int 	MAXITER_NEWT 		= 100;
	static	double 	newton_tolerance	= Math.pow(10,-6);


	// i > iteraton over X coord.
	// j > iteration over y coord.
	public static void main(String[] args) {

		// Model parameters - SI UNITS
		double 			psi_r				= 0.0;			// Right boundary condition for pressure	
		double 			psi_l				= 0.0;			// Left boundary condition for pressure

		// Time initial conditions
		int 			time 				= 0;

		// Working variables
		double[][] 		psis 				= new double[NUM_CV_X][NUM_CV_Y];	
		double[][] 		dpsis 				= new double[NUM_CV_X][NUM_CV_Y];	
		double[][] 		dis 				= new double[NUM_CV_X][NUM_CV_Y];	
		double[][] 		mpsis 				= new double[NUM_CV_X][NUM_CV_Y];	
		double[][] 		psis_outer			= new double[NUM_CV_X][NUM_CV_Y];	// Used for saving \psi of the outer iteration to be plugged inside the inner iteration
		double 			k_r					= 0.0;								
		double 			k_l					= 0.0;
		double[][] 		kappas				= new double[NUM_CV_X][NUM_CV_Y];
		double[][] 		thetas				= new double[NUM_CV_X][NUM_CV_Y];
		double[][] 		a 					= new double[NUM_CV_X][NUM_CV_Y];
		double[][] 		b 					= new double[NUM_CV_X][NUM_CV_Y];
		double[][] 		c 					= new double[NUM_CV_X][NUM_CV_Y];
		double[][] 		fs					= new double[NUM_CV_X][NUM_CV_Y];
		double[][] 		fks					= new double[NUM_CV_X][NUM_CV_Y];
		double[][] 		deltas				= new double[NUM_CV_X][NUM_CV_Y];
		double[][] 		rhss				= new double[NUM_CV_X][NUM_CV_Y];
		double 			k_p					= 0.0;
		double 			k_m					= 0.0;
		double 			outer_residual		= 0.0;
		double 			inner_residual		= 0.0;

		// Domain initial conditions
		for(int i=0; i < NUM_CV_X; i++) {
		    for(int j; j < NUM_CV_Y; j++) {
		        psi[i][j] = -x(i); // Hydrostatic pressure  
		    }
		}

		for(int main_iter = 0; main_inter < MAXITER; main_inter++) {
		    if(time + time_delta > time_end) {
		        time_delta = time_end - time;
		    }

		    // right-upper boundary condition for pressure
		    if(time <= Math.pow(10,5)) {
		        psi_r = -0.05 + 0.03 * Math.sin(2*Math.PI * time/Math.pow(10,5));
		    } else if(time > Math.pow(10,5) && time <= 1.8 * Math.pow(10,5)) {
		        psi_r = +0.1;
		    }
		    else {
		        psi_r = -0.05 + 2952.45 * Math.exp(-time / 18204.8);
		    }
		    // left-lower boundary condition for pressure
		    psi_l = 0.0;

		    k_r = kappa(psi_r); 
		    k_l = kappa(psi_l);
		    for(int i = 0; i < NUM_CV_X; i++) {
		    	for(int j = 0; j < NUM_CV_X; j++) {
				    thetas[i][j] = thetaf(psis[i][j]);
				    kappas[i][j] = kappa(psis[i][j]);           
				}
			}

		    if(time >= time_end) {
		        break;
		    }

		    for(int i = 0; i < NUM_CV_X; i++) {
		    	for(int j = 0; j < NUM_CV_X; j++) {

		            if(i==1) {
		                k_p = 0.5*(kappas[i][j]+kappas[i+1][j]);
		                k_m = 0.5*(kappas[i][j]+k_l);
		                // Dirichlet (pressure) boundary condition
		                rhss[i][j] = theta[i][j] + gridvar*(k_p-k_m) + 2*k_m*gridvarsq*psi_l;
		            } else if(i == NUM_CV_X -1) {
		                k_p = 0.5*(kappas[i][j]] + k_r);
		                k_m = 0.5*(kappas[i][j] + kappas[i-1][j]);
		                // Dirichlet (pressure) boundary condition
		                rhss[i][j] = theta[i][j] + gridvar*(k_p-k_m) + 2*k_p*gridvarsq*psi_r;
		            } else {
		            	k_p = 0.5*(kappas[i][j] + kappas[i+1][j])
		                k_m = 0.5*(kappas[i][j] + kappas[i-1][j])
		                rhss[i,j] = thetas[i][j] + gridvar*(k_p-k_m)
		            }
				}
			}

			// Initial guess for psi
		    for(int i = 0; i < NUM_CV_X; i++) {
		    	for(int j = 0; j < NUM_CV_X; j++) {			
		    		psis[i][j] = Math.min(psis[i][j], psi_crit);
		    	}
		    }



		    ///// NEWTON ////
		    for(int newt_iter = 0; newt_iter < MAXITER_NEWT; newt_iter++) {
		    	// A double gets assigned by default the 0.0 value, but if I want to have the array filled with
		    	// some other value, that could be useful:
				// Fill each row
				for (double[] row: dis)
				    Arrays.fill(row, 0.0);
				mpsis = matop(dis,psis);
			    for(int i = 0; i < NUM_CV_X; i++) {
			    	for(int j = 0; j < NUM_CV_X; j++) {			
			    		fs[i][j] = thetaf(psis[i][j]) + mpsis[i][j] - rhss[i][j];
			    		outer_residual += Math.pow(Math.abs(fs[i][j]*fs[i][j]),1.0/2.0);
			    	}
			    }
			    if(outer_residual < newton_tolerance) {
			    	break;
			    }

			    psis_outer = psis.clone();
			    for(int i = 0; i < NUM_CV_X; i++) {
			    	for(int j = 0; j < NUM_CV_X; j++) {			
			    		psis[i][j] = Math.max(psis[i][j],psi_crit);
			    	}
			    }

			    for(int inner_iter = 0; inner_iter < MAXITER_NEWT; inner_iter++) {
					for (double[] row: dis)
					    Arrays.fill(row, 0.0);
					mpsis = matop(dis,psis);
				    for(int i = 0; i < NUM_CV_X; i++) {
				    	for(int j = 0; j < NUM_CV_X; j++) {			
		                    fks[i][j]=theta1(psi[i][j]) - (theta2(psis_outer[i][j]) + dtheta2(psis_outer[i][j])*(psi[i][j]-psis_outer[i][j])) + mpsis[i][j]-rhs[i][j];
		                    dis[i][j]=dtheta1(psi[i][j]) - dtheta2(psis_outer[i][j]);
		                    inner_residual += Math.pow(Math.abs(fks[i][j]*fks[i][j]),1.0/2.0);
				    	}
				    }
				    if(inner_residual < newton_tolerance) {
				    	break;
				    }
				    dpsis = cgop(fks);
				    for(int i = 0; i < NUM_CV_X; i++) {
				    	for(int j = 0; j < NUM_CV_X; j++) {			
				    		psis[i][j] = psis[i][j] - dpsis[i][j];
				    	}
				    }
				}			 
		    }
		        
		    time += time_delta;

		} //// MAIN CYCLE ENDED ////
	} /// MAIN METHOD ENDED ///


	////////////////////
	// OTHER METHODS //
	///////////////////

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
		kappa = ks * Math.pow(saturation, 1/2 ) * Math.pow(1 - Math.pow(1 - Math.pow(saturation, 1/m), m), 2);
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
		    theta_f = theta_r + (theta_s - theta_r) / Math.pow(1 + Math.pow(Math.abs(alpha*psi), n), m);
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
		    dtheta = alpha*n*m*(theta_s - theta_r) / Math.pow(1 + Math.pow(Math.abs(alpha*psi), n), m+1)*Math.pow(Math.abs(alpha*psi), n-1);
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

		for(int i = 0; i < DIM; i++) {
			gamma = 1 / (b[i] - c[i]*a[i]);
			c[i] = c[i]*gamma;
			d[i] = (d[i] - a[i]*d[i])*gamma;
		}


		solution[DIM-1] = d[DIM-1];

		for(int i = DIM-1; i > 0; i--) {
			solution[i] = d[i] - c[i]*solution[i-1];
		}

		return solution;

	}


	%matrix-free conjugate method for the solution of A*x=b with A symmetric positive definite
	%Input: (we don't need the matrix as input)
	%       b = right hand side
	%output: 
	%       x = solution vector
	%the user has to provide a function "matop2D(x)" wich computes the matrix-vector product A*x

	function x=CGop2D(b)
	x = b;                %initial guess
	N = numel(b);         %get the problem size 
	r = b-matop2D(x);     %residual
	p = r;                %initial search direction
	alpha = sum(sum(r.*r));  %square of r (scalar product) %or r'*r; 
	for k=1:4*N           %2* or 4*N to have more precision
	    if (sqrt(alpha)<1e-14)   %nominal machine operation is 1e-16 so we use something smaller
	        %norm of the residual is sufficently small
	        break
	    end
	    v = matop2D(p);      
	    lambda = alpha/sum(sum(p.*v));
	    x = x + lambda*p;   
	    r = r - lambda*v;
	    alphak = alpha;
	    alpha = sum(sum(r.*r));
	    p = r + alpha/alphak*p;   
	end

	function Apsi = matop2D(psi)
	global di dx dy dt K IMAX JMAX KL KR
	% diagonal part including the derivatives of the nonlinear functions
	Apsi = di.*psi;
	% linear part
	public static double[] matop(double[][] di, double[][] psis) {
		double[][] apsi ==
	    for(int i = 0; i < NUM_CV_X; i++) {
	    	for(int j = 0; j < NUM_CV_X; j++) {
	    		if(i == 1) {

	    		}
	    	}
	    }	
	        if(i==1)
	            Kp = 0.5*(K[i][j]+K(i+1,j));
	            Km = 0.5*(K[i][j]+KL);
	            Apsi[i][j] = Apsi[i][j]-dt/dx^2*( Kp*(psi(i+1,j)-psi[i][j])...
	                                           - 2*Km*(psi[i][j]-0) );
	        elseif(i==IMAX)
	            Kp = 0.5*(K[i][j] + KR);
	            Km = 0.5*(K[i][j] + K(i-1,j));
	            Apsi[i][j] = Apsi[i][j]-dt/dx^2*( 2*Kp*(0-psi[i][j])...
	                                           -Km*(psi[i][j]-psi(i-1,j)) );
	        else
	            Kp = 0.5*(K[i][j] + K(i+1,j));
	            Km = 0.5*(K[i][j] + K(i-1,j));
	            Apsi[i][j] = Apsi[i][j]-dt/dx^2*( Kp*(psi(i+1,j)-psi[i][j])...
	                                           -Km*(psi[i][j]-psi(i-1,j)) );
	        end
	        % y fluxes
	        if (j==1)
	           Kp = 0.5*(K[i][j] + K(i,j+1));
	           Apsi[i][j] = Apsi[i][j] - dt/dy^2*( Kp*(psi(i,j+1)-psi[i][j]) - 0 );
	        elseif (j==JMAX)
	           Km = 0.5*(K[i][j]+K(i,j-1));
	           Apsi[i][j] = Apsi[i][j] - dt/dy^2*( 0 - Km*(psi[i][j]-psi(i,j-1)) );
	        else
	           Kp = 0.5*(K[i][j]+K(i,j+1));
	           Km = 0.5*(K[i][j]+K(i,j-1));
	           Apsi[i][j] = Apsi[i][j] - dt/dy^2*( Kp*(psi(i,j+1)-psi[i][j])...
	                                            -Km*(psi[i][j]-psi(i,j-1)) );
	        end
	    end
	end

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

	private static void printArray(double[] arr) {
		System.out.println(java.util.Arrays.toString(arr));
	}
}