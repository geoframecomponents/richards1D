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
		int 	static 	days				= 24*3600;
		double 			ks 					= 0.062/day;  //[meter/second]
		double 			theta_s				= 0.41;         //[-] saturated water content
		double 			theta_r				= 0.095;        //[-] residuel water content
		double 			n					= 1.31;         // For Van Genuchten
		double 			m					= 1-1/n;        // For Van Genuchten
		double 	static	alpha				= 1.9;          // For Van Genuchten
		double 			psic 				= Math.pow(-1/alpha * (n-1)/n , 1/n);  // Where \frac{\partial\psi(\theta)}{\theta}=0
		double 			psi_r				= 0;					// Right boundary condition for pressure	
		double 			psi_l				= 0;					// Left boundary condition for pressure

		// Domain
		double 	static 	space_bottom		= 0;
		double 	static 	space_top			= 2;
		double 	static 	num_control_volumes	= 100; 
		double 	static 	space_delta			= (space_top - space_bottom) / num_control_volumes; 			// delta
		double 	static 	space_cv_centres[]	= seq(space_bottom + space_delta / 2,space_top - space_delta / 2,num_control_volumes); // Centres of the "control volumes"

		// Cycle variables
		int    	static 	max_iterations 		= 100000;
		double 	static 	time_end 			= Math.exp(3,5);            
		double 	static 	time_initial 		= 0.0;
		int 	static 	time_delta 			= 1000;
		int 			time 				= 0;

		// Working variables
		double 	static	gridvar				= time_delta / space_delta;
		double 	static	gridvarsq			= Math.pow((time_delta / space_delta),2);		
		double[] 		psis 				= new double[num_control_volumes];
		double 			k_r					= 0;
		double 			k_l					= 0;
		double[] 		kappas				= new double[num_control_volumes];
		double[] 		a 					= new double[num_control_volumes];
		double[] 		b 					= new double[num_control_volumes];
		double[] 		c 					= new double[num_control_volumes];
		double[] 		rhss				= new double[num_control_volumes];
		double 			k_p					= 0;
		double 			k_m					= 0;

		// Initial domain conditions
		for(int i = 0; i <= num_control_volumes; i++) {
			psi[i] = -y[i];
		}

		////////////////////
		//// Main cycle ////
		////////////////////
		for(int niter = 0; niter < max_iterations; niter++) {
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

		    // Gets every kappas right
		    k_r = kappa(psi_r); 
		    k_l = kappa(psi_l);
		    for(int i = 0; i < max_iterations; i++) {
			    theta[i] = thetaf(theta_s, theta_r, alpha, psi, n, m, psis[i]);
			    kappas[i] = kappa(alpha, theta_2, theta_r, psi[i]);           
			}

			// Right-hand side of equation
			for(int i = 0; i < max_iterations; i++) {
				if( i==1 ) {

		            k_p = 0.5*(kappas[i] + kappas[i+1]);
		            k_m = 0.5*(kappas[i] + k_l);
		            a[i] = 0;
		            b[i] = gridvarsq * (2 * k_m + k_p);
		            c[i] = -k_p * gridvarsq;
		            rhs[i] = thetas[i] + gridvar * (k_p-k_m) + 2 * k_m * gridvarsq * psi_l; 

				} else if(i == max_iterations) {

		            k_p = 0.5*(kappas[i] + k_r);
		            k_m = 0.5*(kappas[i] + k_l);
		            a[i] = -k_m * time_delta / space_;
		            b[i] = gridvarsq * (k_m + 2*k_p);
		            c[i] = 0;
		            rhs[i] = thetas[i] + gridvar * (k_p-k_m) + 2 * k_p * gridvarsq * psi_r; 

		            Kp = 0.5*(K(i)+KR);
		            Km = 0.5*(K(i)+K(i-1));
		            a(i) = -Km*dt/dx^2;
		            b(i) = dt/dx^2*(Km+2*Kp);
		            c(i) = 0;
		            rhs(i) = theta(i) + dt/dx*(Kp-Km) + 2*Kp*dt/dx^2*psi_r;              

				} else {

		            Kp = 0.5*(K(i)+K(i+1));
		            Km = 0.5*(K(i)+K(i-1));
		            a(i) = -Km*dt/dx^2;
		            b(i) = dt/dx^2*(Km+Kp);
		            c(i) = -Kp*dt/dx^2;
		            rhs(i) = theta(i) + dt/dx*(Kp-Km);           
				}
			}
		
			// "The show must go on"
		   	time =+ time_delta;

		   	// "Your time has come"
		   	if(time > time_end) {
		   		break;
		   	}
		}




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
	 * Returns a K given the parameters
	 * <p>
	 * 
	 *
	 *
	 *
	 * @param  alpha  	a real number, the value of the lower bound of the sequence
	 * @param  theta_s  an integer, the number of equally spaced points in the sequence 
	 * @param  theta_r  an integer, the number of equally spaced points in the sequence 
	 * @return 			a real number storing the control volume's K
	 */	
	public static double k(double alpha, double theta_s, double theta_r, double psi) {

		double kappa;
		double saturation;

		saturation = (thetaf(psi) - theta_r) / (theta_s - theta_r); 
		kappa = Ks * sqrt(saturation) * Math.pow(1 - (1 - Math.pow(Math.pow(sat,1/m)),m), 2);
		return kappa;
	}

	// Computes the \theta_{f}
	public static double thetaf(double theta_s, double theta_r, double alpha, double psi, double n, double m) {

		double theta_f;

		if(psi <= 0) {
		    theta_f = theta_r + (theta_s-theta_r)/( (1+Math.abs(alpha*psi)^n)^m );
		} else {
		    theta_f = theta_s;
		}
		return theta_f;
	}

}