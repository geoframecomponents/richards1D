import java.lang.Math

public class Richards1d {
	public static void main(String[] args) {
		// Model parameters - SI UNITS
		int static days	= 24*3600;
		double ks 		= 0.062/day;  //[meter/second]
		double theta_s 	= 0.41;         //[-] saturated water content
		double theta_r 	= 0.095;        //[-] residuel water content
		double n       	= 1.31;         // For Van Genuchten
		double m       	= 1-1/n;        // For Van Genuchten
		double alpha   	= 1.9;          // For Van Genuchten
		double psic    	= Math.pow(-1/alpha * (n-1)/n , 1/n);  // Where \frac{\partial\psi(\theta)}{\theta}=0
		double psi_r	= 0;					// Right boundary condition for pressure	
		double psi_l	= 0;					// Left boundary condition for pressure
		// Domain
		double static y_bottom = 0;
		double static y_top = 2;
		double static num_control_volumes = 100; 
		double static dy = (y_top-y_bottom)/num_control_volumes; 			// delta
		double static y[] = seq(y_bottom+dy/2,y_top-dy/2,num_control_volumes); // Centres of the "control volumes"

		// Cycle variables
		int static max_iterations 	= 100000;
		double static time_end 		= Math.exp(3,5);            
		double static time_initial 	= 0.0;
		int static time_delta 		= 1000;
		int time 					= 0;

		// Working variables
		double[] psi 				= new double[num_control_volumes];
		double kappa_r				= 0;
		double kappa_l				= 0;

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

		    if(time <= Math.pow(10,5)) {
		        psi_r = -0.05 + 0.03*Math.sin(2*Math.PI * time/Math.pow(10,5));
		    } else if(time>1e5 && time<=1.8e5) {
		        psi_r = +0.1;
		    }
		    else {
		        psi_r = -0.05+2952.45*Math.exp(-time/18204.8);
		    }

		    kappa_r = kappa(psi_r); 
		    kappa_l = kappa(psi_l);
		    for(int i = 0; i < max_iterations; i++) {
			    theta(i)=Thetaf(psi(i));
			    K(i) = kappa(psi(i));           
			}

		
		   	time = time + time_delta;
		   	if(time > time_end) {
		   		
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

	// Computes the Kappa
	public static double k(double alpha, double theta_s, double theta_r, double psis) {
		double kappa;
		double saturation;

		saturation = 
	}

	// Computes the \theta_{f}
	public static double thetaf(double theta_s, double theta_r, double alpha, double psi, double n, double m) {
		double theta_f;

		if(psi <= 0) {
		    theta_f = theta_r + (theta_s-theta_r)/( (1+Math.abs(alpha*psi)^n)^m );
		} else {
		    theta_f = theta_s;
		}
		return 
	}

function K=kappa(psi)
global alpha thetas thetar n m Ks psic
sat = (Thetaf(psi)-thetar)/(thetas-thetar); %saturation is between 0 1 
K = Ks*sqrt(sat)*(1-(1-sat^(1/m))^m)^2;

}