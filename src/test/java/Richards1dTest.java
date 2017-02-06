
import org.junit.Test;

import richards.Richards1d;
import richards_classes.*;
public class Richards1dTest {
	static  int 	days				= 24*3600;
	static	double 	ks 					= 0.062/days;  	// [meter/second]
	static	double 	theta_s				= 0.41;         // Saturated water content
	static	double 	theta_r				= 0.095;        // Residual water content
	static	double 	n					= 1.31;         // For Van Genuchten
	static	double 	m					= 1-1/n;        // For Van Genuchten
	static	double 	alpha				= 1.9;          // For Van Genuchten

	// Space
	static 	double 	space_bottom		= 0.0;
	static 	double 	space_top			= 2.0;
	static 	int 	NUM_CONTROL_VOLUMES	= 10; 
	static 	double 	space_delta			= (space_top - space_bottom) / NUM_CONTROL_VOLUMES; 			// delta
	static 	double[] space_cv_centres	= DomainDiscretization.seq(space_bottom + space_delta / 2,space_top - space_delta / 2,NUM_CONTROL_VOLUMES); // Centres of the "control volumes"

	// Time
	static 	double 	time_end 			= 10000;//3*Math.pow(10,8);            
	static 	double 	time_initial 		= 0.0;
	static 	double 	time_delta 			= 1000.0;

	// Time and space
	static	double 	gridvar				= time_delta / space_delta;
	static	double 	gridvarsq			= time_delta / Math.pow(space_delta,2);		

	// Cycle variables
	static 	int 	MAXITER 			= 1000;
	static  int 	MAXITER_NEWT 		= 100000;
	static	double 	newton_tolerance	= Math.pow(10,-12);
	
	@Test
	public void testTest(){ 
		
		Richards1d richards1d = new Richards1d(days, ks, theta_s, theta_r, n, alpha, space_bottom, space_top, NUM_CONTROL_VOLUMES, space_delta, space_cv_centres, time_end, time_initial, time_delta, gridvar, gridvarsq, MAXITER, MAXITER_NEWT, newton_tolerance);
		
		richards1d.solve();
	}
}
