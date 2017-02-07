
import org.junit.Test;

import richards.Richards1d;
import richards_classes.*;
public class Richards1dTest {

	private int 	days				= 24*3600;
	private double 	ks 					= 0.062/days;  	// [meter/second]
	private double 	theta_s				= 0.41;         // Saturated water content
	private double 	theta_r				= 0.095;        // Residual water content
	private double 	n					= 1.31;         // For Van Genuchten double 
	private double 	alpha				= 1.9;          // For Van Genuchten

	// Space
	private double 	space_bottom		= 0.0;
	private double 	space_top			= 2.0;
	private int 	NUM_CONTROL_VOLUMES	= 10;
	private double 	space_delta			= (space_top - space_bottom) / NUM_CONTROL_VOLUMES; 			// delta
	private double[] space_cv_centres	= DomainDiscretization.seq(space_bottom + space_delta / 2,space_top - space_delta / 2,NUM_CONTROL_VOLUMES); // Centres of the "control volumes"

	// Time
	private double 	time_end 			= 10000;
	private double 	time_initial 		= 0.0;
	private double 	time_delta 			= 1000.0;

	// Time and space
	private double 	gridvar				= time_delta / space_delta;
	private double 	gridvarsq			= time_delta / Math.pow(space_delta,2);

	// Cycle variables
	private int 	MAXITER 			= 1000;
	private int 	MAXITER_NEWT 		= 100000;
	private double 	newton_tolerance	= Math.pow(10,-12);
	
	@Test
	public void testTest(){ 
		
		Richards1d richards1d = new Richards1d(days, ks, theta_s, theta_r, n, alpha, space_bottom, space_top, NUM_CONTROL_VOLUMES, space_delta, space_cv_centres, time_end, time_initial, time_delta, gridvar, gridvarsq, MAXITER, MAXITER_NEWT, newton_tolerance);
		
		richards1d.solve();
	}
}
