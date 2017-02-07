package richards;

import java.lang.Math;
import oms3.annotations.*;

import richards_classes.*;
import richards.Richards1d;

@Author(
		name="Niccolò Tubini, Francesco Serafin, Aaron Iemma",
		org="DICAM - Departement of Environmental and Civil Engineering - Trento, UNITN"
	)
@Description("A class to resolve the water pressures in a moisted soil for a pourous domain, unsing the Richards method")
@Keywords("Richards,1D")


public class Richards1d_OMS {
	
	@In public int days					;
	private double ks 					= 0.062/days;  	// [meter/second]
	@In public double theta_s			;         		// Saturated water content
	@In public double theta_r			;        		// Residual water content
	@In public double n					;         		// For Van Genuchten double 
	@In public double alpha				;          		// For Van Genuchten

	// Space
	@In public double space_bottom		;
	@In public double space_top			;
	@In public int NUM_CONTROL_VOLUMES	;
	private double space_delta			= (space_top - space_bottom) / NUM_CONTROL_VOLUMES; 			// delta
	private double[] space_cv_centres	= DomainDiscretization.seq(space_bottom + space_delta / 2,space_top - space_delta / 2,NUM_CONTROL_VOLUMES); // Centres of the "control volumes"

	// Time
	@In public double time_end 			;
	@In public double time_initial 		;
	@In public double time_delta 		;

	// Time and space
	private double 	gridvar				= time_delta / space_delta;
	private double 	gridvarsq			= time_delta / Math.pow(space_delta,2);

	// Cycle variables
	@In public int 	MAXITER 			;
	@In public int 	MAXITER_NEWT 		;
	private double 	newton_tolerance	= Math.pow(10,-12);
	
	
	@Execute
	public void doit(){ 	
		Richards1d richards1d = new Richards1d(days, ks, theta_s, theta_r, n, alpha, space_bottom, space_top, NUM_CONTROL_VOLUMES, space_delta, space_cv_centres, time_end, time_initial, time_delta, gridvar, gridvarsq, MAXITER, MAXITER_NEWT, newton_tolerance);
		richards1d.solve();
	}
}
