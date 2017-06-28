/*
 * GNU GPL v3 License
 *
 * Copyright 2016 Marialaura Bancheri
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

package Richards1DSolver;

import java.util.HashMap;

import oms3.annotations.*;

import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import richards_classes.*;

@Description("Solve the Richards equation for the 1D domain.")
@Documentation("")
@Author(name = "Aaron Iemma, Niccolo' Tubini, Francesco Serafin, Michael Dumbser and Riccardo Rigon", contact = "tubini.niccolo@gmail.com")
@Keywords("Hydrology, Richards, Infiltration")
@Bibliography("Casulli (2010)")
//@Label(JGTConstants.HYDROGEOMORPHOLOGY)
//@Name("shortradbal")
//@Status(Status.CERTIFIED)
@License("General Public License Version 3 (GPLv3)")
public class Richards1DSolver {
	@Description("The hydraulic conductivity at saturation")
	@In 
	@Unit ("m/s")
	public double ks;

	@Description("Saturated water content")
	@In 
	@Unit ("-")
	public double thetaS;

	@Description("Residual water content")
	@In 
	@Unit ("-")
	public double thetaR;

	@Description("Exponent of Van Genuchten model")
	@In
	@Unit ("-")
	public double n;            

	@Description("Parameter of Van Genuchten model")
	@In 
	@Unit ("m^(-1)")
	public double alpha;

	@Description("Exponent of Brooks and Corey model")
	@In
	@Unit ("-")
	public double lambda;

	@Description("Parameter of Brooks and Corey model")
	@In
	@Unit ("m")
	public double psiE;

	@Description("Median pore radius of the pore size-distribution for Kosugi Model")
	@In
	@Unit ("m")
	public double rMedian;

	@Description("Standard deviation of the pore size-distribution for Kosugi Model")
	@In
	@Unit ("m")
	public double sigma;

	@Description("It is possibile to chose between 3 different models to compute "
			+ "the soil hydraulic properties: Van Genuchten; Brooks and Corey; Kosugi unimodal")
	@In 
	public String soilHydraulicModel;

	@Description("Depth of the soil column")
	@In
	@Unit ("m")
	public double spaceBottom;

	@Description("")
	@Unit ("m")
	public double spaceTop=0;

	@Description("Number of control volume for domain discetrization")
	//@In
	@Unit (" ")
	int NUM_CONTROL_VOLUMES; 

	@Description("It is needed to iterate on the date")
	int step;

	@Description("Time amount at every time-loop")
	@In
	@Unit ("s")
	public double tTimestep;

	@Description("Space step")
	@Unit ("m")
	double[] spaceDelta;

	@Description("Tolerance for Newton iteration")
	@In
	public double newtonTolerance;

	@Description("Initial condition for the soil suction")
	@In 
	@Unit ("m")
	public double[] iC;

	@Description("Slope of the soil")
	@In 
	@Unit ("°")
	public double delta;
	
	@Description("Depth at which the initial condition is defined")
	@In 
	@Unit ("m")
	public double[] depth;

	@Description("The HashMap with the time series of the boundary condition at the top of soil column")
	@In
	@Unit ("m")
	public HashMap<Integer, double[]> inTopBC;

	@Description("It is possibile to chose between 2 different kind "
			+ "of boundary condition at the top of the domain: "
			+ "- Dirichlet boundary condition --> Top Dirichlet"
			+ "- Neumann boundary condition --> Top Neumann")
	@In 
	public String topBCType;

	@Description("The HashMap with the time series of the boundary condition at the bottom of soil column")
	@In
	@Unit ("m")
	public HashMap<Integer, double[]> inBottomBC;

	@Description("It is possibile to chose between 2 different kind "
			+ "of boundary condition at the bottom of the domain: "
			+ "- Dirichlet boundary condition --> Bottom Dirichlet"
			+ "- Neumann boundary condition --> Bottom Neumann")
	@In 
	public String bottomBCType;

	@Description("The first day of the simulation.")
	@In
	public String inCurrentDate;

	@Description("Path of output folder")
	@In
	public String dir;

	@Description("Control parameter for nested Newton algorithm:"
			+"0 --> simple Newton method"
			+"1 --> nested Newton method")
	@In
	public int nestedNewton; 

	@Description("Maximun number of Newton iterations")
	final int MAXITER_NEWT = 50;

	@Description("Top boundary condition according with topBCType")
	@Unit ("")
	double topBC;

	@Description("Bottom boundary condition according with bottomBCType")
	@Unit ("")
	double bottomBC;

	@Description("Psi values")
	@Unit ("m")
	double[] psis;

	@Description("Hydraulic conductivity at the top of the soil column")
	@Unit ("m/s")
	double k_t;

	@Description("Hydraulic conductivity at the bottom of the soil column")
	@Unit ("m/s")
	double k_b;

	@Description("Vector collects the hydraulic conductivity of each cell")
	@Unit ("m/s")
	double[] kappas;

	@Description("Vector collects the adimensional water content of each cell")
	@Unit ("-")
	double[] thetas;

	@Description("Vector collects the lower diagonal entries of the coefficient matrix")
	@Unit ("?")
	double[] lowerDiagonal;

	@Description("Vector collects the main diagonal entries of the coefficient matrix")
	@Unit ("?")
	double[] mainDiagonal;

	@Description("Vector collects the upper diagonal entries of the coefficient matrix")
	@Unit ("?")
	double[] upperDiagonal;

	@Description("Right hand side vector of the scalar equation to solve")
	@Unit ("-")
	double[] rhss;

	@Description("Hydraulic conductivity at the cell interface i+1/2")
	@Unit ("m/s")
	double kP;

	@Description("Hydraulic conductivity at the cell interface i-1/2")
	@Unit ("m/s")
	double kM;

	@Description("Risidual of the outer iteration of Nested Newton method")
	double outerResidual;

	@Description("Risidual of the inner iteration of Nested Newton method")
	double innerResidual;

	@Description("Vector containing the z coordinates of the centres of control volumes")
	double[] zeta;

	double time=0;

	NestedNewton nestedNewtonAlg;
	
	PrintTXT print;

	SoilParametrization soilPar;

	BoundaryCondition topBoundaryCondition;
	BoundaryCondition bottomBoundaryCondition;

	@Execute
	public void solve() {

		System.out.println("RICHARDS 1D "+inCurrentDate);

		if(step==0){

			iC = iC.getClass().cast(checkIC(depth, iC, depth));

			NUM_CONTROL_VOLUMES = iC.length;

			psis 		  = new double[NUM_CONTROL_VOLUMES];		
			k_t			  = 0.0;								
			k_b			  = 0.0;
			kappas 		  = new double[NUM_CONTROL_VOLUMES];
			thetas 		  = new double[NUM_CONTROL_VOLUMES];
			lowerDiagonal = new double[NUM_CONTROL_VOLUMES];
			mainDiagonal  = new double[NUM_CONTROL_VOLUMES];
			upperDiagonal = new double[NUM_CONTROL_VOLUMES];
			rhss 		  = new double[NUM_CONTROL_VOLUMES];
			kP 		      = 0.0;
			kM	          = 0.0;
			zeta 		  = new double[NUM_CONTROL_VOLUMES];
			spaceDelta    = new double[NUM_CONTROL_VOLUMES+1];

			print = new PrintTXT();

			SimpleSoilParametrizationFactory soilParFactory = new SimpleSoilParametrizationFactory();
			soilPar = soilParFactory.createSoilParametrization(soilHydraulicModel,alpha,n,psiE,lambda,rMedian,sigma,thetaR,thetaS,ks);
			
			nestedNewtonAlg = new NestedNewton(nestedNewton, newtonTolerance, MAXITER_NEWT, NUM_CONTROL_VOLUMES, soilPar);

			SimpleBoundaryConditionFactory boundCondFactory = new SimpleBoundaryConditionFactory();
			topBoundaryCondition = boundCondFactory.createBoundaryCondition(topBCType);		
			bottomBoundaryCondition = boundCondFactory.createBoundaryCondition(bottomBCType);	
			// Initial domain conditions
			/* da rivedere: forse è meglio mettere z=0 alla superficie anzichè alla base della colonna di suolo
			*  oltre a rivedere la definizione della variabile spaceDelta è da rivedere anche il file della condizione iniziale
			*  in cui la profondità deve essere data come negativa
			*/
			for(int i = 0; i < NUM_CONTROL_VOLUMES; i++) {
				psis[i] = iC[i];
				zeta[i] = spaceBottom-depth[i];
			}
			for(int i = 0; i <=NUM_CONTROL_VOLUMES; i++) {
				if (i==0){
					spaceDelta[i] = zeta[i];
				}else if (i== NUM_CONTROL_VOLUMES){
					spaceDelta[i] = spaceBottom-zeta[i-1];
				}else{
					spaceDelta[i] = zeta[i]-zeta[i-1];
				}
			}
			
			// conversion from degree to radiant of slope angle
			delta = delta*Math.PI/180;

		} // chiudi step==0

		time = time + tTimestep;

		topBC = 0.0;
		if(inTopBC != null){
			if(topBCType.equalsIgnoreCase("Top Neumann") || bottomBCType.equalsIgnoreCase("TopNeumann")){
				topBC = (inTopBC.get(1)[0]/1000)/tTimestep;
			} else if(topBCType.equalsIgnoreCase("Top Dirichlet") || bottomBCType.equalsIgnoreCase("TopDirichlet")){
				topBC = inTopBC.get(1)[0]/1000;
			}
		}
		bottomBC = 0.0;
		if(inBottomBC != null)
			bottomBC = inBottomBC.get(1)[0];

		// Compute hydraulic conductivity and water content
		k_t = soilPar.hydraulicConductivity(topBC);
		k_b = soilPar.hydraulicConductivity(bottomBC);
		for(int i = 0; i < NUM_CONTROL_VOLUMES; i++) {           
			thetas[i] = soilPar.waterContent(psis[i]);
			kappas[i] = soilPar.hydraulicConductivity(psis[i]);
		}
		
		/* COEFFICIENT MATRIX IS BUILD BY THREE VECTORS COLLECTING ELEMENTS OF THE THREE DIAGONAL:
				   	 a lower diagonal
				   	 b main diagonal
				   	 c upper diagonal
				   	 RIGHT HAND SIDE */
		for(int i = 0; i < NUM_CONTROL_VOLUMES; i++) {
			if( i == 0 ) {

				kP = 0.5*(kappas[i] + kappas[i+1]);
				if(bottomBCType.equalsIgnoreCase("Bottom Free Drainage") || bottomBCType.equalsIgnoreCase("BottomFreeDrainage")){
					kM = kappas[i];
				} else {
					kM = 0.5*(kappas[i] + k_b);
				}
				lowerDiagonal[i] =  bottomBoundaryCondition.lowerDiagonal(-999, kP, kM, spaceDelta[i+1], spaceDelta[i], tTimestep, delta);
				mainDiagonal[i] = bottomBoundaryCondition.mainDiagonal(-999, kP, kM, spaceDelta[i+1], spaceDelta[i], tTimestep, delta);
				upperDiagonal[i] = bottomBoundaryCondition.upperDiagonal(-999, kP, kM, spaceDelta[i+1], spaceDelta[i], tTimestep, delta);
				rhss[i] = thetas[i] + bottomBoundaryCondition.rightHandSide(bottomBC, kP, kM, spaceDelta[i+1], spaceDelta[i], tTimestep, delta);

			} else if(i == NUM_CONTROL_VOLUMES -1) {

				if(bottomBCType.equalsIgnoreCase("Top Neumann") || bottomBCType.equalsIgnoreCase("TopNeumann")){
					kP = kappas[i];
				} else {
					kP = 0.5*(kappas[i] + k_t);
				}
				kM = 0.5*(kappas[i] + kappas[i-1]);
				lowerDiagonal[i] = topBoundaryCondition.lowerDiagonal(-999, kP, kM, spaceDelta[i+1], spaceDelta[i], tTimestep, delta); 
				mainDiagonal[i] = topBoundaryCondition.mainDiagonal(-999, kP, kM, spaceDelta[i+1], spaceDelta[i], tTimestep, delta);
				upperDiagonal[i] = topBoundaryCondition.upperDiagonal(-999, kP, kM,  spaceDelta[i+1], spaceDelta[i], tTimestep, delta);
				rhss[i] = thetas[i] + topBoundaryCondition.rightHandSide(topBC, kP, kM, spaceDelta[i+1], spaceDelta[i], tTimestep, delta); 

			} else {

				kP = 0.5*(kappas[i] + kappas[i+1]);
				kM = 0.5*(kappas[i] + kappas[i-1]);
				lowerDiagonal[i] = -kM * tTimestep/(spaceDelta[i]/2+spaceDelta[i+1]/2)*1/spaceDelta[i]*1/Math.pow(Math.cos(delta),2);
				mainDiagonal[i] = tTimestep/(spaceDelta[i]/2+spaceDelta[i+1]/2)*1/spaceDelta[i]*1/Math.pow(Math.cos(delta),2) * kM + tTimestep/(spaceDelta[i]/2+spaceDelta[i+1]/2)*1/spaceDelta[i+1]*1/Math.pow(Math.cos(delta),2)*kP;
				upperDiagonal[i] = -kP * tTimestep/(spaceDelta[i]/2+spaceDelta[i+1]/2)*1/spaceDelta[i+1]*1/Math.pow(Math.cos(delta),2);
				rhss[i] = thetas[i] + tTimestep/(spaceDelta[i]/2+spaceDelta[i+1]/2) * kP - tTimestep/(spaceDelta[i]/2+spaceDelta[i+1]/2)*kM; 

			}
		}
		
		//// NESTED NEWTON ALGORITHM ////
		nestedNewtonAlg.set(psis, mainDiagonal, upperDiagonal, lowerDiagonal, rhss);
		psis = nestedNewtonAlg.solver();

		//// PRINT OUTPUT FILES ////
		print.setValueFirstVector(depth);
		print.setValueSecondVector(psis);
		print.PrintTwoVectors(dir,"Psi_"+step+".csv", "Psi values at time: "+inCurrentDate, "Depth[m],Psi[m] ");

		for(int i = 0; i < NUM_CONTROL_VOLUMES; i++) {           
			thetas[i] = soilPar.waterContent(psis[i]);
		}
		print.setValueFirstVector(depth);
		print.setValueSecondVector(thetas);
		print.PrintTwoVectors(dir, "Theta_"+step+".csv", "Theta values at time: "+inCurrentDate, "Depth[m],Theta[-] ");

		step++;

	} //// MAIN CYCLE END ////

	private Object checkIC(double[] from_depth, Object iC, double[] to_depth) {
		if (iC instanceof Double || iC instanceof Float) {
			return iC;
		} else if (iC instanceof double[]){
		    return checkIClength(from_depth, (double []) iC, to_depth);
		} else {
			throw new UnsupportedOperationException();
		}
	}

	private Object checkIClength(double[] from_depth, double[] iC, double[] to_depth) {
		if (to_depth.length == iC.length) {
			return iC;
		} else if ((to_depth.length + 1) == iC.length) {
			return linearInterp(from_depth, iC, to_depth);
		} else {
			String msg = "Cannot process " + iC.length + " values of IC";
			msg += " for " + to_depth.length + " vertices";
			throw new IndexOutOfBoundsException(msg);
		}
	}

	private double[] linearInterp(double[] x, double[] y, double[] xi) {
       LinearInterpolator li = new LinearInterpolator(); // or other interpolator
       PolynomialSplineFunction psf = li.interpolate(x, y);

       double[] yi = new double[xi.length];
       for (int i = 0; i < xi.length; i++) {
           yi[i] = psf.value(xi[i]);
       }
       return yi;
    }
}  /// CLOSE Richards1d ///



