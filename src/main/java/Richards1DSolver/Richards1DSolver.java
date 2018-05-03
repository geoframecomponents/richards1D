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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashMap;

import oms3.annotations.*;

import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import richards_classes.*;
import monodimensionalProblemTimeDependent.WriteNetCDFRichardsParameterization;

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
	
	// SOIL PARAMETERS
	@Description("The hydraulic conductivity at saturation")
	@In 
	@Unit ("m/s")
	public double[] ks;

	@Description("Saturated water content")
	@In 
	@Unit ("-")
	public double[] thetaS;

	@Description("Residual water content")
	@In 
	@Unit ("-")
	public double[] thetaR;
	
	@Description("First parameter of SWRC")
	@In 
	@Unit ("-")
	public double[] par1SWRC;
	
	@Description("Second parameter of SWRC")
	@In 
	@Unit ("-")
	public double[] par2SWRC;

	@Description("It is possibile to chose between 3 different models to compute "
			+ "the soil hydraulic properties: Van Genuchten; Brooks and Corey; Kosugi unimodal")
	@In 
	public String soilHydraulicModel;
	
	/////////////////////////////////////////////

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

	@Description("Time step of integration")
	@In
	@Unit ("s")
	public double timeDelta;

	@Description("Control variable for the integration time loop ")
	@In
	@Unit ("s")
	public double sumTimeDelta;

	@Description("Space step")
	@Unit ("m")
	double[] spaceDelta;

	@Description("Tolerance for Newton iteration")
	@In
	public double newtonTolerance;
	
	@Description("Control parameter for nested Newton algorithm:"
			+"0 --> simple Newton method"
			+"1 --> nested Newton method")
	@In
	public int nestedNewton; 

	@Description("Slope of the soil")
	@In 
	@Unit ("°")
	public double delta;

	// BOUNDARY CONDITIONS
	
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

	@Description("It is possibile to chose among 3 different kind "
			+ "of boundary condition at the bottom of the domain: "
			+ "- Dirichlet boundary condition --> Bottom Dirichlet"
			+ "- Neumann boundary condition --> Bottom Neumann"
			+ "- Impervious boundary condition --> Bottom Impervious")
	@In 
	public String bottomBCType;

	@Description("The first day of the simulation.")
	@In
	@Out
	public String inCurrentDate;

	@Description("Path of output files")
	@In
	public String dir;
	
	@Description("ArrayList of variable to be stored in the buffer writer")
	@Out
	public ArrayList<double[]> outputToBuffer;
	

	//////////////////////////////////////////
	//////////////////////////////////////////

	@Description("Maximun number of Newton iterations")
	final int MAXITER_NEWT = 50;

	@Description("Top boundary condition according with topBCType")
	@Unit ("")
	double topBC;

	@Description("Bottom boundary condition according with bottomBCType")
	@Unit ("")
	double bottomBC;

	@Description("Top boundary condition according with topBCType")
	@Unit ("")
	double tBC;

	@Description("Bottom boundary condition according with topBCType")
	@Unit ("")
	double bBC;

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

	@Description("Vector collects the dimensional water content of each cell at time level n")
	@Unit ("-")
	double[] volumes;
	
	@Description("Vector collects the dimensional water content of each cell at time level n+1")
	@Unit ("-")
	double[] volumesNew;
	
	@Description("Vector collects the adimensional water content of each soil cell and the water depth at soil surface at time level n+1")
	@Unit ("-")
	double[] thetasNew;
	
	@Description("Total volume at time level n")
	@Unit ("-")
	double volume;
	
	@Description("Total volume at time level n+1")
	@Unit ("-")
	double volumeNew;
	
	@Description("Volume error between time levels n+1 and n")
	@Unit ("-")
	double errorVolume;

	@Description("Vector collects velocities at cells' interfaces")
	@Unit ("m/s")
	double[] velocities;

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
	
	@Description("Vector containing the length of each control volume")
	double[] dx;
	
	double time=0;
	
	@Description("Variable containing theta-psi, hydraulic conductivity-psi, and moisture capacity-psi")
	LinkedHashMap<String,double[]> hydraulicParametrization;

	NestedNewton nestedNewtonAlg;

	

	SoilParametrization soilPar;
	TotalDepth totalDepth;

	BoundaryCondition topBoundaryCondition;
	BoundaryCondition bottomBoundaryCondition;

	WriteNetCDFRichardsParameterization writeSoilPar;

	
    ///////////////////////////////
	@Description("Initial condition for water head read from grid NetCDF file")
	@In
	@Unit("m")
	public double[] psiIC;

	@Description("z coordinate read from grid NetCDF file")
	@In
	@Unit("m")
	public double[] z;
	
	@Description("Space delta to compute gradients read from grid NetCDF file")
	@In 
	@Unit("m")
	public double[] spaceDeltaZ;
	
	@Description("Length of control volumes read from grid NetCDF file")
	@In 
	@Unit("m")
	public double[] deltaZ;
	
	
	@Execute
	public void solve() {

		System.out.println("RICHARDS 1D "+inCurrentDate);

		if(step==0){


			
			NUM_CONTROL_VOLUMES = z.length;

			psis 		  = new double[NUM_CONTROL_VOLUMES];
			k_t			  = 0.0;								
			k_b			  = 0.0;
			kappas 		  = new double[NUM_CONTROL_VOLUMES];
			volumes		  = new double[NUM_CONTROL_VOLUMES];
			volumesNew    = new double[NUM_CONTROL_VOLUMES];
			thetasNew     = new double[NUM_CONTROL_VOLUMES];
			velocities    = new double[NUM_CONTROL_VOLUMES+1];
			lowerDiagonal = new double[NUM_CONTROL_VOLUMES];
			mainDiagonal  = new double[NUM_CONTROL_VOLUMES];
			upperDiagonal = new double[NUM_CONTROL_VOLUMES];
			rhss 		  = new double[NUM_CONTROL_VOLUMES];
			kP 		      = 0.0;
			kM	          = 0.0;
			zeta 		  = new double[NUM_CONTROL_VOLUMES];
			spaceDelta    = new double[NUM_CONTROL_VOLUMES];
			dx			  = new double[NUM_CONTROL_VOLUMES];
			outputToBuffer= new ArrayList<double[]>();
			
			
			writeSoilPar = new WriteNetCDFRichardsParameterization();
			
			SimpleSoilParametrizationFactory soilParFactory = new SimpleSoilParametrizationFactory();
			//soilPar = soilParFactory.createSoilParametrization(soilHydraulicModel,alpha,n,psiE,lambda,rMedian,sigma,thetaR,thetaS,ks);
			soilPar = soilParFactory.createSoilParametrization(soilHydraulicModel);
			totalDepth = new TotalDepth();
			

			SimpleBoundaryConditionFactory boundCondFactory = new SimpleBoundaryConditionFactory();
			topBoundaryCondition = boundCondFactory.createBoundaryCondition(topBCType);		
			bottomBoundaryCondition = boundCondFactory.createBoundaryCondition(bottomBCType);	

			for(int i = 0; i < NUM_CONTROL_VOLUMES; i++) {
				psis[i] = psiIC[i];
				zeta[i] = z[i];
				spaceDelta[i] = spaceDeltaZ[i];
			}
			for(int i = 0; i < NUM_CONTROL_VOLUMES-1; i++) {
				dx[i] = deltaZ[i];
			}

			nestedNewtonAlg = new NestedNewton(nestedNewton, newtonTolerance, MAXITER_NEWT, NUM_CONTROL_VOLUMES, dx, soilPar, totalDepth);
			
			// conversion from degree to radiant of slope angle
			delta = delta*Math.PI/180;

			// Create and print a matrxi with data necessary to plot SWRC, hydraulic conductivity and moisture capacity parametrization
			//hydraulicParametrization = soilPar.hydraulicModelCurves1();
			//writeSoilPar.writeNetCDF(hydraulicParametrization, dir+"/HydraulicParameterization", soilHydraulicModel);
			
		} // chiudi step==0


		//time = time + tTimestep;

		topBC = 0.0;
		topBC = (inTopBC.get(1)[0]/1000)/tTimestep;
		
		bottomBC = 0.0;
		if(inBottomBC != null)
			bottomBC = inBottomBC.get(1)[0];
		
		outputToBuffer.clear();

		// Hydraulic conductivity are computed at time level n
		//k_b = soilPar.hydraulicConductivity(bottomBC);
		//for(int i = 0; i < NUM_CONTROL_VOLUMES; i++) {  
		//	kappas[i] = soilPar.hydraulicConductivity(psis[i]);
		//}
		
		sumTimeDelta = 0;
		while(sumTimeDelta < tTimestep) {

			if(sumTimeDelta + timeDelta>tTimestep) {
				timeDelta = tTimestep - sumTimeDelta;
			}
			sumTimeDelta = sumTimeDelta + timeDelta;

			// Compute volumes and hydraulic conductivity
			
			for(int i = 0; i < NUM_CONTROL_VOLUMES; i++) {  
				if(i==1) {
					soilPar.set(par1SWRC[i],par2SWRC[i],thetaR[i],thetaS[i],ks[i]);
					volumes[i] = soilPar.waterContent(psis[i])*dx[i];
					kappas[i] = soilPar.hydraulicConductivity(psis[i]);
					k_b = soilPar.hydraulicConductivity(bottomBC);  // I use the same parameters of the bottom cell
				} else if(i==NUM_CONTROL_VOLUMES-1) {
					volumes[i] = totalDepth.totalDepth(psis[i]);
					kappas[i] = soilPar.hydraulicConductivity(psis[i]);
				} else {
				soilPar.set(par1SWRC[i],par2SWRC[i],thetaR[i],thetaS[i],ks[i]);
				volumes[i] = soilPar.waterContent(psis[i])*dx[i];
				kappas[i] = soilPar.hydraulicConductivity(psis[i]);
				}
			}
			//k_b = soilPar.hydraulicConductivity(bottomBC);
			/* COEFFICIENT MATRIX IS BUILD BY THREE VECTORS COLLECTING ELEMENTS OF THE THREE DIAGONAL:
				   	 a lower diagonal psi_(i+1)
				   	 b main diagonal  psi_i
				   	 c upper diagonal psi_(i-1)
				   	 RIGHT HAND SIDE */
			for(int i = 0; i < NUM_CONTROL_VOLUMES; i++) {
				if( i == 0 ) {

					kP = 0.5*(kappas[i] + kappas[i+1]);
					if(bottomBCType.equalsIgnoreCase("Bottom Free Drainage") || bottomBCType.equalsIgnoreCase("BottomFreeDrainage")){
						kM = kappas[i];
					} else {
						kM = 0.5*(kappas[i] + k_b);
					}
					lowerDiagonal[i] =  bottomBoundaryCondition.lowerDiagonal(-999, kP, kM, spaceDelta[i+1], spaceDelta[i], timeDelta, delta);
					mainDiagonal[i] = bottomBoundaryCondition.mainDiagonal(-999, kP, kM, spaceDelta[i+1], spaceDelta[i], timeDelta, delta);
					upperDiagonal[i] = bottomBoundaryCondition.upperDiagonal(-999, kP, kM, spaceDelta[i+1], spaceDelta[i], timeDelta, delta);
					rhss[i] = volumes[i] + bottomBoundaryCondition.rightHandSide(bottomBC, kP, kM, spaceDelta[i+1], spaceDelta[i], timeDelta, delta);// + timeDelta*et[i]*(volumes[i] - thetaR[i]*dx[i]);;

				} else if(i == NUM_CONTROL_VOLUMES -1) {
					kP = 0;	
						// Water flux has to assigned as the minimum between rainfall rate and the maximum infiltrability of the soil
						//tBC = Math.min(topBC, kP);
					kM = 0.5*(kappas[i] + kappas[i-1]);
					lowerDiagonal[i] = topBoundaryCondition.lowerDiagonal(-999, kP, kM, spaceDelta[i], spaceDelta[i-1], timeDelta, delta); 
					mainDiagonal[i] = topBoundaryCondition.mainDiagonal(-999, kP, kM, spaceDelta[i], spaceDelta[i-1], timeDelta, delta);
					upperDiagonal[i] = topBoundaryCondition.upperDiagonal(-999, kP, kM,  spaceDelta[i], spaceDelta[i-1], timeDelta, delta);
					rhss[i] = volumes[i] + topBoundaryCondition.rightHandSide(topBC, kP, kM,  spaceDelta[i], spaceDelta[i-1], timeDelta, delta);// + timeDelta*et[i]*(volumes[i]);

				} else {

					kP = 0.5*(kappas[i] + kappas[i+1]);
					kM = 0.5*(kappas[i] + kappas[i-1]);
					lowerDiagonal[i] = -kM*timeDelta/spaceDelta[i]*1/Math.pow(Math.cos(delta),2);
					mainDiagonal[i] = kM*timeDelta/spaceDelta[i]*1/Math.pow(Math.cos(delta),2)  + kP*timeDelta/spaceDelta[i+1]*1/Math.pow(Math.cos(delta),2);
					upperDiagonal[i] = -kP*timeDelta/spaceDelta[i+1]*1/Math.pow(Math.cos(delta),2);
					rhss[i] = volumes[i] + timeDelta*(kP - kM);// + timeDelta*et[i]*(volumes[i] - thetaR[i]*dx[i]); 

				}
				
			}

			//// NESTED NEWTON ALGORITHM ////
			//nestedNewtonAlg.set(psis, mainDiagonal, upperDiagonal, lowerDiagonal, rhss);
			nestedNewtonAlg.set(psis, mainDiagonal, upperDiagonal, lowerDiagonal, rhss, thetaS, thetaR, par1SWRC, par2SWRC);
			psis = nestedNewtonAlg.solver();

			/* COMPUTE velocities AT CELL INTERFACES at time level n+1
			 * with hydraulic conductivity at time level n 
			 */ 
			volume = 0.0;
			volumeNew =0.0;
			for(int i = 0; i < NUM_CONTROL_VOLUMES; i++) {
				if( i == 0 ) {

					kP = 0.5*(kappas[i] + kappas[i+1]);
					if(bottomBCType.equalsIgnoreCase("Bottom Free Drainage") || bottomBCType.equalsIgnoreCase("BottomFreeDrainage")){
						kM = kappas[i];
						velocities[i] =  -kM;
					} else if(bottomBCType.equalsIgnoreCase("Bottom Impervious") || bottomBCType.equalsIgnoreCase("BottomImpervious")) {
						velocities[i] = + 0;
					}
					else {
						kM = 0.5*(kappas[i] + k_b);
						velocities[i] =  -kM * (psis[i]-bottomBC)/spaceDelta[i] - kM;
						
						soilPar.set(par1SWRC[i],par2SWRC[i],thetaR[i],thetaS[i],ks[i]);
						volumesNew[i] = soilPar.waterContent(psis[i])*dx[i];
						thetasNew[i] = soilPar.waterContent(psis[i]);
					}

				} else if(i == NUM_CONTROL_VOLUMES-1) {
					kP = kappas[i];
					velocities[i] =  -kP * (psis[i]-psis[i-1])/spaceDelta[i] - kP;
					
					volumesNew[i] = totalDepth.totalDepth(psis[i]);
					thetasNew[i] = totalDepth.totalDepth(psis[i]);
				} else {

					kP = 0.5*(kappas[i] + kappas[i+1]);
					velocities[i+1] =  -kP * (psis[i+1]-psis[i])/spaceDelta[i+1] - kP;
					
					soilPar.set(par1SWRC[i],par2SWRC[i],thetaR[i],thetaS[i],ks[i]);
					volumesNew[i] = soilPar.waterContent(psis[i])*dx[i];
					thetasNew[i] = soilPar.waterContent(psis[i]);
				}
				volume +=volumes[i];
				volumeNew +=volumesNew[i];
			}
			errorVolume = volumeNew - volume - timeDelta*(topBC - velocities[0]);
			//System.out.println("    errorMass: "+errorVolume);
			

		}
		outputToBuffer.add(psis);
		outputToBuffer.add(thetasNew);
		outputToBuffer.add(psiIC);
		outputToBuffer.add(velocities);
		outputToBuffer.add(new double[] {errorVolume});
		outputToBuffer.add(new double[] {topBC});
		outputToBuffer.add(new double[] {bottomBC});
		

		
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



