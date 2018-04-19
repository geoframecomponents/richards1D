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

import Energy1DSolver.*;
import richards_classes.*;
import thermalParameterization.SimpleThermalSoilParameterizationFactory;
import thermalParameterization.ThermalSoilParameterization;

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

	@Description("Time interval of boundary conditions")
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

	@Description("Initial condition for the soil suction")
	@In 
	@Unit ("m")
	public double[] iC;

	@Description("Source/sink value")
	@In 
	@Unit ("s^(-1)")
	public double[] sourceSink;
	
	@Description("Source/sink value at each time step of integration")
	@Unit ("s^(-1)")
	public double sS;

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
	
	@Description("Initial condition for the soil temperature")
	@In 
	@Unit ("C")
	public double[] temperatureIC;
	
	@Description("Initial condition for the soil temperature")
	@In 
	@Unit ("m")
	public double[] temperatureDepth;

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
	public String inCurrentDate;

	@Description("Path of output folder")
	@In
	public String dir;

	@Description("Control parameter for nested Newton algorithm:"
			+"0 --> simple Newton method"
			+"1 --> nested Newton method")
	@In
	public int nestedNewton; 
///////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////
	@Description("Type of energy integrator:"
			+" no energy --> do not solve energy equation"
			+" pure diffusion --> consider only diffusion fluxes"
			+" diffusion convection --> consider both diffusion and convection")
	@In
	public String energySolver; 
	
	@In
	public String tempTopBCType;
	@In
	public String tempBottomBCType;
	@In
	public String tempConvectionTopBCType;
	@In
	public String tempConvectionBottomBCType;
	
	@In
	public String soilThermalModel;
	
	@In
	public String soilType;
	
	@In
	public String kappaWithTemperatureModel;
	
	@In
	public double sandFraction;
	
	@In
	public double clayFraction;
	
	@In
	public double lambda0;
	
	@In
	public double quartzFraction;
	
	@Description("The HashMap with the time series of the boundary condition at the bottom of soil column")
	@In
	@Unit ("m")
	public HashMap<Integer, double[]> inTemperatureBottomBC;
	
	@Description("The HashMap with the time series of the boundary condition at the bottom of soil column")
	@In
	@Unit ("m")
	public HashMap<Integer, double[]> inTemperatureTopBC;
		
	double[] temperatures;
	 
	double temperatureBottomBC;
	
	double temperatureTopBC;
	
	EnergyIntegrator energyIntegrator;
	
	BoundaryCondition tempTopBoundaryCondition;
	BoundaryCondition tempBottomBoundaryCondition;
	BoundaryCondition tempConvectionTopBoundaryCondition;
	BoundaryCondition tempConvectionBottomBoundaryCondition;
	
	ThermalSoilParameterization thermalSoilPar;
	HydraulicConductivityTemperature kappaWithTemperature;
	
	SimpleThermalSoilParameterizationFactory thermalSoilParFactory;
	
	double[][] thermalParametrization;
	double[] peclets;

	
	///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
	@Description("Maximun number of Newton iterations")
	final int MAXITER_NEWT = 50;

	@Description("Top boundary condition according with topBCType")
	@Unit ("")
	double topBC;

	@Description("Bottom boundary condition according with bottomBCType")
	@Unit ("")
	double bottomBC;

	@Description("Top boundary condition for each time step of integration")
	@Unit ("")
	double tBC;

	@Description("Bottom boundary condition according with topBCType")
	@Unit ("")
	double bBC;

	@Description("Psi values")
	@Unit ("m")
	double[] psis;
	
	@Description("Psi values")
	@Unit ("m")
	double[] psisNew;

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
	
	@Description("Vector collects the dimensional water content of each cell at time level n")
	@Unit ("-")
	double[] thetas;
	
	@Description("Vector collects the dimensional water content of each cell at time level n+1")
	@Unit ("-")
	double[] volumesNew;
	
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

	NestedNewton nestedNewtonAlg;

	PrintTXT print;

	SoilParametrization soilPar;
	TotalDepth totalDepth;

	BoundaryCondition topBoundaryCondition;
	BoundaryCondition bottomBoundaryCondition;


	double[][] hydraulicParametrization;

	@Execute
	public void solve() {

		System.out.println("RICHARDS 1D "+inCurrentDate);

		if(step==0){

			iC = iC.getClass().cast(checkIC(depth, iC, depth));
			temperatureIC = temperatureIC.getClass().cast(checkIC(temperatureDepth, temperatureIC, temperatureDepth));
			sourceSink = sourceSink.getClass().cast(checkIC(depth, sourceSink, depth));

			NUM_CONTROL_VOLUMES = iC.length;

			psis 		  = new double[NUM_CONTROL_VOLUMES];
			psis 		  = new double[NUM_CONTROL_VOLUMES];
			temperatures  = new double[NUM_CONTROL_VOLUMES-1];
			k_t			  = 0.0;								
			k_b			  = 0.0;
			kappas 		  = new double[NUM_CONTROL_VOLUMES];
			volumes		  = new double[NUM_CONTROL_VOLUMES];
			thetas		  = new double[NUM_CONTROL_VOLUMES];
			volumesNew    = new double[NUM_CONTROL_VOLUMES];
			velocities    = new double[NUM_CONTROL_VOLUMES];
			peclets       = new double[NUM_CONTROL_VOLUMES];
			lowerDiagonal = new double[NUM_CONTROL_VOLUMES];
			mainDiagonal  = new double[NUM_CONTROL_VOLUMES];
			upperDiagonal = new double[NUM_CONTROL_VOLUMES];
			rhss 		  = new double[NUM_CONTROL_VOLUMES];
			kP 		      = 0.0;
			kM	          = 0.0;
			zeta 		  = new double[NUM_CONTROL_VOLUMES];
			spaceDelta    = new double[NUM_CONTROL_VOLUMES];
			dx			  = new double[NUM_CONTROL_VOLUMES];
			print = new PrintTXT();
			
			SimpleHydraulicConductivityTemperatureFactory hydraulicConductivityTemperatureFactory = new SimpleHydraulicConductivityTemperatureFactory();
			kappaWithTemperature = hydraulicConductivityTemperatureFactory.createHydraulicConductivityTemperature(kappaWithTemperatureModel);
			
			SimpleSoilParametrizationFactory soilParFactory = new SimpleSoilParametrizationFactory();
			soilPar = soilParFactory.createSoilParametrization(soilHydraulicModel,kappaWithTemperature,alpha,n,psiE,lambda,rMedian,sigma,thetaR,thetaS,ks);
			totalDepth = new TotalDepth();
			
			thermalSoilParFactory = new SimpleThermalSoilParameterizationFactory();
			thermalSoilPar = thermalSoilParFactory.createThermalSoilParameterization( soilThermalModel, thetaS, sandFraction, clayFraction, soilPar, soilType, lambda0, quartzFraction);
			
			
			SimpleBoundaryConditionFactory boundCondFactory = new SimpleBoundaryConditionFactory();
			topBoundaryCondition = boundCondFactory.createBoundaryCondition(topBCType);		
			bottomBoundaryCondition = boundCondFactory.createBoundaryCondition(bottomBCType);	
			tempTopBoundaryCondition = boundCondFactory.createBoundaryCondition(tempTopBCType);		
			tempBottomBoundaryCondition = boundCondFactory.createBoundaryCondition(tempBottomBCType);
			tempConvectionTopBoundaryCondition = boundCondFactory.createBoundaryCondition(tempConvectionTopBCType);		
			tempConvectionBottomBoundaryCondition = boundCondFactory.createBoundaryCondition(tempConvectionBottomBCType);
			
     		// Initial domain conditions
			/* da rivedere: forse è meglio mettere z=0 alla superficie anzichè alla base della colonna di suolo
			 *  oltre a rivedere la definizione della variabile spaceDelta è da rivedere anche il file della condizione iniziale
			 *  in cui la profondità deve essere data come negativa
			 */
			for(int i = 0; i < NUM_CONTROL_VOLUMES; i++) {
				psis[i] = iC[i];
				zeta[i] = spaceBottom-depth[i];
				//temperatures[i] = temperatureIC[i];
			}
			
			for(int i = 0; i < NUM_CONTROL_VOLUMES-1; i++) {
				temperatures[i] = temperatureIC[i];
			}
	
			for(int i = 0; i <NUM_CONTROL_VOLUMES; i++) {
				if (i==0){
					spaceDelta[i] = zeta[i];
				}//else if (i== NUM_CONTROL_VOLUMES-1){
					//spaceDelta[i] = zeta[i+1]-zeta[i];
				//}
				else{
					spaceDelta[i] = zeta[i]-zeta[i-1];
					//System.out.println("i:"+i+"  "+zeta[i]+"  "+zeta[i-1]);
				}
			}
			
			for(int i = 0; i <NUM_CONTROL_VOLUMES; i++) {
				if (i==0){
					dx[i] = spaceDelta[i]+spaceDelta[i+1]/2;
					//System.out.println("i:"+i+"  "+spaceDelta[i]+"  "+spaceDelta[i+1]);
				} else if (i== NUM_CONTROL_VOLUMES-1) {
					dx[i] =0;
				}
				else if (i== NUM_CONTROL_VOLUMES-2){
					dx[i] = spaceDelta[i]/2+spaceDelta[i+1];
					//System.out.println("i:"+i+"  "+spaceDelta[i]+"  "+spaceDelta[i+1]+"  "+dx[i]);
				}else{
					dx[i] = (spaceDelta[i]+spaceDelta[i+1])/2;
					//System.out.println("i:"+i+"  "+spaceDelta[i]+"  "+spaceDelta[i+1]+"  "+dx[i]);
				}
			}
			nestedNewtonAlg = new NestedNewton(nestedNewton, newtonTolerance, MAXITER_NEWT, NUM_CONTROL_VOLUMES, dx, soilPar, totalDepth);
			
			
			SimpleEnergyIntegratorFactory energyIntegratorFactory = new SimpleEnergyIntegratorFactory();
			energyIntegrator = energyIntegratorFactory.createEnergyIntegrator(energySolver,thermalSoilPar, NUM_CONTROL_VOLUMES, dx, spaceDelta, tempTopBoundaryCondition, tempBottomBoundaryCondition,
																					tempConvectionTopBoundaryCondition, tempConvectionBottomBoundaryCondition);
	
			
			// conversion from degree to radiant of slope angle
			delta = delta*Math.PI/180;

			// Create and print a matrix with data necessary to plot SWRC, hydraulic conductivity and moisture capacity parametrization
			hydraulicParametrization = soilPar.hydraulicModelCurves();
			print.setValueMatrix(hydraulicParametrization);
			print.PrintMatrix(dir, "Hydraulic Parametrization"+".csv", soilHydraulicModel, "Psi[m], Se[-], Theta[-], dTheta[1/m], K[m/s]");
			
			thermalParametrization = this.thermalSoilPar.thermalSoilPropertiesCurves();
			print.setValueMatrix(thermalParametrization);
			print.PrintMatrix(dir, "Thermal Parametrization"+".csv", soilThermalModel, "Psi[m], Heat Capacity [J K^-1 m^-3], Thermal Conductivity[W m^-1 K^-1], NaN, NaN");

		} // chiudi step==0


		//time = time + tTimestep;

		topBC = 0.0;
		//if(inTopBC != null)
		topBC = (inTopBC.get(1)[0]/1000)/tTimestep;
		
		bottomBC = 0.0;
		//if(inBottomBC != null)
		bottomBC = inBottomBC.get(1)[0];

		temperatureTopBC = 0.0;
		//if(inTemperatureBottomBC != null)
		temperatureTopBC = inTemperatureTopBC.get(1)[0];
		
		temperatureBottomBC = 0.0;
		//if(inTemperatureBottomBC != null)
		temperatureBottomBC = inTemperatureBottomBC.get(1)[0];

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
					volumes[i] = soilPar.waterContent(psis[i])*dx[i];
					//kappas[i] = soilPar.hydraulicConductivity(psis[i]);
					kappas[i] = soilPar.hydraulicConductivity(psis[i], temperatures[i]);
				} else if(i==NUM_CONTROL_VOLUMES-1) {
					volumes[i] = totalDepth.totalDepth(psis[i]);
					//kappas[i] = soilPar.hydraulicConductivity(psis[i]);
					kappas[i] = soilPar.hydraulicConductivity(psis[i], temperatureTopBC);
				} else {
					volumes[i] = soilPar.waterContent(psis[i])*dx[i];
					//kappas[i] = soilPar.hydraulicConductivity(psis[i]);
					kappas[i] = soilPar.hydraulicConductivity(psis[i], temperatures[i]);
				}
			}
			//k_b = soilPar.hydraulicConductivity(bottomBC);
			k_b = soilPar.hydraulicConductivity(bottomBC, temperatureBottomBC);
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
					sS = sourceSink[i]/tTimestep *timeDelta;
					lowerDiagonal[i] =  bottomBoundaryCondition.lowerDiagonal(-999, kP, kM, spaceDelta[i+1], spaceDelta[i], timeDelta, delta);
					mainDiagonal[i] = bottomBoundaryCondition.mainDiagonal(-999, kP, kM, spaceDelta[i+1], spaceDelta[i], timeDelta, delta);
					upperDiagonal[i] = bottomBoundaryCondition.upperDiagonal(-999, kP, kM, spaceDelta[i+1], spaceDelta[i], timeDelta, delta);
					rhss[i] = volumes[i] + bottomBoundaryCondition.rightHandSide(bottomBC, kP, kM, spaceDelta[i+1], spaceDelta[i], timeDelta, delta)+ timeDelta*sS;


				} else if(i == NUM_CONTROL_VOLUMES -1) {
					kP = 0;	
						// Water flux has to assigned as the minimum between rainfall rate and the maximum infiltrability of the soil
						//tBC = Math.min(topBC, kP);
					kM = 0.5*(kappas[i] + kappas[i-1]);

					sS = sourceSink[i]/tTimestep *timeDelta;
					lowerDiagonal[i] = topBoundaryCondition.lowerDiagonal(-999, kP, kM, spaceDelta[i], spaceDelta[i-1], timeDelta, delta); 
					mainDiagonal[i] = topBoundaryCondition.mainDiagonal(-999, kP, kM, spaceDelta[i], spaceDelta[i-1], timeDelta, delta);
					upperDiagonal[i] = topBoundaryCondition.upperDiagonal(-999, kP, kM,  spaceDelta[i], spaceDelta[i-1], timeDelta, delta);
					rhss[i] = volumes[i] + topBoundaryCondition.rightHandSide(topBC, kP, kM,  spaceDelta[i], spaceDelta[i-1], timeDelta, delta)+ timeDelta*sS;


				} else {

					kP = 0.5*(kappas[i] + kappas[i+1]);
					kM = 0.5*(kappas[i] + kappas[i-1]);

					sS = sourceSink[i]/tTimestep *timeDelta;
					lowerDiagonal[i] = -kM*timeDelta/spaceDelta[i]*1/Math.pow(Math.cos(delta),2);
					mainDiagonal[i] = kM*timeDelta/spaceDelta[i]*1/Math.pow(Math.cos(delta),2)  + kP*timeDelta/spaceDelta[i+1]*1/Math.pow(Math.cos(delta),2);
					upperDiagonal[i] = -kP*timeDelta/spaceDelta[i+1]*1/Math.pow(Math.cos(delta),2);
					rhss[i] = volumes[i] + timeDelta*(kP - kM) + timeDelta*sS; 

				}
				
			}

			//// NESTED NEWTON ALGORITHM ////
			nestedNewtonAlg.set(psis, mainDiagonal, upperDiagonal, lowerDiagonal, rhss);
			psisNew = nestedNewtonAlg.solver();
	
			
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
						velocities[i] =  -kM * (psisNew[i]-bottomBC)/spaceDelta[i] - kM;
						
						volumesNew[i] = soilPar.waterContent(psisNew[i])*dx[i];
						thetas[i] = soilPar.waterContent(psisNew[i]);
					}

				} else if(i == NUM_CONTROL_VOLUMES-1) {
					kP = kappas[i];
					velocities[i] =  -kP * (psisNew[i]-psisNew[i-1])/spaceDelta[i] - kP;
					
					volumesNew[i] = totalDepth.totalDepth(psisNew[i]);
					thetas[i] = totalDepth.totalDepth(psisNew[i]);
				} else {

					kP = 0.5*(kappas[i] + kappas[i+1]);
					velocities[i+1] =  -kP * (psisNew[i+1]-psisNew[i])/spaceDelta[i+1] - kP;
					
					volumesNew[i] = soilPar.waterContent(psisNew[i])*dx[i];
					thetas[i] = soilPar.waterContent(psisNew[i]);
				}
				volume +=volumes[i];
				volumeNew +=volumesNew[i];
			}
			
			errorVolume = volumeNew - volume - timeDelta*(topBC - velocities[0]);
			System.out.println("    errorMass: "+errorVolume);
						
			/**
			 * Solve energy equation 
			 */
			energyIntegrator.set(psis, psisNew, velocities, temperatures, sumTimeDelta, temperatureBottomBC, temperatureTopBC);
			temperatures = energyIntegrator.solver1D();
			//DA RIVEDERE
			energyIntegrator.computeConvectionFluxes(velocities, temperatures, temperatureBottomBC, temperatureTopBC);
			energyIntegrator.computeDiffusionFluxes(psis, temperatures, temperatureBottomBC, temperatureTopBC);
			peclets=energyIntegrator.computePeclets();

			psis = psisNew.clone();
			
		}
		

		//// PRINT OUTPUT FILES ////
		print.setValueFirstVector(depth);
		print.setValueSecondVector(psis);
		print.setValueThirdVector(thetas);
		print.PrintThreeVectors(dir,"Richards_"+step+".csv", inCurrentDate, "Depth[m],Psi[m],Theta[-] ");

		//print.setValueFirstVector(depth);
		//print.setValueSecondVector(velocities);
		//print.PrintTwoVectors(dir,"Flux_"+step+".csv", inCurrentDate, "Depth[m],Velocity[m/s] ");

		//print.setValueFirstVector(depth);
		//print.setValueSecondVector(volumesNew);
		//print.PrintTwoVectors(dir, "Volume_"+step+".csv", inCurrentDate, "Depth[m],Volumes[m] ");
		
		// da sistemare il controllo su quando stampare la temperatura
		if(energySolver.equalsIgnoreCase("PureDiffusion") || energySolver.equalsIgnoreCase("Pure Diffusion") || energySolver.equalsIgnoreCase("diffusion convection")){
			
			print.setValueFirstVector(temperatureDepth);
			print.setValueSecondVector(temperatures);
			print.setValueThirdVector(peclets);
			print.PrintThreeVectors(dir, "Temperature_"+step+".csv", inCurrentDate, "Depth[m],Temperatures[C], Peclets[-] ");
		
		}
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



