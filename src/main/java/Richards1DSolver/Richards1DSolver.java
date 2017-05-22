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

import richards_classes.*;

@Description("Solve the Richards equation for the 1D domain.")
@Documentation("")
@Author(name = "Aaron Iemma, Niccolò Tubini, Michael Dumbser and Riccardo Rigon", contact = "tubini.niccolo@gmail.com")
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
	
	@Description("Time step over mesh space")
	double[] gridvar;

	@Description("Time step over the square of mesh space")
	double[] gridvarsq;

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

	double[] dpsis;

	double[] psis_outer;

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

	@Description("This is the vector which element are those of the principal diagonal and it is used in the Thomas algorithm")
	double[] bb; 

	@Description("This is the vector which element are those of the upper diagonal and it is used in the Thomas algorithm")
	double[] cc;

	@Description("Vector of residuals of the inner iteration")
	double[] fs;

	@Description("Vector of residuals of the inner iteration")
	double[] fks;

	@Description("Correction of the solution after the Newton iterations")
	double[] dis;

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
	double[] topFlux;
	
	Thomas thomasAlg;

	PrintTXT print;

	SoilParametrization soilPar;

	JordanDecomposition jordanDecomposition;

	BoundaryCondition topBoundaryCondition;
	BoundaryCondition bottomBoundaryCondition;

	@Execute
	public void solve() {

		System.out.println("   RICHARDS 1D SOLVER");

		if(step==0){
			NUM_CONTROL_VOLUMES = iC.length;
			//spaceDelta= Math.abs((spaceTop - spaceBottom) / NUM_CONTROL_VOLUMES); 			
			//space_cv_centres= DomainDiscretization.seq(spaceBottom + spaceDelta / 2,spaceTop - spaceDelta / 2,NUM_CONTROL_VOLUMES);
			//gridvar=tTimestep / spaceDelta;
			//gridvarsq=tTimestep / Math.pow(spaceDelta,2);		

			psis 		  = new double[NUM_CONTROL_VOLUMES];	
			dpsis 		  = new double[NUM_CONTROL_VOLUMES];	
			psis_outer    = new double[NUM_CONTROL_VOLUMES];	
			k_t			  = 0.0;								
			k_b			  = 0.0;
			kappas 		  = new double[NUM_CONTROL_VOLUMES];
			thetas 		  = new double[NUM_CONTROL_VOLUMES];
			lowerDiagonal = new double[NUM_CONTROL_VOLUMES];
			mainDiagonal  = new double[NUM_CONTROL_VOLUMES];
			bb            = new double[NUM_CONTROL_VOLUMES]; 
			upperDiagonal = new double[NUM_CONTROL_VOLUMES];
			cc 			  = new double[NUM_CONTROL_VOLUMES];
			fs			  = new double[NUM_CONTROL_VOLUMES];
			fks			  = new double[NUM_CONTROL_VOLUMES];
			dis			  = new double[NUM_CONTROL_VOLUMES];
			rhss 		  = new double[NUM_CONTROL_VOLUMES];
			kP 		      = 0.0;
			kM	          = 0.0;
			outerResidual = 0.0;
			innerResidual = 0.0;
			zeta 		  = new double[NUM_CONTROL_VOLUMES];
			spaceDelta    = new double[NUM_CONTROL_VOLUMES];
			gridvar       = new double[NUM_CONTROL_VOLUMES];
			gridvarsq     = new double[NUM_CONTROL_VOLUMES];
			
			topFlux       = new double[1461];
			thomasAlg = new Thomas();
			print = new PrintTXT();

			SimpleSoilParametrizationFactory soilParFactory = new SimpleSoilParametrizationFactory();
			soilPar = soilParFactory.createSoilParametrization(soilHydraulicModel,alpha,n,psiE,lambda,rMedian,sigma,thetaR,thetaS,ks);

			jordanDecomposition = new JordanDecomposition(soilPar);
			SimpleBoundaryConditionFactory boundCondFactory = new SimpleBoundaryConditionFactory();
			topBoundaryCondition = boundCondFactory.createBoundaryCondition(topBCType);		
			bottomBoundaryCondition = boundCondFactory.createBoundaryCondition(bottomBCType);	
			// Initial domain conditions
			for(int i = 0; i < NUM_CONTROL_VOLUMES; i++) {
				psis[i] = iC[i];
				zeta[i] = spaceBottom-depth[i];
			}
			for(int i = 0; i < NUM_CONTROL_VOLUMES; i++) {
				if (i==0){
					spaceDelta[i] = zeta[i];
				}else if (i== NUM_CONTROL_VOLUMES-1){
					spaceDelta[i] = spaceBottom-zeta[i];
				}else{
					spaceDelta[i] = zeta[i+1]-zeta[i];
				}
			}
			for(int i = 0; i < NUM_CONTROL_VOLUMES; i++) {
				gridvar[i] = tTimestep / spaceDelta[i];
				gridvarsq[i] = tTimestep / Math.pow(spaceDelta[i],2);	
			}
		}
		time = time + tTimestep;

		topBC = 0.0;
		if(inTopBC != null)
			topBC = inTopBC.get(1)[0];///1000;

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
		print.setValueFirstVector(depth);
		print.setValueSecondVector(thetas);
		print.PrintTwoVectors(dir, "Theta_"+step+".csv", "Theta values at time: "+inCurrentDate, "Depth[m],Theta[-] ");
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
				lowerDiagonal[i] =  bottomBoundaryCondition.lowerDiagonal(-999, kP, kM, gridvarsq[i+1], gridvarsq[i], gridvar[i+1], gridvar[i]);
				mainDiagonal[i] = bottomBoundaryCondition.mainDiagonal(-999, kP, kM, gridvarsq[i+1], gridvarsq[i], gridvar[i+1], gridvar[i]);
				upperDiagonal[i] = bottomBoundaryCondition.upperDiagonal(-999, kP, kM, gridvarsq[i+1], gridvarsq[i], gridvar[i+1], gridvar[i]);//-kP * gridvarsq;
				rhss[i] = thetas[i] + bottomBoundaryCondition.rightHandSide(bottomBC, kP, kM, gridvarsq[i+1], gridvarsq[i], gridvar[i+1], gridvar[i]);//  ?gridvar * (kP-kM) + 2 * kM * gridvarsq * psiBottom; 

			} else if(i == NUM_CONTROL_VOLUMES -1) {

				kP = 0.5*(kappas[i] + k_t);
				kM = 0.5*(kappas[i] + kappas[i-1]);
				lowerDiagonal[i] = topBoundaryCondition.lowerDiagonal(-999, kP, kM, gridvarsq[i], gridvarsq[i-1], gridvar[i], gridvar[i-1]); //-kM * gridvarsq;
				mainDiagonal[i] = topBoundaryCondition.mainDiagonal(-999, kP, kM, gridvarsq[i], gridvarsq[i-1], gridvar[i], gridvar[i-1]);//* (kM + 2*kP);
				upperDiagonal[i] = topBoundaryCondition.upperDiagonal(-999, kP, kM, gridvarsq[i], gridvarsq[i-1], gridvar[i], gridvar[i-1]);
				rhss[i] = thetas[i] + topBoundaryCondition.rightHandSide(topBC, kP, kM,  gridvarsq[i], gridvarsq[i-1], gridvar[i], gridvar[i-1]); //gridvar * (kP-kM) + 2 * kP * gridvarsq * psiTop; 
				
				//topFlux[step] = kP*(topBC-psis[i])/gridvar[i];
				
			} else {

				kP = 0.5*(kappas[i] + kappas[i+1]);
				kM = 0.5*(kappas[i] + kappas[i-1]);
				lowerDiagonal[i] = -kM * gridvarsq[i];
				mainDiagonal[i] = gridvarsq[i] * kM + gridvarsq[i+1]*kP;
				upperDiagonal[i] = -kP * gridvarsq[i+1];
				rhss[i] = thetas[i] + gridvar[i+1] * kP - gridvar[i]*kM; 

			}
		}

		// Initial guess of psis
		for(int i = 0; i < NUM_CONTROL_VOLUMES; i++) {
			psis[i] = Math.min(psis[i], soilPar.getPsiStar());
		}

		//// NESTED NEWTON ////
		//// OUTER CYCLE, linearizes one of q_{1}, q_{2} ////
		for(int i = 0; i < MAXITER_NEWT; i++) {
			// I have to assign 0 to outerResidual otherwise I will take into account of the previous error
			outerResidual = 0.0;
			for(int j = 0; j < NUM_CONTROL_VOLUMES; j++) {
				fs[j] = soilPar.waterContent(psis[j]) - rhss[j];
				//fs[j] = jordanDecomposition.waterContent1(psis[j]) - rhss[j];
				if(j == 0) {
					fs[j] = fs[j] + mainDiagonal[j]*psis[j] + upperDiagonal[j]*psis[j+1];
				}else if(j == NUM_CONTROL_VOLUMES -1) {
					fs[j] = fs[j] + lowerDiagonal[j]*psis[j-1] + mainDiagonal[j]*psis[j];
				}else {
					fs[j] = fs[j] + lowerDiagonal[j]*psis[j-1] + mainDiagonal[j]*psis[j] + upperDiagonal[j]*psis[j+1];
				}
				dis[j] = soilPar.dWaterContent(psis[j]);
				outerResidual += fs[j]*fs[j];
			}
			outerResidual = Math.pow(outerResidual,0.5);  
			//System.out.println("Outer iteration " + i + " with residual " +  outerResidual);
			if(outerResidual < newtonTolerance) {
				break;
			}
			if(nestedNewton == 0){
				bb = mainDiagonal.clone();
				cc = upperDiagonal.clone();
				for(int y = 0; y < NUM_CONTROL_VOLUMES; y++) {
					bb[y] += dis[y];
				}
				thomasAlg.set(cc,bb,lowerDiagonal,fs);
				dpsis = thomasAlg.solver();

				//// PSIS UPDATE ////
				for(int s = 0; s < NUM_CONTROL_VOLUMES; s++) {
					psis[s] = psis[s] - dpsis[s];
				}
			}else{
			psis_outer = psis.clone();

			// Initial guess of psis
			for(int j = 0; j < NUM_CONTROL_VOLUMES; j++) {
				psis[j] = Math.max(psis[j], soilPar.getPsiStar());
			}

			//// INNER CYCLE ////
			for(int j = 0; j < MAXITER_NEWT; j++) {
				// I have to assign 0 to innerResidual otherwise I will take into account of the previous error
				innerResidual = 0.0; 
				for(int l=0; l < NUM_CONTROL_VOLUMES; l++) {
					fks[l] = jordanDecomposition.waterContent1(psis[l]) - (jordanDecomposition.waterContent2(psis_outer[l]) + jordanDecomposition.dWaterContent2(psis_outer[l])*(psis[l] - psis_outer[l])) - rhss[l];
					if(l == 0) {
						fks[l] = fks[l] + mainDiagonal[l]*psis[l] + upperDiagonal[l]*psis[l+1];
					}else if(l == NUM_CONTROL_VOLUMES -1) {
						fks[l] = fks[l] + lowerDiagonal[l]*psis[l-1] + mainDiagonal[l]*psis[l];
					}else {
						fks[l] = fks[l] + lowerDiagonal[l]*psis[l-1] + mainDiagonal[l]*psis[l] + upperDiagonal[l]*psis[l+1];
					}
					dis[l] = jordanDecomposition.dWaterContent1(psis[l]) - jordanDecomposition.dWaterContent2(psis_outer[l]);
					innerResidual += fks[l]*fks[l];
				}
				innerResidual = Math.pow(innerResidual,0.5);

				//System.out.println("	Inner iteration " + j + " with residual " +  innerResidual);    

				if(innerResidual < newtonTolerance) {
					break;
				}

				//// THOMAS ////
				// Attention: the main diagonal of the coefficient matrix must not change!! The same for the upper diagonal

				bb = mainDiagonal.clone();
				cc = upperDiagonal.clone();
				for(int y = 0; y < NUM_CONTROL_VOLUMES; y++) {
					bb[y] += dis[y];
				}
				thomasAlg.set(cc,bb,lowerDiagonal,fks);
				dpsis = thomasAlg.solver();

				//// PSIS UPDATE ////
				for(int s = 0; s < NUM_CONTROL_VOLUMES; s++) {
					psis[s] = psis[s] - dpsis[s];
				}
			}
			} //// INNER CYCLE END ////
		} //// OUTER CYCLE END ////


		//// PRINT  ////
		print.setValueFirstVector(depth);
		print.setValueSecondVector(psis);
		//print.PrintTwoVectors(dir,"Psi_"+step+".csv", "Psi values at time: "+inCurrentDate, "Depth[m] Psi[m] ");
		print.PrintTwoVectors(dir,"Psi_"+step+".csv", "Psi values at time: "+inCurrentDate, "Depth[m],Psi[m] ");
		step++;
		
		//print.setValueFirstVector(topFlux);
		//print.PrintOneVector(dir,"TopFlux.csv", "Flux at the top", "Flux");
	} //// MAIN CYCLE END ////
		
}  /// CLOSE Richards1d ///



