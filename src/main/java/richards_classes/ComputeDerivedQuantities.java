/*
 * GNU GPL v3 License
 *
 * Copyright 2017  Niccolo` Tubini
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

package richards_classes;


/**
 * This class computes derived quantities for Richards problem such as: 
 * 	- velocities at interfaces;
 *  - adimensional water content;
 *  - total water depth;
 *  - water volume
 *  
 *  @author Niccolo' Tubini
 */


public class ComputeDerivedQuantities {

	int NUM_CONTROL_VOLUMES;
	double kM;
	double kP;
	double bottomBC;
	double k_b;
	double volume;

	double[] psis;
	double[] kappas;
	double[] spaceDelta;
	double[] dx;

	double[] thetas;
	double[] volumes;
	double[] darcyVelocities;
	double[] darcyVelocitiesCapillary;
	double[] darcyVelocitiesGravity;
	double[] poreVelocities;
	double[] celerities;       // Rasmussen et al. 2000
	double[] kinematicRatios; // Rasmussen et al. 2000

	String bottomBCType;

	SoilParametrization soilPar;
	TotalDepth totalDepth;
	InterfaceHydraulicConductivity interfaceHydraulicConductivity;



	/**
	 * @param NUM_CONTROL_VOLUMES number of control volumes
	 * @param dx vector containing control volume length
	 * @param spaceDelta vector containing distances between control volumes centroids
	 * @param par1SWRC vector containing first parameter of SWRC model
	 * @param par2SWRC vector containing second parameter of SWRC model
	 * @param thetaR vector containing thetaR
	 * @param thetaS vector containing thetaS
	 * @param soilPar is the class to compute the soil hydraulic properties
	 * @param totalDepth is the class to compute the total water depth
	 */
	public ComputeDerivedQuantities(int NUM_CONTROL_VOLUMES, double[] dx, double[] spaceDelta, SoilParametrization soilPar, TotalDepth totalDepth, InterfaceHydraulicConductivity interfaceHydraulicConductivity, String bottomBCType){

		this.NUM_CONTROL_VOLUMES = NUM_CONTROL_VOLUMES;
		this.spaceDelta = spaceDelta;
		this.dx = dx;
		this.soilPar = soilPar;
		this.totalDepth = totalDepth;
		this.interfaceHydraulicConductivity = interfaceHydraulicConductivity;
		this.bottomBCType = bottomBCType; 

		this.thetas = new double[NUM_CONTROL_VOLUMES];
		this.volumes = new double[NUM_CONTROL_VOLUMES];
		this.kappas = new double[NUM_CONTROL_VOLUMES];
		this.darcyVelocities = new double[NUM_CONTROL_VOLUMES];
		this.darcyVelocitiesCapillary = new double[NUM_CONTROL_VOLUMES];
		this.darcyVelocitiesGravity = new double[NUM_CONTROL_VOLUMES];
		this.poreVelocities = new double[NUM_CONTROL_VOLUMES];
		this.celerities = new double[NUM_CONTROL_VOLUMES];
		this.kinematicRatios = new double[NUM_CONTROL_VOLUMES];
		
	}



	/**
	 * This method is used to update the following variables at each time step.
	 * 
	 * @param psis vector containing all psi values
	 * @param kappas vector containing all hydraulic conductivity values
	 * @param bottomBC water head at the bottom
	 * @param k_b hydraulic conductivity at the bottom
	 */
	public void setComputeDerivedQuantities(double[] psis, double[] kappas, double bottomBC, double k_b) {

		this.bottomBC = bottomBC;
		this.k_b = k_b;
		this.psis = psis;
		this.kappas = kappas;

	}



	/**
	 * This method computes the adimensional water content for those control volumes within the soil
	 * and the total water depth at soil surface
	 * 
	 * @return thetas
	 */
	public double[] computeThetas() {

		for(int i = 0; i < NUM_CONTROL_VOLUMES; i++) {
			if( i == 0 ) {
				thetas[i] = soilPar.waterContent(psis[i],i);
			} else if(i == NUM_CONTROL_VOLUMES-1) {
				thetas[i] = totalDepth.totalDepth(psis[i]);
			} else {
				thetas[i] = soilPar.waterContent(psis[i],i);
			}
		}
		return thetas;	

	}	


	
	/**
	 * This method computes the hydraulic conductivity
	 * 
	 * @return kappas
	 */
	public double[] computeKappas() {

		for(int i = 0; i < NUM_CONTROL_VOLUMES; i++) {

			if(i==NUM_CONTROL_VOLUMES-1) {
				kappas[i] = soilPar.hydraulicConductivity(this.psis[i],i-1);
			} else {
				kappas[i] = soilPar.hydraulicConductivity(this.psis[i],i);
			}

		}
		return kappas;	

	}



	/**
	 * This method computes the water volumes for each control volumes
	 * note that for the 1D case the water volume is a length (volume per unit area)
	 * 
	 * @return volumes
	 */
	public double[] computeWaterVolumes() {

		for(int i = 0; i < NUM_CONTROL_VOLUMES; i++) {
			if( i == 0 ) {
				volumes[i] = soilPar.waterContent(this.psis[i],i)*this.dx[i];
			} else if(i == NUM_CONTROL_VOLUMES-1) {
				volumes[i] = totalDepth.totalDepth(this.psis[i]);
			} else {
				volumes[i] = soilPar.waterContent(this.psis[i],i)*this.dx[i];
			}
		}
		return volumes;	

	}	



	/**
	 * This method computes total water volume 
	 * note that for the 1D case the water volume is a length (volume per unit area)
	 * 
	 * @return volume
	 */
	public double computeTotalWaterVolumes(double[] volumes) {

		volume = 0.0;
		for(int i = 0; i < NUM_CONTROL_VOLUMES; i++) {
			volume += volumes[i];
		}
		return volume;	

	}	



	/**
	 * This method computes velocities at each control volume interface
	 * velocities are computed with Darcy's formula 
	 * 
	 * @return velocities
	 */
	public double[] computeDarcyVelocity() {

		for(int i = 0; i < NUM_CONTROL_VOLUMES; i++) {
			if( i == 0 ) {

				//kP = interfaceHydraulicConductivity.compute(kappas[i],kappas[i+1],dx[i],dx[i+1]);
				if(this.bottomBCType.equalsIgnoreCase("Bottom Free Drainage") || this.bottomBCType.equalsIgnoreCase("BottomFreeDrainage")){
					kM = kappas[i];
					darcyVelocities[i] =  -kM;
				} else if (this.bottomBCType.equalsIgnoreCase("Bottom Impervious") || this.bottomBCType.equalsIgnoreCase("BottomImpervious")) {
					darcyVelocities[i] = + 0;
				} else {
					kM = interfaceHydraulicConductivity.compute(kappas[i],k_b,dx[i],dx[i]);
					darcyVelocities[i] =  -kM * (psis[i]-bottomBC)/spaceDelta[i] - kM;
				}

			} else if(i == NUM_CONTROL_VOLUMES-1) {
				kP = interfaceHydraulicConductivity.compute(kappas[i],kappas[i-1],dx[i],dx[i-1]);
				darcyVelocities[i] =  -kP * (psis[i]-psis[i-1])/spaceDelta[i] - kP;

			} else {
				kM = interfaceHydraulicConductivity.compute(kappas[i],kappas[i-1],dx[i],dx[i-1]);
				darcyVelocities[i] =  -kM * (psis[i]-psis[i-1])/spaceDelta[i] - kM;

			}
		}
		return darcyVelocities;

	}

	

	/**
	 * This method computes velocities due to capillary forces at each control volume interface
	 * 
	 * @return darcyVelocitiesCapillary
	 */
	public double[] computeDarcyVelocityCapillary() {

		for(int i = 0; i < NUM_CONTROL_VOLUMES; i++) {
			if( i == 0 ) {

				//kP = interfaceHydraulicConductivity.compute(kappas[i],kappas[i+1],dx[i],dx[i+1]);
				if(this.bottomBCType.equalsIgnoreCase("Bottom Free Drainage") || this.bottomBCType.equalsIgnoreCase("BottomFreeDrainage")){
					darcyVelocitiesCapillary[i] =  0.0;
				} else if (this.bottomBCType.equalsIgnoreCase("Bottom Impervious") || this.bottomBCType.equalsIgnoreCase("BottomImpervious")) {
					darcyVelocitiesCapillary[i] = + 0.0;
				} else {
					kM = interfaceHydraulicConductivity.compute(kappas[i],k_b,dx[i],dx[i]);
					darcyVelocitiesCapillary[i] =  -kM * (psis[i]-bottomBC)/spaceDelta[i];
				}

			} else if(i == NUM_CONTROL_VOLUMES-1) {
				kP = interfaceHydraulicConductivity.compute(kappas[i],kappas[i-1],dx[i],dx[i-1]);
				darcyVelocitiesCapillary[i] =  -kP * (psis[i]-psis[i-1])/spaceDelta[i];

			} else {
				kM = interfaceHydraulicConductivity.compute(kappas[i],kappas[i-1],dx[i],dx[i-1]);
				darcyVelocitiesCapillary[i] =  -kM * (psis[i]-psis[i-1])/spaceDelta[i];

			}
		}
		return darcyVelocitiesCapillary;

	}

	
	/**
	 * This method computes velocities due to gravity gradient at each control volume interface
	 * 
	 * @return darcyVelocitiesGravity
	 */
	public double[] computeDarcyVelocityGravity() {

		for(int i = 0; i < NUM_CONTROL_VOLUMES; i++) {
			if( i == 0 ) {

				//kP = interfaceHydraulicConductivity.compute(kappas[i],kappas[i+1],dx[i],dx[i+1]);
				if(this.bottomBCType.equalsIgnoreCase("Bottom Free Drainage") || this.bottomBCType.equalsIgnoreCase("BottomFreeDrainage")){
					kM = kappas[i];
					darcyVelocitiesGravity[i] =  -kM;
				} else if (this.bottomBCType.equalsIgnoreCase("Bottom Impervious") || this.bottomBCType.equalsIgnoreCase("BottomImpervious")) {
					darcyVelocitiesGravity[i] = + 0.0;
				} else {
					kM = interfaceHydraulicConductivity.compute(kappas[i],k_b,dx[i],dx[i]);
					darcyVelocitiesGravity[i] =  -kM;
				}

			} else if(i == NUM_CONTROL_VOLUMES-1) {
				kP = interfaceHydraulicConductivity.compute(kappas[i],kappas[i-1],dx[i],dx[i-1]);
				darcyVelocitiesGravity[i] =  -kP;

			} else {
				kM = interfaceHydraulicConductivity.compute(kappas[i],kappas[i-1],dx[i],dx[i-1]);
				darcyVelocitiesGravity[i] =  -kM;

			}
		}
		return darcyVelocitiesGravity;

	}
	
	
	
	/**
	 * This method computes pore velocity at each control volume interface
	 * 
	 * The pore velocity is the Darcy velocity divided by porosity.
	 * FIXME: at the interface of two different soils what is porosity?.
	 * 
	 * @return poreVelocities
	 */
	public double[] computePoreVelocity() {

		for(int i = 0; i < NUM_CONTROL_VOLUMES; i++) {
			if( i == 0 ) {

				//kP = interfaceHydraulicConductivity.compute(kappas[i],kappas[i+1],dx[i],dx[i+1]);
				if(this.bottomBCType.equalsIgnoreCase("Bottom Free Drainage") || this.bottomBCType.equalsIgnoreCase("BottomFreeDrainage")){
					kM = kappas[i];
					poreVelocities[i] =  darcyVelocities[i]/(soilPar.thetaS[i]-soilPar.thetaR[i]);
				} else if (this.bottomBCType.equalsIgnoreCase("Bottom Impervious") || this.bottomBCType.equalsIgnoreCase("BottomImpervious")) {
					poreVelocities[i] = + 0;
				} else {
					//kM = interfaceHydraulicConductivity.compute(thetas[i],soilPar.waterContent(bottomBC,i),dx[i],dx[i]);
					poreVelocities[i] =  darcyVelocities[i]/(soilPar.thetaS[i]-soilPar.thetaR[i]);
				}

			} else if(i == NUM_CONTROL_VOLUMES-1) {
				//kP = interfaceHydraulicConductivity.compute(kappas[i],kappas[i-1],dx[i],dx[i-1]);
				poreVelocities[i] =  darcyVelocities[i]/(soilPar.thetaS[i-1]-soilPar.thetaR[i-1]);

			} else {
				//kM = interfaceHydraulicConductivity.compute(thetas[i],thetas[i-1],dx[i],dx[i-1]);
				poreVelocities[i] = darcyVelocities[i]/(soilPar.thetaS[i]-soilPar.thetaR[i]);

			}
		}
		return poreVelocities;

	}
	
	
	
	/**
	 * This method computes celerity of pressure wave at each control volume interface
	 * 
	 * @return celerities
	 */
	public double[] computeCelerity() {

		for(int i = 0; i < NUM_CONTROL_VOLUMES; i++) {
			if( i == 0 ) {

				//kP = interfaceHydraulicConductivity.compute(kappas[i],kappas[i+1],dx[i],dx[i+1]);
				if(this.bottomBCType.equalsIgnoreCase("Bottom Free Drainage") || this.bottomBCType.equalsIgnoreCase("BottomFreeDrainage")){
					//kM = kappas[i];
					celerities[i] =  -soilPar.dHydraulicConductivity(psis[i], i);
				} else if (this.bottomBCType.equalsIgnoreCase("Bottom Impervious") || this.bottomBCType.equalsIgnoreCase("BottomImpervious")) {
					celerities[i] = + 0.0;
				} else {
					kM = interfaceHydraulicConductivity.compute(soilPar.dHydraulicConductivity(psis[i], i),soilPar.dHydraulicConductivity(bottomBC, i),dx[i],dx[i]);
					celerities[i] =  -kM * ( (psis[i]-bottomBC)/spaceDelta[i] +1 );
				}

			} else if(i == NUM_CONTROL_VOLUMES-1) {
				kP = interfaceHydraulicConductivity.compute(soilPar.dHydraulicConductivity(psis[i-1], i-1),soilPar.dHydraulicConductivity(psis[i-1], i-1),dx[i],dx[i-1]);
				celerities[i] =  -kP * ( (psis[i]-psis[i-1])/spaceDelta[i] +1 );

			} else {
				kM = interfaceHydraulicConductivity.compute(soilPar.dHydraulicConductivity(psis[i], i),soilPar.dHydraulicConductivity(psis[i-1], i-1),dx[i],dx[i-1]);
				celerities[i] =  -kM * ( (psis[i]-psis[i-1])/spaceDelta[i] +1 ) ;

			}
		}
		return celerities;

	}	
	
	/**
	 * This method computes the ratio between celerity and pore velocity 
	 * at each control volume interface
	 * 
	 * @return kinematicRatios
	 */
	public double[] computeKinematicRatio() {

		for(int i = 0; i < NUM_CONTROL_VOLUMES; i++) {
			if( i == 0 ) {

				if (this.bottomBCType.equalsIgnoreCase("Bottom Impervious") || this.bottomBCType.equalsIgnoreCase("BottomImpervious")) {
					kinematicRatios[i] = Double.NaN;
				} else {
					kinematicRatios[i] =  celerities[i]/poreVelocities[i];
				}

			} else {
				kinematicRatios[i] =  celerities[i]/poreVelocities[i];
			}
		}
		return kinematicRatios;

	}

}






