/*
 * GNU GPL v3 License
 *
 * Copyright 2017 
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

	double[] thetaS;
	double[] thetaR;
	double[] par1SWRC;
	double[] par2SWRC;
	double[] alphaSpecificStorage;
	double[] betaSpecificStorage;
	double[] ks;

	double[] thetas;
	double[] volumes;
	double[] velocities;

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
	public ComputeDerivedQuantities(int NUM_CONTROL_VOLUMES, double[] dx, double[] spaceDelta, double[] par1SWRC, double[] par2SWRC, double[] alphaSpecificStorage, double[] betaSpecificStorage, double[] thetaR, double[] thetaS, double[] ks,
			SoilParametrization soilPar, TotalDepth totalDepth, InterfaceHydraulicConductivity interfaceHydraulicConductivity, String bottomBCType){

		this.NUM_CONTROL_VOLUMES = NUM_CONTROL_VOLUMES;
		this.spaceDelta = spaceDelta;
		this.dx = dx;
		this.par1SWRC = par1SWRC;
		this.par2SWRC = par2SWRC;
		this.alphaSpecificStorage = alphaSpecificStorage;
		this.betaSpecificStorage = betaSpecificStorage;
		this.thetaR = thetaR;
		this.thetaS = thetaS;
		this.ks = ks;
		this.soilPar = soilPar;
		this.totalDepth = totalDepth;
		this.interfaceHydraulicConductivity = interfaceHydraulicConductivity;
		this.bottomBCType = bottomBCType; 

		this.thetas = new double[NUM_CONTROL_VOLUMES];
		this.volumes = new double[NUM_CONTROL_VOLUMES];
		this.kappas = new double[NUM_CONTROL_VOLUMES];
		this.velocities = new double[NUM_CONTROL_VOLUMES+1];
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
				soilPar.set(par1SWRC[i],par2SWRC[i],alphaSpecificStorage[i],betaSpecificStorage[i],thetaR[i],thetaS[i],-999);
				thetas[i] = soilPar.waterContent(psis[i]);
			} else if(i == NUM_CONTROL_VOLUMES-1) {
				thetas[i] = totalDepth.totalDepth(psis[i]);
			} else {
				soilPar.set(par1SWRC[i],par2SWRC[i],alphaSpecificStorage[i],betaSpecificStorage[i],thetaR[i],thetaS[i],-999);
				thetas[i] = soilPar.waterContent(psis[i]);
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
				soilPar.set(par1SWRC[i-1],par2SWRC[i-1],alphaSpecificStorage[i-1],betaSpecificStorage[i-1],thetaR[i-1],thetaS[i-1],ks[i-1]);
				kappas[i] = soilPar.hydraulicConductivity(this.psis[i]);
			} else {
				soilPar.set(par1SWRC[i],par2SWRC[i],alphaSpecificStorage[i],betaSpecificStorage[i],thetaR[i],thetaS[i],ks[i]);
				kappas[i] = soilPar.hydraulicConductivity(this.psis[i]);
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
				soilPar.set(par1SWRC[i],par2SWRC[i],alphaSpecificStorage[i],betaSpecificStorage[i],thetaR[i],thetaS[i],-999);
				volumes[i] = soilPar.waterContent(this.psis[i])*this.dx[i];
			} else if(i == NUM_CONTROL_VOLUMES-1) {
				volumes[i] = totalDepth.totalDepth(this.psis[i]);
			} else {
				soilPar.set(par1SWRC[i],par2SWRC[i],alphaSpecificStorage[i],betaSpecificStorage[i],thetaR[i],thetaS[i],-999);
				volumes[i] = soilPar.waterContent(this.psis[i])*this.dx[i];
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
	public double[] computeVelocities() {

		for(int i = 0; i < NUM_CONTROL_VOLUMES; i++) {
			if( i == 0 ) {

				//kP = 0.5*(kappas[i] + kappas[i+1]);
				kP = interfaceHydraulicConductivity.compute(kappas[i],kappas[i+1],dx[i],dx[i+1]);
				if(this.bottomBCType.equalsIgnoreCase("Bottom Free Drainage") || this.bottomBCType.equalsIgnoreCase("BottomFreeDrainage")){
					kM = kappas[i];
					velocities[i] =  -kM;
				} else if(this.bottomBCType.equalsIgnoreCase("Bottom Impervious") || this.bottomBCType.equalsIgnoreCase("BottomImpervious")) {
					velocities[i] = + 0;

				}
				else {
					//kM = 0.5*(kappas[i] + k_b);
					kM = interfaceHydraulicConductivity.compute(kappas[i],k_b,dx[i],dx[i]);
					velocities[i] =  -kM * (psis[i]-bottomBC)/spaceDelta[i] - kM;

				}

			} else if(i == NUM_CONTROL_VOLUMES-1) {
				kP = kappas[i];
				velocities[i] =  -kP * (psis[i]-psis[i-1])/spaceDelta[i] - kP;

			} else {

				//kP = 0.5*(kappas[i] + kappas[i+1]);
				kP = interfaceHydraulicConductivity.compute(kappas[i],kappas[i+1],dx[i],dx[i+1]);
				velocities[i] =  -kP * (psis[i+1]-psis[i])/spaceDelta[i+1] - kP;

			}
		}
		return velocities;

	}




}






