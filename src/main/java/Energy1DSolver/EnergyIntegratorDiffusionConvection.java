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

package Energy1DSolver;


/**
 * This class allows to solve energy equation considering both conduction and convection.
 * Conduction (diffusion part) is discretized in an implicit way, whilst convection is explicit in time.
 * This require the definition of a time step restriction but allow to have a matrix A that is symmetric. 
 * 
 * @author Niccolo' Tubini
 */

import richards_classes.BoundaryCondition;
import richards_classes.SoilParametrization;
import thermalParameterization.ThermalSoilParameterization;

public class EnergyIntegratorDiffusionConvection extends EnergyIntegrator{
	
	
	public EnergyIntegratorDiffusionConvection (ThermalSoilParameterization thermalSoilPar, int NUM_CONTROL_VOLUMES, double[] dx, double[] spaceDelta,
			BoundaryCondition tempTopBoundaryCondition, BoundaryCondition tempBottomBoundaryCondition,
			BoundaryCondition tempConvectionTopBoundaryCondition, BoundaryCondition tempConvectionBottomBoundaryCondition) {
		
		super( thermalSoilPar, NUM_CONTROL_VOLUMES, dx, spaceDelta, tempTopBoundaryCondition, tempBottomBoundaryCondition,
				tempConvectionTopBoundaryCondition, tempConvectionBottomBoundaryCondition);
		System.out.println("Energy equation with conduction and convection.");
		
	}
	
	@Override
	public double[] solver1D() {
		internalEnergy = super.computeTotalInternalEnergy(psis, temperatures);
		
		for(int i=0; i<NUM_CONTROL_VOLUMES; i++) {
			if(i==0) {
				
				lambdaP = 0.5*( super.thermalSoilPar.thermalConductivity(super.psis[i]) + super.thermalSoilPar.thermalConductivity(super.psis[i+1]) );
				lambdaM = 0.5*( super.thermalSoilPar.thermalConductivity(super.psis[i]) + super.thermalSoilPar.thermalConductivity(super.psis[i]) );

				upperDiagonal[i] = tempBottomBoundaryCondition.upperDiagonal(-999, lambdaP, lambdaM, spaceDelta[i+1], spaceDelta[i],timeDelta,-999)
									+tempConvectionBottomBoundaryCondition.upperDiagonal(-999, velocities[i+1], velocities[i], spaceDelta[i+1],spaceDelta[i],timeDelta,-999);
				mainDiagonal[i] = super.thermalSoilPar.heatCapacity(super.psisNew[i])*dx[i] + tempBottomBoundaryCondition.mainDiagonal(-999, lambdaP, lambdaM, spaceDelta[i+1], spaceDelta[i],timeDelta,-999)
									+tempConvectionBottomBoundaryCondition.mainDiagonal(-999, velocities[i+1], velocities[i], spaceDelta[i+1],spaceDelta[i],timeDelta,-999);
				lowerDiagonal[i] = tempBottomBoundaryCondition.lowerDiagonal(-999, lambdaP, lambdaM, spaceDelta[i+1], spaceDelta[i],timeDelta,-999)
									+tempConvectionBottomBoundaryCondition.lowerDiagonal(-999, velocities[i+1], velocities[i], spaceDelta[i+1],spaceDelta[i],timeDelta,-999);
				rhss[i] = super.thermalSoilPar.heatCapacity(super.psis[i])*temperatures[i]*dx[i] + tempBottomBoundaryCondition.rightHandSide(temperatureBottomBC, lambdaP, lambdaM, spaceDelta[i+1], spaceDelta[i],timeDelta,-999)
									+tempConvectionBottomBoundaryCondition.rightHandSide(temperatureBottomBC, velocities[i+1], velocities[i], spaceDelta[i+1],spaceDelta[i],timeDelta,-999);;
			
				
			} else if(i==NUM_CONTROL_VOLUMES-1) {
		
				lambdaP = 0.5*( super.thermalSoilPar.thermalConductivity(super.psis[i]) + super.thermalSoilPar.thermalConductivity(super.psis[i]) );
				lambdaM = 0.5*( super.thermalSoilPar.thermalConductivity(super.psis[i]) + super.thermalSoilPar.thermalConductivity(super.psis[i-1]) );

				upperDiagonal[i] = tempTopBoundaryCondition.upperDiagonal(-999, lambdaP, lambdaM, spaceDelta[i+1], spaceDelta[i],timeDelta,-999)
									+tempConvectionTopBoundaryCondition.upperDiagonal(-999, velocities[i+1], velocities[i], spaceDelta[i+1], spaceDelta[i],timeDelta,-999);
				mainDiagonal[i] = super.thermalSoilPar.heatCapacity(super.psisNew[i])*dx[i] + tempTopBoundaryCondition.mainDiagonal(-999, lambdaP, lambdaM, spaceDelta[i+1], spaceDelta[i],timeDelta,-999)
									+tempConvectionTopBoundaryCondition.mainDiagonal(-999, velocities[i+1], velocities[i], spaceDelta[i+1], spaceDelta[i],timeDelta,-999);
				lowerDiagonal[i] = tempTopBoundaryCondition.lowerDiagonal(-999, lambdaP, lambdaM, spaceDelta[i+1], spaceDelta[i],timeDelta,-999)
									+tempConvectionTopBoundaryCondition.lowerDiagonal(-999, velocities[i+1], velocities[i], spaceDelta[i+1], spaceDelta[i],timeDelta,-999);
				rhss[i] = super.thermalSoilPar.heatCapacity(super.psis[i])*temperatures[i]*dx[i] + tempTopBoundaryCondition.rightHandSide(temperatureTopBC, lambdaP, lambdaM, spaceDelta[i+1], spaceDelta[i],timeDelta,-999)
									+tempConvectionTopBoundaryCondition.rightHandSide(temperatureTopBC, velocities[i+1], velocities[i], spaceDelta[i+1], spaceDelta[i],timeDelta,-999);
			} else {
			
				lambdaP = 0.5*( super.thermalSoilPar.thermalConductivity(super.psis[i]) + super.thermalSoilPar.thermalConductivity(super.psis[i+1]) );
				lambdaM = 0.5*( super.thermalSoilPar.thermalConductivity(super.psis[i]) + super.thermalSoilPar.thermalConductivity(super.psis[i-1]) );

				upperDiagonal[i] = -timeDelta/spaceDelta[i+1]*lambdaP + timeDelta*4188000*( 0.5*( velocities[i+1]-Math.abs(velocities[i+1]) )  );
				mainDiagonal[i] = super.thermalSoilPar.heatCapacity(super.psisNew[i])*dx[i] + timeDelta/spaceDelta[i]*lambdaM + timeDelta/spaceDelta[i+1]*lambdaP + timeDelta*4188000*( 0.5*( velocities[i+1]+Math.abs(velocities[i+1])  - velocities[i]+Math.abs(velocities[i]) )  );
				lowerDiagonal[i] = -timeDelta/spaceDelta[i]*lambdaM + timeDelta*4188000*( 0.5*( -velocities[i]-Math.abs(velocities[i]) )  );
				rhss[i] = super.thermalSoilPar.heatCapacity(super.psis[i])*temperatures[i]*dx[i];
				
			}
		}
				
		tempThomasAlg.set(upperDiagonal, mainDiagonal, lowerDiagonal, rhss);
		super.temperatures = tempThomasAlg.solver();
		internalEnergyNew = super.computeTotalInternalEnergy(psisNew, temperatures);
		
		errorEnergy = internalEnergyNew - internalEnergy + timeDelta*( 4188000*0.5*( velocities[NUM_CONTROL_VOLUMES]*(temperatureTopBC+temperatures[NUM_CONTROL_VOLUMES-1]) - Math.abs(velocities[NUM_CONTROL_VOLUMES])*(temperatureTopBC-temperatures[NUM_CONTROL_VOLUMES-1]) )
																			-super.thermalSoilPar.thermalConductivity(super.psis[NUM_CONTROL_VOLUMES-1])*(temperatureTopBC-temperatures[NUM_CONTROL_VOLUMES-1])/spaceDelta[NUM_CONTROL_VOLUMES] 
																	   -4188000*0.5*( velocities[0]*(temperatures[0]+temperatureBottomBC) - Math.abs(velocities[0])*(temperatures[0]-temperatureTopBC) )					
																	   + super.thermalSoilPar.thermalConductivity(super.psis[0])*(temperatures[0]-temperatureBottomBC)/spaceDelta[0]);
		
		/*for(int i=0; i<NUM_CONTROL_VOLUMES; i++) {
			System.out.println(i+" "+super.temperatures[i]);
		}*/	
		
		System.out.println(" \n====\n   errorEnergy: "+errorEnergy);
		
		return super.temperatures;
	}

}
