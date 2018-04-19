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

import richards_classes.BoundaryCondition;
import thermalParameterization.ThermalSoilParameterization;

/**
 * A simple design factory for creating a SoilParametrization objects.
 * 
 * @author Niccolo' Tubini
 */

public class SimpleEnergyIntegratorFactory {
	
	/**
	 * Creates a new SoilParametrization object.
	 * 
	 * @param type name of the soil thermal model
	 * @param NUM_CONTROL_VOLUMES number of control volume of Richards equation coupled with surface flow
	 * @param dx vector containing length of each control volume
	 * @param spaceDelta vector containing distances between each centroids
	 * @param temperatureTopBoundaryCondition temperature top boundary condition
	 * @param temperatureBottomBoundaryCondition  temperature bottom boundary condition
	 * @return energyIntegrator
	 */

	public EnergyIntegrator createEnergyIntegrator (String type, ThermalSoilParameterization thermalSoilPar, int NUM_CONTROL_VOLUMES, double[] dx, double[] spaceDelta,
			BoundaryCondition temperatureTopBoundaryCondition, BoundaryCondition temperatureBottomBoundaryCondition,
			BoundaryCondition tempConvectionTopBoundaryCondition, BoundaryCondition tempConvectionBottomBoundaryCondition) {

		EnergyIntegrator energyIntegrator = null;
		if(type.equalsIgnoreCase("NoEnergy") || type.equalsIgnoreCase("No energy")){
			energyIntegrator = new NoEnergyIntegrator ();	
		}
		else if(type.equalsIgnoreCase("Pure Diffusion") || type.equalsIgnoreCase("Pure Conduction")){
			energyIntegrator = new EnergyIntegratorPureDiffusion (thermalSoilPar, NUM_CONTROL_VOLUMES, dx, spaceDelta, temperatureTopBoundaryCondition, temperatureBottomBoundaryCondition,
																	tempConvectionTopBoundaryCondition, tempConvectionBottomBoundaryCondition);
		}
		else if(type.equalsIgnoreCase("Conduction convection") || type.equalsIgnoreCase("Diffusion convection")){
			energyIntegrator = new EnergyIntegratorDiffusionConvection (thermalSoilPar, NUM_CONTROL_VOLUMES, dx, spaceDelta, temperatureTopBoundaryCondition, temperatureBottomBoundaryCondition,
																			tempConvectionTopBoundaryCondition, tempConvectionBottomBoundaryCondition);
			
		} else {
			System.out.println("Invalid string for energy integrator");
		}

		return energyIntegrator;
		}
}
