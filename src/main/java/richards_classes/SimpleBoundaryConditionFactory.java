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
 * A simple design factory for creating a BoundaryCondition objects.
 */

public class SimpleBoundaryConditionFactory {
	
	/**
	 * Creates a new BoundaryCondition object.
	 * 
	 * @param type of boundary condition
	 * @return boundCond
	 */
	
	public BoundaryCondition createBoundaryCondition (String type) {

		BoundaryCondition boundaryCondition = null;
		if(type.equalsIgnoreCase("Psi Top Dirichlet") || type.equalsIgnoreCase("PsiTopDirichlet")){
			boundaryCondition = new TopBoundaryConditionDirichlet();
		}
		else if(type.equalsIgnoreCase("Psi Top Neumann") || type.equalsIgnoreCase("PsiTopNeumann")){
			boundaryCondition = new TopBoundaryConditionNeumann();
		}
		else if(type.equalsIgnoreCase("Psi Bottom Free Drainage") || type.equalsIgnoreCase("PsiBottomFreeDrainage")){
			boundaryCondition = new BottomBoundaryConditionFreeDrainage();
		}
		else if(type.equalsIgnoreCase("Psi Bottom Dirichlet") || type.equalsIgnoreCase("PsiBottomDirichlet")){
			boundaryCondition = new BottomBoundaryConditionDirichlet();
		}
		else if(type.equalsIgnoreCase("Psi Bottom Impervious") || type.equalsIgnoreCase("PsiBottomImpervious")){
			boundaryCondition = new BottomBoundaryConditionImpervious();
		}
		else if(type.equalsIgnoreCase("Temperature Diffusion Bottom Dirichlet") || type.equalsIgnoreCase("TemperatureDiffusionBottomDirichelet")){
			boundaryCondition = new TemperatureBottomBoundaryConditionDirichlet();
		}
		else if(type.equalsIgnoreCase("Temperature Diffusion Top Dirichlet") || type.equalsIgnoreCase("TemperatureDiffusionTopDirichlet")){
			boundaryCondition = new TemperatureTopBoundaryConditionDirichlet();
		}
		else if(type.equalsIgnoreCase("Temperature Convection Bottom Dirichlet") || type.equalsIgnoreCase("TemperatureConvectionBottomDirichlet")){
			boundaryCondition = new TemperatureConvectionBottomBoundaryConditionDirichlet();
		}
		else if(type.equalsIgnoreCase("Temperature Convection Top Dirichlet") || type.equalsIgnoreCase("TemperatureConvectionTopDirichlet")){
			boundaryCondition = new TemperatureConvectionTopBoundaryConditionDirichlet();
		}
		return boundaryCondition;
		}
}
