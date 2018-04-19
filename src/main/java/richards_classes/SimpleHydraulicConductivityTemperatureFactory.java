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
 * A simple design factory for creating a SoilParametrization objects.
 */

public class SimpleHydraulicConductivityTemperatureFactory {
	
	/**
	 * Creates a new SoilParametrization object.
	 * 
	 * @param type name of the model to describe the dependance of hydraulic conductivity on temperature
	 * @return soilPar
	 */
	public HydraulicConductivityTemperature createHydraulicConductivityTemperature (String type) {

		HydraulicConductivityTemperature kappaTemperature = null;
		if(type.equalsIgnoreCase("Hornberger") || type.equalsIgnoreCase("Hornberger1998")){
			kappaTemperature = new Hornberger1998();
		}
		else if(type.equalsIgnoreCase("Isothermal flow") || type.equalsIgnoreCase("Isothermalflow")){
			kappaTemperature = new IsothermalFlow();
		}
		else if(type.equalsIgnoreCase("Ronan")|| type.equalsIgnoreCase("Ronan1998")){
			kappaTemperature = new Ronan1998();
		}
		return kappaTemperature;
		}
}
