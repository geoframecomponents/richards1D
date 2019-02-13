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
 * A simple design factory to create an InterfaceHydraulicConductivity object.
 * @author Niccolo' Tubini
 */

public class SimpleInterfaceHydraulicConductivityFactory {


	/**
	 * Creates a new InterfaceHydraulicConductivity object.
	 * 
	 * @param type of method to compute the hydraulic conductivity at interfaces
	 * @return InterfaceHydraulicConductivity
	 */

	public InterfaceHydraulicConductivity createInterfaceHydraulicConductivity (String type) {

		InterfaceHydraulicConductivity interfaceHydraulicConductivity = null;
		if(type.equalsIgnoreCase("Average") || type.equalsIgnoreCase("Mean")){
			interfaceHydraulicConductivity = new InterfaceHydraulicConductivityAverage();
		}
		else if(type.equalsIgnoreCase("Max") || type.equalsIgnoreCase("Maximum")){
			interfaceHydraulicConductivity = new InterfaceHydraulicConductivityMax();
		}
		else if(type.equalsIgnoreCase("Min") || type.equalsIgnoreCase("Minimum")){
			interfaceHydraulicConductivity = new InterfaceHydraulicConductivityMin();
		}
		else if(type.equalsIgnoreCase("WeightedAverage") || type.equalsIgnoreCase("Weighted Average")){
			interfaceHydraulicConductivity = new InterfaceHydraulicConductivityWeightedAverage();
		}

		return interfaceHydraulicConductivity;

	}

}
