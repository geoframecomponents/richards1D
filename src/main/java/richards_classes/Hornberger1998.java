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
 * @author Niccolo` Tubini
 *
 */
public class Hornberger1998 extends HydraulicConductivityTemperature{
	
	
	@Override
	public double temperatureCorrection(double temperature) {
		
		corrector =  1.5869*Math.pow(10, -4)*Math.pow(temperature, 2) + 2.5263*Math.pow(10,-2)*temperature + 0.7315; 
		
		return corrector;
	}

}
