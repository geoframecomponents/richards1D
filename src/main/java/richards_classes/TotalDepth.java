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
 * The soil prametrization abstract class.
 * @author Niccolo' Tubini
 *
 */

public class TotalDepth {

	public double H;
	public double dH;
	/**
	 * This method compute the total water depth given the suction value 
	 * @param suction
	 * @return
	 */
	public double totalDepth(double suction) {
		
		if (suction>0) {
			H = suction;
		} else {
			H = 0;
		}
		
		return H;
	};
	
	
	/**
	 * This method compute the value of the derivative of theta given the suction value
	 * @param suction
	 * @return
	 */
	public double dTotalDepth(double suction) {
		
		if (suction>0) {
			dH = 1;
		} else {
			dH = 0;
		}
		
		return dH;
	}
	
	

}