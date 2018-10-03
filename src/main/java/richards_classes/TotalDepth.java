/*
 * GNU GPL v3 License
 *
 * Copyright 2017 Niccolo Tubini
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
 * 
 * IMPROVEMENT: 'whereas hS is usually set equal to zero; if positive, 
 *              hS represents a small layer of water ponded which can form on 
 *              top of the soil surface during heavy rains before initiation of runoff' (Hydrus1D manual)
 * @author Niccolo' Tubini
 *
 */

public class TotalDepth {

	protected double H;
	protected double dH;
	protected double f1;
	protected double f2;
	protected double df1;
	protected double df2;
	
	
	
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
	
	
	
	/**
	 * @param suction
	 * @return H1 
	 */
	public double pIntegral(double suction){
		this.f1 = totalDepth(suction);

		return this.f1;
	}
	
	
	
	/**
	 * @param suction
	 * @return H2
	 */
	public double qIntegral(double suction){
		this.f2 = 0;

		return this.f2;
	}
	
	
	
	/**
	 * @param suction
	 * @return dH1
	 */
	public double p(double suction){
		this.df1 = dTotalDepth(suction);

		return this.df1;
	}
	
	
	
	/**
	 * @param suction
	 * @return dH2
	 */
	public double q(double suction){
		this.df2 = 0;

		return this.df2;
	}

	
}