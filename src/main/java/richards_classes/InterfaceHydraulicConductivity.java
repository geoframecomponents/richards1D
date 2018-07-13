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
 * Compute hydraulic conductivity at control volume interface
 * It is possible to compute hydraulic conductivity K(i+1/2):
 * - average of k[i] and k[i+1]
 * - max(k[i],k[i+1])
 * - min(k[i],k[i+1])
 * - weighted average of k[i] and k[i+1]
 * @author Niccolo' Tubini
 */
public abstract class InterfaceHydraulicConductivity {
	
	double kappa;
	
	public abstract double compute(double kappa1, double kappa2, double dx1, double dx2);

}
