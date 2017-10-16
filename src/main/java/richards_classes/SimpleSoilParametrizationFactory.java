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

public class SimpleSoilParametrizationFactory {
	
	/**
	 * Creates a new SoilParametrization object.
	 * 
	 * @param type name of the SWRC model
	 * @param alpha SWRC parameter
	 * @param n SWRC parameter
	 * @param psi_e SWRC parameter
	 * @param lambda SWRC parameter
	 * @param rMedian SWRC parameter
	 * @param sigma SWRC parameter
	 * @param theta_r dimensionless residual water content
	 * @param theta_s dimensionless water content at saturation
	 * @param ks hydraulic conductivity at saturation
	 * @return soilPar
	 */
	public SoilParametrization createSoilParametrization (String type, double alpha, double n, double psi_e, double lambda, double rMedian, double sigma, double theta_r, double theta_s, double ks) {

		SoilParametrization soilPar = null;
		if(type.equalsIgnoreCase("BrooksCorey") || type.equalsIgnoreCase("Brooks Corey")){
			soilPar = new BrooksCorey(lambda,psi_e,theta_r,theta_s,ks);
		}
		else if(type.equalsIgnoreCase("VanGenuchten") || type.equalsIgnoreCase("Van Genuchten")){
			soilPar = new VanGenuchten(n,alpha,theta_r,theta_s,ks);
		}
		else if(type.equalsIgnoreCase("Kosugi")){
			soilPar = new KosugiUnimodal(rMedian,sigma,theta_r,theta_s,ks);
		}
		else if(type.equalsIgnoreCase("Gardner")){
			soilPar = new Gardner(alpha,theta_r,theta_s,ks);
		}
		return soilPar;
		}
}
