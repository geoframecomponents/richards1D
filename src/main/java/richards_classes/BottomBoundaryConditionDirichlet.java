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
 * This class compute the element of the coefficient matrix and the right-hand side
 * when a Dirichlet boundary condition is applied at the bottom of the domain.
 * @author Niccolò Tubini
 *
 */
public class BottomBoundaryConditionDirichlet extends BoundaryCondition {

	
	
	public double upperDiagonal(double bC, double kP, double kM, double spaceDeltaP, double spaceDeltaM, double tTimestep, double delta) {
		this.bC = bC;
		this.kP = kP;
		this.kM = kM;
		this.spaceDeltaP = spaceDeltaP;
		this.spaceDeltaM = spaceDeltaM;
		this.tTimestep = tTimestep;
		this.delta = delta;
		
		this.term =-this.kP * this.tTimestep/(this.spaceDeltaM+this.spaceDeltaP/2)*1/this.spaceDeltaP*1/Math.pow(Math.cos(this.delta),2);
		
		return term;
	}
	
	
	
	public double mainDiagonal(double bC, double kP, double kM, double spaceDeltaP, double spaceDeltaM, double tTimestep, double delta) {
		this.bC = bC;
		this.kP = kP;
		this.kM = kM;
		this.spaceDeltaP = spaceDeltaP;
		this.spaceDeltaM = spaceDeltaM;
		this.tTimestep = tTimestep;
		this.delta = delta;
		
		this.term = this.tTimestep/(this.spaceDeltaM+this.spaceDeltaP/2)*1/this.spaceDeltaM*1/Math.pow(Math.cos(this.delta),2)* this.kM + this.tTimestep/(this.spaceDeltaM+this.spaceDeltaP/2)*1/this.spaceDeltaP*1/Math.pow(Math.cos(this.delta),2)*this.kP;

		return term;

	}
	
	
	
	public double lowerDiagonal(double bC, double kP, double kM, double spaceDeltaP, double spaceDeltaM, double tTimestep, double delta) {
		this.bC = bC;
		this.kP = kP;
		this.kM = kM;
		this.spaceDeltaP = spaceDeltaP;
		this.spaceDeltaM = spaceDeltaM;
		this.tTimestep = tTimestep;
		this.delta = delta;
		
		this.term = 0;

		return term;

	}

	
	
	public double rightHandSide(double bC, double kP, double kM, double spaceDeltaP, double spaceDeltaM, double tTimestep, double delta) {
		this.bC = bC;
		this.kP = kP;
		this.kM = kM;
		this.spaceDeltaP = spaceDeltaP;
		this.spaceDeltaM = spaceDeltaM;
		this.tTimestep = tTimestep;
		this.delta = delta;
		
		this.term = this.tTimestep/(this.spaceDeltaM+this.spaceDeltaP/2)* this.kP - this.tTimestep/(this.spaceDeltaM+this.spaceDeltaP/2)*this.kM + this.kM * this.tTimestep/(this.spaceDeltaM+this.spaceDeltaP/2)*1/this.spaceDeltaM*1/Math.pow(Math.cos(this.delta),2)* this.bC;

		return term;

	}
}
