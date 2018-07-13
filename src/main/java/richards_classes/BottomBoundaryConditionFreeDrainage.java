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
 * when a free drainage boundary condition is applied at the bottom of the domain.
 * @author Niccolo' Tubini
 *
 */

public class BottomBoundaryConditionFreeDrainage extends BoundaryCondition {

	
	
	public double upperDiagonal(double bC, double kP, double kM, double spaceDeltaP, double spaceDeltaM, double timeDelta, double delta) {
		this.bC = bC;
		this.kP = kP;
		this.kM = kM;
		this.spaceDeltaP = spaceDeltaP;
		this.spaceDeltaM = spaceDeltaM;
		this.timeDelta = timeDelta;
		this.delta = delta;
		
		this.term = -this.kP*this.timeDelta/this.spaceDeltaP*1/Math.pow(Math.cos(this.delta),2);
		
		return term;
	}
	
	
	
	public double mainDiagonal(double bC, double kP, double kM, double spaceDeltaP, double spaceDeltaM, double timeDelta, double delta) {
		this.bC = bC;
		this.kP = kP;
		this.kM = kM;
		this.spaceDeltaP = spaceDeltaP;
		this.spaceDeltaM = spaceDeltaM;
		this.timeDelta = timeDelta;
		this.delta = delta;
		
		this.term = this.kP*this.timeDelta/this.spaceDeltaP*1/Math.pow(Math.cos(this.delta),2);
		
		return term;

	}
	
	
	
	public double lowerDiagonal(double bC, double kP, double kM, double spaceDeltaP, double spaceDeltaM, double timeDelta, double delta) {
		this.bC = bC;
		this.kP = kP;
		this.kM = kM;
		this.spaceDeltaP = spaceDeltaP;
		this.spaceDeltaM = spaceDeltaM;
		this.timeDelta = timeDelta;
		this.delta = delta;
		
		this.term = 0;

		return term;

	}

	
	
	public double rightHandSide(double bC, double kP, double kM, double spaceDeltaP, double spaceDeltaM, double timeDelta, double delta) {
		this.bC = bC;
		this.kP = kP;
		this.kM = kM;
		this.spaceDeltaP = spaceDeltaP;
		this.spaceDeltaM = spaceDeltaM;
		this.timeDelta = timeDelta;
		this.delta = delta;
		
		this.term = this.timeDelta*(this.kP - this.kM); 
		
		return term;

	}	

}
