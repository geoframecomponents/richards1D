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
 * when a Dirichlet boundary condition is applied at the top of the domain.
 * @author Niccolò Tubini
 *
 */

public class TemperatureConvectionTopBoundaryConditionDirichlet extends BoundaryCondition{

	
	public double upperDiagonal(double temperatureBC, double velocityP, double velocityM, double spaceDeltaP, double spaceDeltaM, double timeDelta, double delta) {
		this.bC = temperatureBC;
		this.kP = velocityP;
		this.kM = velocityM;
		this.spaceDeltaP = spaceDeltaP;
		this.spaceDeltaM = spaceDeltaM;
		this.timeDelta = timeDelta;
		
		this.term = 0;

		return term;
	}
	
	
	
	public double mainDiagonal(double temperatureBC, double velocityP, double velocityM, double spaceDeltaP, double spaceDeltaM, double timeDelta, double delta) {
		this.bC = temperatureBC;
		this.kP = velocityP;
		this.kM = velocityM;
		this.spaceDeltaP = spaceDeltaP;
		this.spaceDeltaM = spaceDeltaM;
		this.timeDelta = timeDelta;
				
		this.term = this.timeDelta*4188000*0.5*( this.kP+Math.abs(this.kP) - this.kM+Math.abs(this.kM)  ); 
		
		return term;

	}
	
	
	
	public double lowerDiagonal(double temperatureBC, double velocityP, double velocityM, double spaceDeltaP, double spaceDeltaM, double timeDelta, double delta) {
		this.bC = temperatureBC;
		this.kP = velocityP;
		this.kM = velocityM;
		this.spaceDeltaP = spaceDeltaP;
		this.spaceDeltaM = spaceDeltaM;
		this.timeDelta = timeDelta;
		
		
		this.term = timeDelta*4188000*0.5*( -this.kM-Math.abs(this.kM) );

		return term;

	}

	
	
	public double rightHandSide(double temperatureBC, double velocityP, double velocityM, double spaceDeltaP, double spaceDeltaM, double timeDelta, double delta) {
		this.bC = temperatureBC;
		this.kP = velocityP;
		this.kM = velocityM;
		this.spaceDeltaP = spaceDeltaP;
		this.spaceDeltaM = spaceDeltaM;
		this.timeDelta = timeDelta;


		this.term =  -this.timeDelta*4188000*0.5*this.bC*( this.kP-Math.abs(this.kP) );

		return term;

	}
}
