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
 * The boundary condition abstract class.
 * 
 * For the 1D problem the matrix is tridiagonal and thus it is necessary to compute and store
 * only the main, the lower, and the upper diagonal of the matrix.
 * @author Niccolo' Tubini
 *
 */
public abstract class BoundaryCondition {
	
	protected double bC;
	protected double kP;
	protected double kM;
	protected double spaceDeltaP;
	protected double spaceDeltaM;
	protected double timeDelta;
	protected double delta;
	protected double term;
	
	/**
	 * 
	 * This method computes the upper diagonal term of the coefficient matrix T of the system
	 * @param bC value of the boundary condition
	 * @param kP value of the hydraulic conductivity at the interface of volumes (i) and (i+1)
	 * @param kM value of the hydraulic conductivity at the interface of volumes (i) and (i-1)
	 * @param spaceDeltaP distance between nodes of volumes (i) and (i+1), this quantity in notes is referred to as z_n
	 * @param spaceDeltaM distance between nodes of volumes (i) and (i-1), this quantity in notes is referred to as z_s
	 * @param tTimestep time step of the simulation
	 * @param delta length of the control volume
	 * @return
	 */
	public abstract double upperDiagonal(double bC, double kP, double kM, double spaceDeltaP, double spaceDeltaM, double timeDelta, double delta);
	
	/**
	 * 
	 * This method computes the main diagonal term of the coefficient matrix T of the system
	 * @param bC value of the boundary condition
	 * @param kP value of the hydraulic conductivity at the interface of volumes (i) and (i+1)
	 * @param kM value of the hydraulic conductivity at the interface of volumes (i) and (i-1)
	 * @param spaceDeltaP distance between nodes of volumes (i) and (i+1), this quantity in notes is referred to as z_n
	 * @param spaceDeltaM distance between nodes of volumes (i) and (i-1), this quantity in notes is referred to as z_s
	 * @param tTimestep time step of the simulation
	 * @param delta length of the control volume
	 * @return
	 */
	public abstract double mainDiagonal(double bC, double kP, double kM, double spaceDeltaP, double spaceDeltaM, double timeDelta, double delta);
	
	/**
	 *
	 * This method computes the lower diagonal term of the coefficient matrix T of the system
	 * @param bC value of the boundary condition
	 * @param kP value of the hydraulic conductivity at the interface of volumes (i) and (i+1)
	 * @param kM value of the hydraulic conductivity at the interface of volumes (i) and (i-1)
	 * @param spaceDeltaP distance between nodes of volumes (i) and (i+1), this quantity in notes is referred to as z_n
	 * @param spaceDeltaM distance between nodes of volumes (i) and (i-1), this quantity in notes is referred to as z_s
	 * @param tTimestep time step of the simulation
	 * @param delta length of the control volume
	 * @return
	 */
	public abstract double lowerDiagonal(double bC, double kP, double kM, double spaceDeltaP, double spaceDeltaM, double timeDelta, double delta);
	
	/**
	 *
	 * This method computes the right-hand side term of the system
	 * @param bC value of the boundary condition
	 * @param kP value of the hydraulic conductivity at the interface of volumes (i) and (i+1)
	 * @param kM value of the hydraulic conductivity at the interface of volumes (i) and (i-1)
	 * @param spaceDeltaP distance between nodes of volumes (i) and (i+1), this quantity in notes is referred to as z_n
	 * @param spaceDeltaM distance between nodes of volumes (i) and (i-1), this quantity in notes is referred to as z_s
	 * @param tTimestep time step of the simulation
	 * @param delta length of the control volume
	 * @return
	 */
	public abstract double rightHandSide(double bC, double kP, double kM, double spaceDeltaP, double spaceDeltaM, double timeDelta, double delta);
}
