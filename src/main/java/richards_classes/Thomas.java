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
 * In the case where the matrix A is tridiagonal, i.e. with the special structure there is a very fast and efficient special case of the Gauﬂ-algorithm, called the Thomas-algorithm.
 * Since it is a variation of the Gauﬂ-algorithm, the Thomas algorithm is a direct method to solve general
 * linear tridiagonal systems. As the original Gauﬂ-algorithm, it proceeds in two stages, one forward
 * elimination and one back-substitution.
 * (Numerical Methods for Free Surface Hydrodynamics Solution of linear algebraic equation systems, Dumbser, M., Universit‡ degli studi di Trento)
 * 
 * @author Niccolo' Tubini
 */

public class Thomas {
	
	double[] mainDiagonal;
	double[] upperDiagonal;
	double[] lowerDiagonal;
	double[] rhss;
	double[] solution;
	int DIM;
	

	public Thomas(){}
	
	/**
	 * 
	 * @param upperDiagonal upper diagonal of the coefficient matrix A of the linear system, it is a vector of length NUM_CONTROL_VOLUMES
	 * @param mainDiagonal main diagonal of the coefficient matrix A of the linear system, it is a vector of length NUM_CONTROL_VOLUMES
	 * @param lowerDiagonal lower diagonal of the coefficient matrix A of the linear system, it is a vector of length NUM_CONTROL_VOLUMES
	 * @param rhss right hand side term of the linear system, it is a vector of length NUM_CONTROL_VOLUMES
	 */
	public void set(double[] upperDiagonal, double[] mainDiagonal, double[] lowerDiagonal, double[] rhss){
			this.upperDiagonal = upperDiagonal;
			this.mainDiagonal = mainDiagonal;
			this.lowerDiagonal = lowerDiagonal;
			this.rhss = rhss;
			
			this.DIM = this.mainDiagonal.length;
			this.solution = new double[this.DIM];
	}
	
	
	
	/**
	 * The method solver computes the solution of the linear system with the Thomas algorithm
	 */
	public double[] solver(){
				
		double gamma = 0.0;
		if(this.mainDiagonal.length!= this.upperDiagonal.length | this.mainDiagonal.length!= this.lowerDiagonal.length | this.mainDiagonal.length!= this.rhss.length){
			throw new IllegalArgumentException( "System size error! |n"
					+ "Check the length of diagonal vectors and the right hand side term of the system ");
		}
		

		this.upperDiagonal[0] = this.upperDiagonal[0]/this.mainDiagonal[0];
		this.rhss[0] = this.rhss[0]/this.mainDiagonal[0];
		

		for(int i = 1; i < this.DIM; i++) {
			gamma = 1 / (this.mainDiagonal[i] - this.upperDiagonal[i-1]*this.lowerDiagonal[i]);
			this.upperDiagonal[i] = this.upperDiagonal[i]*gamma;
			this.rhss[i] = (this.rhss[i] - this.lowerDiagonal[i]*this.rhss[i-1])*gamma;
		}

		this.solution[this.DIM-1] = this.rhss[this.DIM-1];

		for(int i=this.DIM-2; i > -1; i--) {
			this.solution[i] = this.rhss[i] - this.upperDiagonal[i]*this.solution[i+1];
		}
		
		return this.solution;

	}

	
	
	/**
	 * This method prints the solution on the screen
	 */
	public void printSolution(){
		System.out.println("The solution computed with Thomas algorithm is: \n");
		for(int i=0;i<this.solution.length;i++){
			System.out.println(this.solution[i]);
		}
	}
}