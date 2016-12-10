import cern.colt.Arrays;
import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.algo.solver.DoubleBiCG;
import cern.colt.matrix.tdouble.algo.solver.DoubleCG;
import cern.colt.matrix.tdouble.algo.solver.DoubleCGICCTest;
import cern.colt.matrix.tdouble.algo.solver.DoubleCGLS;
import cern.colt.matrix.tdouble.algo.solver.DoubleCGS;
import cern.colt.matrix.tdouble.algo.solver.IterativeSolverDoubleNotConvergedException;
import cern.colt.matrix.tdouble.algo.solver.preconditioner.DoubleDiagonal;
import cern.colt.matrix.tdouble.algo.solver.preconditioner.DoubleICC;
import cern.colt.matrix.tdouble.algo.solver.preconditioner.DoubleILU;
import cern.colt.matrix.tdouble.algo.solver.preconditioner.DoubleILUT;
import cern.colt.matrix.tdouble.algo.solver.preconditioner.DoublePreconditioner;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix1D;
import cern.colt.matrix.tdouble.impl.SparseDoubleMatrix1D;
import cern.colt.matrix.tdouble.impl.SparseRCDoubleMatrix2D;


public class MatrixOps {
	public static void main(String[] args) {

		public DoubleMatrix1D matSol;
		DoublePreconditioner dd;
		DoubleCG conjugateGradient;
		
		/*
		* Instantiates a new conjugate gradient 
		*/
		public ConjugateGradient(int SIZE) {

			mat_x = new DenseDoubleMatrix1D(SIZE);
			matSol = new DenseDoubleMatrix1D(SIZE);
			dd = new DoubleDiagonal(SIZE);

			conjugateGradient = new DoubleCG(mat_x);

		}

		/**
		 * Method solverCG
		 * 
		 * @param matrix_a DoubleMatrix1D 
		 * @param matrix_b SParseDoubleMatrix1D
		 * @throws IterativeSolverDoubleNotConvergedException
		 *
		 */
		public void solverCG(DoubleMatrix1D mat_b, SparseDoubleMatrix2D mat_a) throws IterativeSolverDoubleNotConvergedException {

			dd.setMatrix(matrix_a);
			conjugateGradient.setPreconditioner(dd);
			sol = conjugateGradient.solve(mat_a, mat_b, mat_x);
			
		}

		/**
		 * The main method.
		 * 
		 * @param args
		 *            the arguments
		 * @throws IterativeSolverDoubleNotConvergedException 
		 */
		public static void main(String[] args) throws IterativeSolverDoubleNotConvergedException {

			double[] b = { 1, 2 };
			int[] Mp = { 0, 2, 4 };
			int[] Mi = { 0, 1, 0, 1 };
			double[] Ml = { 4, 1, 1, 3 };

			DenseDoubleMatrix1D matrix_b = new DenseDoubleMatrix1D(b);
			SparseDoubleMatrix2D matrix_A = new SparseDoubleMatrix2D(b.length, b.length, Mp, Mi, Ml);
			ConjugateGradient cg = new RCConjugateGradient(b.length);

			cg.solverCG(matrix_b, matrix_A);
			
			System.out.println(cg.matSol.get(0));

		}
}		
