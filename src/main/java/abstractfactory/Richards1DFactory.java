public class Richards1DFactory implements RichardsFactory {

	// Variable handling
	public Theta createTheta() {
		return new Theta1D();
	}

	public DTheta createDTheta() {
		return new DTheta1D();
	}

	public Theta1 createTheta1() {
		return new Theta11D();
	}

	public DTheta1 createDTheta1() {
		return new DTheta11D();
	}
	public Theta2 createTheta2() {
		return new Theta21D();
	}
	public DTheta2 createDTheta2() {
		return new DTheta21D();
	}


	// "Real" methods!
	public Kappa createKappa() {
		return new Kappa1D();
	}
	public MatOP createMatOP() {
		return new MatOP1D();
	}
	public ConjGrad createConjGrad() {
		return new ConjGrad1D();
	}

}