public class Richards1D extends Richards {
	Richards1DFactory richardsfactory;

	public Richards1D(Richards1DFactory richardsfactory) {
		this.richardsfactory = richardsfactory;
	}

	void solve() {
		System.out.println("Solving for " + name);
		theta = richardsfactory.createTheta();
		dtheta = richardsfactory.createDTheta();
		theta1 = richardsfactory.createTheta1();
		dtheta1 = richardsfactory.createDTheta1();
		theta2 = richardsfactory.createTheta2();
		dtheta2 = richardsfactory.createDTheta2();

		kappa = richardsfactory.createKappa();
		matop = richardsfactory.createMatOP();
		conjgrad = richardsfactory.createConjGrad();

	}
}