public interface RichardsFactory {
	public Theta createTheta();
	public DTheta createDTheta();
	public Theta1 createTheta1();
	public DTheta1 createDTheta1();
	public Theta2 createTheta2();
	public DTheta2 createDTheta2();

	public Kappa createKappa();
	public MatOP createMatOP();
	public ConjGrad createConjGrad();

}