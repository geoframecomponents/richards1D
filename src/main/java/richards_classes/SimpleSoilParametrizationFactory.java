package richards_classes;

public class SimpleSoilParametrizationFactory {
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
