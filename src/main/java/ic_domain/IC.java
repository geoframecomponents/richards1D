package ic_domain;

/**
 *	Reads ICs from file, or from given function that must be applied on a domain
 */
public interface IC {
	public void read(String filePathOrFunction, boolean func);
	public void domainSet(double[] domain);
	public void parse(boolean func);
	public void show();
	public double[] get();
	
}
