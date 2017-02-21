package ic_domain;

public interface Domain {
	public void read(String filePath);
	public void parse();
	public void show();
	public double[] get();
}
