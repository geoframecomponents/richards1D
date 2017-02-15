package ic;

public abstract class Domain {
	public ReadDomain readdomain;
	public ReadIC readic;
	String filepath;
	String filePathOrFunction;
	
	public void performReadDomain(String filePath) {
		readdomain.read(filePath);
		readdomain.parse();
		readdomain.get();
	}
	public void performReadIC(String filePathOrFunction, boolean func) {
		readic.read(filePathOrFunction, func);
		readic.parse(func);
		readic.get();
	}

}
