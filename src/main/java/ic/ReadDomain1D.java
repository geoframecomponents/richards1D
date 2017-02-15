package ic;
import java.util.ArrayList;

import java.util.List;

import richards_utils.Expr;
import richards_utils.TextIO;

public class ReadDomain1D implements ReadDomain {
	public String filePathOrFunction;
	private double[] domain = null;
	private int domainLength;
	
	public void read(String filePath) {
		
		this.filePathOrFunction = filePath;
			try{
				TextIO.readFile(filePath);
			} catch (IllegalArgumentException e) {
				System.out.println("Can't open file " + filePath + " for reading!");
				System.exit(1);  
			}
	}
	
	public void parse() {
		List<Double> centres = new ArrayList<Double>();

		try{
			TextIO.readFile(filePathOrFunction);
		}
		catch (IllegalArgumentException e) {
			System.out.println("You must set a file into stream with read() before parsing it!");
			System.exit(1);  
		}
		
		while ( ! TextIO.eof() ) {
	    	centres.add(TextIO.getlnDouble());  
	    }
	    double[] arr = centres.stream().mapToDouble(Double::doubleValue).toArray(); // method reference!
	    this.domain = arr;
	}
	
	public void show() {
		try{
			for(int i=0;i<domain.length;i++) {
				System.out.println(domain[i]);
			}
		}
		catch (NullPointerException e) {
			System.out.println("You must parse() the given file first!");
			System.exit(1);  			
		}		
	}
	public double[] get(){
		return domain;			
	}
	
}
