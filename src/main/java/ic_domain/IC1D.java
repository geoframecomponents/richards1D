package ic_domain;
import java.util.*;


import richards_utils.Expr;
import richards_utils.TextIO;

public class IC1D implements IC {
	public String filePathOrFunction	= null;
	public double[] domain				= null;
	private Expr fun 					= null;
	public double[] ic 					= null;
	
	public void domainSet(double[] domain) {
		this.domain = domain;
	}
	
	public void read(String filePathOrFunction, boolean func) {
		
		this.filePathOrFunction = filePathOrFunction;

		if(func) {
				try {
					fun = new Expr(filePathOrFunction);
				}
				catch (IllegalArgumentException e) {
					System.out.println("Error! The definition of f(x) is not valid.");
					System.out.println(e.getMessage());
				}
		} else {
			try{
				TextIO.readFile(filePathOrFunction);
			}
			catch (IllegalArgumentException e) {
				System.out.println("Can't open file " + filePathOrFunction + " for reading!");
				System.exit(1);  
			}		
		}
	}
	
	public void parse(boolean func) {		
		if(func) {
			try{
				for(int i=0;i<domain.length;i++) {
					this.ic[i] = fun.value(domain[i]);
				}				
			} catch (IllegalArgumentException e) {
				System.out.println("Function or domain undefined!");
				System.exit(1);  	
			}			
		} else {
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
		    double[] arr = new double[centres.size()];
		    for(int i=0;i<centres.size();i++) {
		    	arr[i] = centres.get(i);
		    }
		    
		    if(this.ic.length==this.domain.length) {
			    this.ic = arr;		    	
		    } else {
				System.out.println("IC and domain files are not of equal length!");
				System.exit(1);  		    	
		    }
		    
		}    
	}
	
	public void show() {
		try{
			for(int i=0;i<ic.length;i++) {
				System.out.println(ic[i]);
			}
		}
		catch (NullPointerException e) {
			System.out.println("You must parse() the given file first!");
			System.exit(1);  			
		}		
	}
	
	public double[] get() {
		
		return ic;
	}
	
	public double[] applyIC() {
		return null;
	}
}
