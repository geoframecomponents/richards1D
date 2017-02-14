package ic;
import java.util.*;
import richards_utils.TextIO;

public class ReadICFile1D implements ReadICFromFile {
	public String filepath;
	public double[] IC;
	
	public void read(String filepath) {
		
		this.filepath = filepath;
		
		try{
			TextIO.readFile(filepath);
		}
		catch (IllegalArgumentException e) {
			System.out.println("Can't open file " + filepath + " for reading!");
			System.exit(1);  
		}		
	}
	
	public double[] parse() {
		List<Double> centres = new ArrayList<Double>();

		try{
			TextIO.readFile(filepath);
		}
		catch (IllegalArgumentException e) {
			System.out.println("You must set a file into stream with read() before parsing it!");
			System.exit(1);  
		}
		
		while ( ! TextIO.eof() ) {
	    	centres.add(TextIO.getDouble());  
	    }
	    double[] arr = centres.stream().mapToDouble(Double::doubleValue).toArray(); // method reference!
	    return arr;
	}
	
	public void show() {
		
	}
}
