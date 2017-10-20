package richards_classes;

import oms3.annotations.*;
import java.io.*;
import java.util.*;

/**
 * Read csv two columns
 *
 * @author Niccolo' Tubini, Francesco Serafin
 */
public class ReadCsvTwoColumns {
	@Description("Path of the input file")
	@In
	public String filePath;
	
	@Description("Depth")
	@Out
	public double[] depth;
	@Description("Variable: generic term to describe either initial condition for psi or the source/sink term")
	@Out
	public double[] variable;
	
	@Execute
	public void process() throws IOException {
		
		List<Double> l1 = new ArrayList<>();
		List<Double> l2 = new ArrayList<>();
		
		BufferedReader fileReader = new BufferedReader(new FileReader(filePath));
		String eachLine = "";
		while ((eachLine = fileReader.readLine()) != null) {
            String[] values = eachLine.split(",");
            l1.add(Double.parseDouble(values[0]));
		    l2.add(Double.parseDouble(values[1]));
		}

		fileReader.close();
        
		depth = new double[l1.size()];
		variable = new double[l2.size()];
		
		for(int i = 0;i < l1.size();i++){
			depth[i] = l1.get(i);
			variable[i] = l2.get(i);
		}

	}

	public void setFilePath(String filePath) {
		this.filePath = filePath;
	}

	public double[] getDepth() {
		return depth;
	}

	public double[] getVariable() {
		return variable;
	}

}

