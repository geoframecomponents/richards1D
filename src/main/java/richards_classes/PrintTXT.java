/**
 * 
 */
package richards_classes;
import java.io.*;

/**
 * @author Niccolò
 *	This class allows to print a .txt file with one or two or three vectors
 *	Calling methods Print you can pass the name of the file, the header and the name of variables
 */
public class PrintTXT {
	// Vectors to print
	private double v1[];
	private double v2[];
	private double v3[];
	private double v4[];
	
	// Matrix to print
	private double matrix[][];
	
	// File name, header, name of variables
	public String name;
	public String header;
	public String variable;
	
	public void setValueFirstVector(double a[]){
		this.v1 = a;
	}
	
	public void setValueSecondVector(double a[]){
		this.v2 = a;
	}
	
	public void setValueThirdVector(double a[]){
		this.v3 = a;
	}
	
	public void setValueFourthVector(double a[]){
		this.v4 = a;
	}
	
	public void setValueMatrix(double a[][]){
		this.matrix = a;
	}
	
	public void PrintOneVector(String dir, String name, String header, String variable){
		this.name = name;
		this.header = header;
		this.variable = variable;
		
		try {
		      //FileWriter fileout = new FileWriter(name);
		      FileWriter fileout = new FileWriter(new File(dir, name));
		      fileout.write(this.header +"\n");
		      fileout.write(this.variable +"\n");
		      for(int i=0;i<v1.length;i++)
		      {
		    	  	fileout.write(v1[i] +"\n");
		    	  	
		      }
		            
		      fileout.close(); // chiude il file
		  }
		      catch (IOException e)
		      {
		          System.out.println("Errore: " + e);
		          System.exit(1);
		      }
	}
	
	
	public void PrintTwoVectors(String dir, String name, String header, String variable){
		this.name = name;
		this.header = header;
		this.variable = variable;
		
		try {
		      //FileWriter fileout = new FileWriter(name);
		      FileWriter fileout = new FileWriter(new File(dir, name));
		      fileout.write(this.header +"\n");
		      fileout.write(this.variable +"\n");
		      for(int i=0;i<v1.length;i++)
		      {
		    	  	fileout.write(v1[i] + "," + v2[i] +"\n");
		    	  	
		      }
		            
		      fileout.close(); // chiude il file
		  }
		      catch (IOException e)
		      {
		          System.out.println("Errore: " + e);
		          System.exit(1);
		      }
	}
	
	
	public void PrintThreeVectors(String dir, String name, String header, String variable){
		this.name = name;
		this.header = header;
		this.variable = variable;
		
		try {
		      //FileWriter fileout = new FileWriter(name);
		      FileWriter fileout = new FileWriter(new File(dir, name));
		      fileout.write(this.header +"\n");
		      fileout.write(this.variable +"\n");
		      for(int i=0;i<v1.length;i++)
		      {
		    	  	fileout.write(v1[i] + "," + v2[i] + "," + v3[i] +"\n");
		    	
		      }
		            
		      fileout.close(); // chiude il file
		  }
		      catch (IOException e)
		      {
		          System.out.println("Errore: " + e);
		          System.exit(1);
		      }
	}
	
	
	
	public void PrintFourVectors(String dir, String name, String header, String variable){
		this.name = name;
		this.header = header;
		this.variable = variable;
		
		try {
		      //FileWriter fileout = new FileWriter(name);
		      FileWriter fileout = new FileWriter(new File(dir, name));
		      fileout.write(this.header +"\n");
		      fileout.write(this.variable +"\n");
		      for(int i=0;i<v1.length;i++)
		      {
		    	  	fileout.write(v1[i] + "," + v2[i] + "," + v3[i] + "," + v4[i]+"\n");
		    	
		      }
		            
		      fileout.close(); // chiude il file
		  }
		      catch (IOException e)
		      {
		          System.out.println("Errore: " + e);
		          System.exit(1);
		      }
	}
	
	
	
	
	public void PrintMatrix(String dir, String name, String header, String variable){
		this.name = name;
		this.header = header;
		this.variable = variable;
		int j;
		try
	    {
	      FileWriter writer = new FileWriter(new File(dir, name));  
	      writer.write(this.header+"\n");
	      writer.write(this.variable+"\n");
	      //writer.flush();
	         for(int i = 0; i < matrix.length; i++)
	         {
	            for (j=0; j<4; j++)
	             {
	                 writer.append(String.valueOf(matrix[i][j]));
	                 writer.append(',');
	             }
	               writer.append(String.valueOf(matrix[i][j])+"\n");
	               //writer.append('\n');
	              // writer.flush();
	         }
	         writer.close();
	      }        
	    catch(Exception e)
	    {
	    	//System.out.println("Errore: " + e);
	        //System.exit(1);
	        e.printStackTrace();
	    }
	}
}