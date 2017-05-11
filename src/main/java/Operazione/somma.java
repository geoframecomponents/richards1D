package Operazione;

import oms3.annotations.*;

public class somma {
	
	@In
	public double value;
	
    private final double addendo=1000;
	private double somma;
    
    @Execute
    public void process(){
    	somma = value+addendo;
    	System.out.println("somma: " +somma);
    }
}
