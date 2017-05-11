package Operazione;

//import java.awt.List;

public class somma2 {
	
	public int[] vettore;
	public int[] somma;
	
	public int[] sum(int[] input){
		this.vettore = input;
		this.somma = new int[this.vettore.length];
		for(int i=0; i<this.vettore.length; i++) {
			somma[i] = vettore[i] + vettore[i] +100;
		}
		
		return somma;
	}
}
