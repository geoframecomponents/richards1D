package ic;

public class Domain1D extends Domain{
	public Domain1D() {
		readdomain = new ReadDomain1D();
		readic = new ReadIC1D();	
	}
	
	public void display() {
		System.out.println("This is from Domain1D, extension of Domain class");
	}
}
