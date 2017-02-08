package ic;

public class Domain1D extends Domain{
	public Domain1D() {
		readdomainfromfile = new ReadDomainFile1D();
		readdomainfromuser = new ReadDomainUser1D();
		readicfromfile = new ReadICFile1D();
		readicfromuser = new ReadICUser1D();
		
	}
	
	public void display() {
		System.out.println("Works?");
	}
}
