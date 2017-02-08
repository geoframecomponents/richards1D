package ic;

public abstract class Domain {
	ReadDomainFromFile readdomainfromfile;
	ReadDomainFromUser readdomainfromuser;
	ReadICFromFile readicfromfile;
	ReadICFromUser readicfromuser;
	
	public Domain() {
		
	}
	
	public void DisplayDomain() {
		
	}
	
	public void performReadDomainFile() {
		readdomainfromfile.read();
	}
	public void performReadDomainUser() {
		readdomainfromuser.read();
	}
	public void performReadICFile() {
		readicfromfile.read();
	}
	public void performReadICUser() {
		readicfromuser.read();
	}

}
