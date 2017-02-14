package ic;

public abstract class Domain {
	ReadDomainFromFile readdomainfromfile;
	ReadDomainFromUser readdomainfromuser;
	ReadICFromFile readicfromfile;
	ReadICFromUser readicfromuser;
	String filepath;
	
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
	public void performReadICFile(String filepath) {
		readicfromfile.read(filepath);
	}
	public void performReadICUser() {
		readicfromuser.read();
	}

}
