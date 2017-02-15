package ic;

public abstract class Domain {
	ReadDomainFromFile readdomainfromfile;
	ReadDomainFromUser readdomainfromuser;
	public ReadICFromFile readicfromfile;
	ReadICFromUser readicfromuser;
	String filepath;
	
	public Domain() {
		
	}
	
	public void DisplayDomain() {
		
	}
	
	public void performReadDomainFile(String domainfile) {
		readdomainfromfile.read(domainfile);
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
