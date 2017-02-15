package ic;
import richards_utils.TextIO;

public interface ReadIC {
	public void read(String filePathOrFunction, boolean func);
	public void parse(boolean func);
	public void show();
	public double[] get();
	
}
