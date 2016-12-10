public class RichardsSolver {
 
	protected Richards createRichards(String item) {
		Richards richards = null;

		if (item.equals("1D")) {
			Richards1DFactory richardsfactory = new Richards1DFactory();
			richards = new Richards1D(richardsfactory);
			richards.setName("one dimensional problem");

		}/* else if (item.equals("2D")) { */

		return richards;
	}

	public Richards solveRichards(String type) {
		Richards richards = createRichards(type);
		System.out.println("Solving a " + richards.getName() + " Richards problem. ");
		richards.solve();
		return richards;
	}
}
