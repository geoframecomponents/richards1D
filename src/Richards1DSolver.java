public class Richards1DSolver extends RichardsSolver {

	protected Richards createRichards(String item) {
		Richards richards = null;

		if (item.equals("1D")) {
			Richards1DFactory richardsfactory = new Richards1DFactory();
			richards = new Richards1D(richardsfactory);
			richards.setName("Richards 1-D solver");
			richards.setDimension(1);

		}/* else if (item.equals("2D")) { */

		return richards;
	}
}
