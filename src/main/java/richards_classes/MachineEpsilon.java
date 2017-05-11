package richards_classes;

public class MachineEpsilon {

	/**
	 * Calculate machine epsilon double.
	 *
	 * this method compute the tolerance of the machine.
	 * @see <a href="https://en.wikipedia.org/wiki/Machine_epsilon#Approximation_using_Java">
	 *     https://en.wikipedia.org/wiki/Machine_epsilon#Approximation_using_Java</a>
	 * <blockquote>In such languages as C or C++ when you do something like
	 * while(1.0 + eps &gt; 1.0) the expression is calculated not with 64 bits
	 * (double) but with the processor precision (80 bits or more depends on
	 * the processor and compile options). Below program calculates exactly on
	 * 32 bits (float) and 64 bits (double)</blockquote>
	 *
	 * @return the tolerance of the machine for double
	 */
	public static double doublePrecision() {

		double machEps = 1.0d;

		do
			machEps /= 2.0d;
		while ((double) (1.0 + (machEps / 2.0)) != 1.0);

		return machEps;
	}

	/**
	 * Calculate machine epsilon float.
	 *
	 * this method compute the tolerance of the machine.
	 * @see <a href="https://en.wikipedia.org/wiki/Machine_epsilon#Approximation_using_Java">
	 *     https://en.wikipedia.org/wiki/Machine_epsilon#Approximation_using_Java</a>
	 * <blockquote>In such languages as C or C++ when you do something like
	 * while(1.0 + eps &gt; 1.0) the expression is calculated not with 64 bits
	 * (double) but with the processor precision (80 bits or more depends on
	 * the processor and compile options). Below program calculates exactly on
	 * 32 bits (float) and 64 bits (double)</blockquote>
	 *
	 * @return the tolerance of the machine for float
	 */
	public static float floatPrecision() {

		float machEps = 1.0f;

		do
			machEps /= 2.0f;
		while ((float) (1.0 + (machEps / 2.0)) != 1.0);

		return machEps;
	}

}