public abstract class Richards {
	String name;
	int dimension;

	Theta theta;
	DTheta dtheta;
	Theta1 theta1;
	DTheta1 dtheta1;
	Theta2	theta2;
	DTheta2 dtheta2;

	ConjGrad conjgrad;
	MatOP matop;
	Kappa kappa;

	abstract void solve();

	void setName(String name) {
		this.name = name;
	}
	void setDimension(int dimension) {
		this.dimension = dimension;
	}

	String getName() {
		return name;
	}
	int getDimension() {
		return dimension;
	}

	public String toString() {
		StringBuffer result = new StringBuffer();
		result.append("---- " + dimension + " ----\n");
		if (theta1 != null) {
			result.append(theta1);
			result.append("\n");
		}
		if (theta2 != null) {
			result.append(theta2);
			result.append("\n");
		}
		if (dtheta1 != null) {
			result.append(dtheta1);
			result.append("\n");
		}
		if (dtheta2 != null) {
			result.append(dtheta2);
			result.append("\n");
		}
		if (theta != null) {
			result.append(theta);
			result.append("\n");
		}
		return result.toString();
	}	
}