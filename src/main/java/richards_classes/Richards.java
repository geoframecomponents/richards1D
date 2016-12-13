package richards_classes;

/*
 * GNU GPL v3 License
 *
 * Copyright 2015 Marialaura Bancheri
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

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