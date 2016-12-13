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

public class Richards1D extends Richards {
	Richards1DFactory richardsfactory;

	public Richards1D(Richards1DFactory richardsfactory) {
		this.richardsfactory = richardsfactory;
	}

	void solve() {
		System.out.println("Solving for " + name);
		theta = richardsfactory.createTheta();
		dtheta = richardsfactory.createDTheta();
		theta1 = richardsfactory.createTheta1();
		dtheta1 = richardsfactory.createDTheta1();
		theta2 = richardsfactory.createTheta2();
		dtheta2 = richardsfactory.createDTheta2();

		kappa = richardsfactory.createKappa();
		matop = richardsfactory.createMatOP();
		conjgrad = richardsfactory.createConjGrad();

	}
}