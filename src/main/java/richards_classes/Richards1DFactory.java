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

public class Richards1DFactory implements RichardsFactory {

	// Variable handling
	public Theta createTheta() {
		return new Theta1D();
	}

	public DTheta createDTheta() {
		return new DTheta1D();
	}

	public Theta1 createTheta1() {
		return new Theta11D();
	}

	public DTheta1 createDTheta1() {
		return new DTheta11D();
	}
	public Theta2 createTheta2() {
		return new Theta21D();
	}
	public DTheta2 createDTheta2() {
		return new DTheta21D();
	}


	// "Real" methods!
	public Kappa createKappa() {
		return new Kappa1D();
	}
	public MatOP createMatOP() {
		return new MatOP1D();
	}
	public ConjGrad createConjGrad() {
		return new ConjGrad1D();
	}

}