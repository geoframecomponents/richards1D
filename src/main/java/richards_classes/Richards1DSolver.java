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
