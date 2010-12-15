
/*
 * File ExpQT.java
 *
 * Copyright (C) 2010 Remco Bouckaert, David Bryant remco@cs.auckland.ac.nz
 *
 * This file is part of SnAP.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * SnAP is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  SnAP is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SnAP; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */
package snap.likelihood;

import snap.matrix.Array2d;

public class ExpQT {

	native static void expQTtx(int N, double u, double v, double coalescenceRate, double t, double [] x);
	native static void expM(double [] x, int n);

	static void expM(Array2d A, Array2d F) {
		int n = A.getNrOfRows();
		double [] x = new double[n*n];
		System.arraycopy(A.asZeroBasedArray(), 0, x, 0, n*n);
		expM(x, n);
		F.set(x, n);
	} // expM

}
