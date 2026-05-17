
/*
 * File Dgemm.java
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
package snap.matrix;

public class Dgemm {
//	static {
//		System.loadLibrary("SSS");
//	}
//	public static native void dgemmX(char arg0, char arg1, int arg2, int arg3, int arg4, double arg5, double[] arg6, int arg7, int arg8, double[] arg9, int arg10, int arg11, double arg12, double[] arg13, int arg14, int arg15);

	public static void dgemm(char arg0, char arg1, int arg2, int arg3, int arg4, double arg5, double[] arg6, int arg7, int arg8, double[] arg9, int arg10, int arg11, double arg12, double[] arg13, int arg14, int arg15) {
		org.netlib.blas.Dgemm.dgemm (""+arg0, ""+arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11, arg12, arg13, arg14, arg15);
	}
} // class Dgemm 
