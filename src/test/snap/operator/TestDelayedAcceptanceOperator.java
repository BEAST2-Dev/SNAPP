package test.snap.operator;

import org.junit.Test;

import beast.core.Distribution;
import beast.core.parameter.RealParameter;
import beast.core.util.CompoundDistribution;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Sequence;
import beast.evolution.operators.ScaleOperator;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.tree.Tree;
import beast.util.TreeParser;

import snap.likelihood.SnapSubstitutionModel;
import snap.operators.DelayedAcceptanceOperator;
//import test.beast.BEASTTestCase;

import junit.framework.TestCase;

public class TestDelayedAcceptanceOperator extends TestCase {

	@Test
	public void testMomentGeneratingFunction() throws Exception {
		// check that moment generating function = 1 at x = 0 = -2*(u+v)
		snap.Data data = getTwoTaxaNoData();
		Tree tree = getTree(data, "(A:" + 0.1 +",B:" + 0.1 +")");
		ScaleOperator dummyoperator = new ScaleOperator();
		Distribution prior = new CompoundDistribution();
		
		
		RealParameter uParameter = new RealParameter(0.0 + "");
		RealParameter vParameter = new RealParameter(0.0 + "");
		Double [] coalescenceRate = new Double[]{0.1,0.2,0.3};
		RealParameter coalescenceRateParameter = new RealParameter(coalescenceRate);
		
		SnapSubstitutionModel substModel = new SnapSubstitutionModel();
		substModel.initByName("mutationRateU", uParameter,
				"mutationRateV", vParameter,
				"coalescenceRate", coalescenceRateParameter);
		SiteModel siteModel = new SiteModel();
		siteModel.initByName("substModel", substModel);
		
		DelayedAcceptanceOperator operator = new DelayedAcceptanceOperator();
		operator.initByName("tree", tree,
				"data", data,
				"operator", dummyoperator,
				"prior", prior,
				"siteModel", siteModel);

		double [] M = new double[3];
		// u = v = 0
		operator.calcMomentGeneratingFunction(M, tree.getRoot(), 0.0, 0.0, coalescenceRate);
		assertEquals(1.0, M[0]);
		assertEquals(1.0, M[1]);
		assertEquals(1.0, M[2]);
	}
	
    Tree getTree(Alignment data, String tree) throws Exception {
        TreeParser t = new TreeParser();
        t.initByName("taxa", data,
                "newick", tree);
        return t;
    }

	@Test
	public void testApproximateLikelihood() throws Exception {
/*
 		double [][] mu = new double[2][2];
 

		mu[0][0] = 0.5544555487941784;
		mu[0][1] = 0.5591458864149094;
		mu[1][0] = mu[0][1];
		mu[1][1] = 0.5498091674862808;
		
		testSingleCase(0.1,  // t
			1, // u
			1, // v
			new Double[]{0.1,  0.2, 0.3}, // 2/theta 
			mu, 
			0.1);
*/
		
		double u;
		double v;
		double t;
		Double [] coalescenceRate = new Double[3];
		double [][] expectedMu = new double[2][2];
		
		for (int i = 0; i < data.length; i += 9) {
			u = data[i];
			v = data[i+1];
			t = data[i+2];
			coalescenceRate[0] = 2.0/data[i+3];
			coalescenceRate[1] = 2.0/data[i+4];
			coalescenceRate[2] = 2.0/data[i+5];
			expectedMu[0][0] = data[i+6];
			expectedMu[0][1] = data[i+7];
			expectedMu[1][0] = data[i+7];
			expectedMu[1][1] = data[i+8];
			testSingleCase(t, u, v, coalescenceRate, expectedMu, 0.0);
		}
		
	}
	
    snap.Data getTwoTaxaNoData() throws Exception {
        Sequence a = new Sequence("A", "0");
        Sequence b = new Sequence("B", "1");
 
        snap.Data data = new snap.Data();
        data.initByName("sequence", a, "sequence", b, "dataType", "binary");
        return data;
    }
	
	void testSingleCase(double t, double u, double v, Double [] coalescenceRate, double[][] expectedMu, double expectedLogL) throws Exception {
		snap.Data data = getTwoTaxaNoData();
		Tree tree = getTree(data, "(A:" + t +",B:" + t +")");
		ScaleOperator dummyoperator = new ScaleOperator();
		Distribution prior = new CompoundDistribution();
		
		
		RealParameter uParameter = new RealParameter(u + "");
		RealParameter vParameter = new RealParameter(v + "");
		RealParameter coalescenceRateParameter = new RealParameter(coalescenceRate);
		
		SnapSubstitutionModel substModel = new SnapSubstitutionModel();
		substModel.initByName("mutationRateU", uParameter,
				"mutationRateV", vParameter,
				"coalescenceRate", coalescenceRateParameter);
		SiteModel siteModel = new SiteModel();
		siteModel.initByName("substModel", substModel);
		
		DelayedAcceptanceOperator operator = new DelayedAcceptanceOperator();
		operator.initByName("tree", tree,
				"data", data,
				"operator", dummyoperator,
				"prior", prior,
				"siteModel", siteModel);
		
		operator.useMatLabFormulae = true;
		
		double [][] mu = operator.calcMu(u, v, coalescenceRate);
		for (int i = 0; i < mu.length; i++) {
			for (int j = 0; j < mu.length; j++) {
				assertEquals(expectedMu[i][j], mu[i][j], 1e-5);
			}
		}
		
		//double approxLogL = operator.approxLikelihood();
		//assertEquals(expectedLogL, approxLogL, 1e-10);
	}
	

double [] data = new double [] {
		    //  Us =,
			//	Vs =,
			//	Ts =,
			//	Thetas =,
			//	Distances =,

				0.83436,
				1.24769,
				0.58225,
				6.4432e-03, 5.9619e-04, 4.2289e-03,
				6.3579e-03, 4.3815e-01, 5.9545e-04,
				2.09857,
				0.65639,
				0.54074,
				 3.7861e-03, 6.8197e-03, 9.4229e-04,
				 3.7470e-03, 3.4458e-01, 6.6940e-03,
				1.13243,
				0.8953,
				0.869942,
				 8.1158e-03, 4.2431e-04, 5.9852e-03,
				 7.9844e-03, 4.7886e-01, 4.2395e-04,
				1.56707,
				0.73429,
				0.26478,
				 5.3283e-03, 7.1445e-04, 4.7092e-03,
				 5.2637e-03, 3.0745e-01, 7.1328e-04,
				0.8413,
				1.2325,
				0.318075,
				 3.5073e-03, 5.2165e-03, 6.9595e-03,
				 3.4819e-03, 3.5513e-01, 5.1607e-03,
				1.71396,
				0.70594,
				0.119216,
				 9.3900e-03, 9.6730e-04, 6.9989e-03,
				 9.1814e-03, 1.8503e-01, 9.6504e-04,
				1.03594,
				0.96647,
				0.93983,
				 8.7594e-03, 8.1815e-03, 6.3853e-03,
				 8.6084e-03, 4.8796e-01, 8.0496e-03,
				1.81816,
				0.68966,
				0.645553,
				 5.5016e-03, 8.1755e-03, 3.3604e-04,
				 5.4267e-03, 3.8312e-01, 8.0112e-03,
				1.88843,
				0.68006,
				0.479464,
				 6.2248e-03, 7.2244e-03, 6.8806e-04,
				 6.1268e-03, 3.5623e-01, 7.0928e-03,
				2.0063,
				0.66597,
				0.639318,
				 5.8704e-03, 1.4987e-03, 3.1960e-03,
				 5.7798e-03, 3.6204e-01, 1.4927e-03,
				1.41108,
				0.7744,
				0.544717,
				 2.0774e-03, 6.5961e-03, 5.3086e-03,
				 2.0680e-03, 4.1574e-01, 6.5023e-03,
				0.67764,
				1.90732,
				0.647312,
				 3.0125e-03, 5.1859e-03, 6.5445e-03,
				 2.9892e-03, 3.7346e-01, 5.1173e-03,
				0.96795,
				1.03424,
				0.543887,
				 4.7092e-03, 9.7297e-03, 4.0762e-03,
				 4.6652e-03, 4.4333e-01, 9.5438e-03,
				2.33667,
				0.63612,
				0.721048,
				 2.3049e-03, 6.4899e-03, 8.1998e-03,
				 2.2892e-03, 3.3187e-01, 6.3671e-03,
				0.81476,
				1.29427,
				0.522496,
				 8.4431e-03, 8.0033e-03, 7.1836e-03,
				 8.2954e-03, 4.2260e-01, 7.8705e-03,
				2.16163,
				0.65045,
				0.993706,
				 1.9476e-03, 4.5380e-03, 9.6865e-03,
				 1.9370e-03, 3.5431e-01, 4.4808e-03,
				1.58668,
				0.73006,
				0.218678,
				 2.2592e-03, 4.3239e-03, 5.3133e-03,
				 2.2475e-03, 2.7684e-01, 4.2810e-03,
				2.50227,
				0.62486,
				0.105799,
				 1.7071e-03, 8.2531e-03, 3.2515e-03,
				 1.6980e-03, 1.5644e-01, 8.0455e-03,
				0.66635,
				2.00285,
				0.109698,
				 2.2766e-03, 8.3470e-04, 1.0563e-03,
				 2.2629e-03, 1.6664e-01, 8.3284e-04,
				1.39536,
				0.77922,
				0.063592,
				 4.3570e-03, 1.3317e-03, 6.1096e-03,
				 4.3161e-03, 1.1568e-01, 1.3279e-03,
				0.72331,
				1.61954,
				0.404581,
				 3.1110e-03, 1.7339e-03, 7.7880e-03,
				 3.0885e-03, 3.6387e-01, 1.7269e-03,
				2.4338,
				0.62928,
				0.448374,
				 9.2338e-03, 3.9094e-03, 4.2345e-03,
				 8.9798e-03, 3.0580e-01, 3.8631e-03,
				0.51927,
				13.47458,
				0.365817,
				 4.3021e-03, 8.3138e-03, 9.0823e-04,
				 4.0578e-03, 7.1457e-02, 7.4474e-03,
				2.05982,
				0.66027,
				0.763506,
				 1.8482e-03, 8.0336e-03, 2.6647e-03,
				 1.8389e-03, 3.6190e-01, 7.8618e-03,
				2.14461,
				0.65201,
				0.627897,
				 9.0488e-03, 6.0471e-04, 1.5366e-03,
				 8.8255e-03, 3.4695e-01, 6.0369e-04,
				2.24739,
				0.64307,
				0.771981,
				 9.7975e-03, 3.9926e-03, 2.8101e-03,
				 9.5277e-03, 3.4201e-01, 3.9470e-03,
				0.67887,
				1.89765,
				0.932855,
				 4.3887e-03, 5.2688e-03, 4.4009e-03,
				 4.3396e-03, 3.8498e-01, 5.1982e-03,
				1.30957,
				0.80881,
				0.972742,
				 1.1112e-03, 4.1680e-03, 5.2714e-03,
				 1.1086e-03, 4.6449e-01, 4.1315e-03,
				1.02974,
				0.97193,
				0.192029,
				 2.5806e-03, 6.5686e-03, 4.5742e-03,
				 2.5674e-03, 2.7008e-01, 6.4834e-03,
				2.11014,
				0.65527,
				0.138875,
				 4.0872e-03, 6.2797e-03, 8.7537e-03,
				 4.0415e-03, 1.9783e-01, 6.1725e-03,
				1.37283,
				0.78643,
				0.696267,
				 5.9490e-03, 2.9198e-03, 5.1805e-03,
				 5.8735e-03, 4.4048e-01, 2.9015e-03,
				2.3313,
				0.63652,
				0.093821,
				 2.6221e-03, 4.3165e-03, 9.4362e-03,
				 2.6019e-03, 1.4914e-01, 4.2619e-03,
				0.87369,
				1.169,
				0.525405,
				 6.0284e-03, 1.5487e-04, 6.3771e-03,
				 5.9551e-03, 4.3306e-01, 1.5482e-04,
				1.03761,
				0.96502,
				0.530345,
				 7.1122e-03, 9.8406e-03, 9.5769e-03,
				 7.0123e-03, 4.4078e-01, 9.6505e-03,
				0.80108,
				1.33035,
				0.861141,
				 2.2175e-03, 1.6717e-03, 2.4071e-03,
				 2.2070e-03, 4.5729e-01, 1.6657e-03,
				0.78214,
				1.38609,
				0.484854,
				 1.1742e-03, 1.0622e-03, 6.7612e-03,
				 1.1712e-03, 4.0569e-01, 1.0597e-03,
				2.24858,
				0.64297,
				0.393457,
				 2.9668e-03, 3.7241e-03, 2.8906e-03,
				 2.9415e-03, 3.1059e-01, 3.6844e-03,
				1.66941,
				0.71378,
				0.671432,
				 3.1878e-03, 1.9812e-03, 6.7181e-03,
				 3.1637e-03, 4.0278e-01, 1.9719e-03,
				1.60972,
				0.72528,
				0.741259,
				 4.2417e-03, 4.8969e-03, 6.9514e-03,
				 4.2001e-03, 4.1504e-01, 4.8415e-03,
				0.79991,
				1.33358,
				0.520053,
				 5.0786e-03, 3.3949e-03, 6.7993e-04,
				 5.0241e-03, 4.1783e-01, 3.3705e-03,
				2.21606,
				0.64568,
				0.347714,
				 8.5516e-04, 9.5163e-03, 2.5479e-03,
				 8.5307e-04, 3.0202e-01, 9.2640e-03,
				1.75411,
				0.69934,
				0.149998,
				 2.6248e-03, 9.2033e-03, 2.2404e-03,
				 2.6080e-03, 2.1342e-01, 9.0001e-03,
				1.2119,
				0.85117,
				0.586093,
				 8.0101e-03, 5.2677e-04, 6.6783e-03,
				 7.8799e-03, 4.4213e-01, 5.2620e-04,
				1.5365,
				0.7412,
				0.262146,
				 2.9220e-04, 7.3786e-03, 8.4439e-03,
				 2.9201e-04, 3.0854e-01, 7.2566e-03,
				1.31362,
				0.80727,
				0.044455,
				 9.2885e-03, 2.6912e-03, 3.4446e-03,
				 9.1088e-03, 8.3862e-02, 2.6759e-03,
				0.66193,
				2.04384,
				0.754934,
				 7.3033e-03, 4.2284e-03, 7.8052e-03,
				 7.1618e-03, 3.6349e-01, 4.1805e-03,
				0.98983,
				1.01038,
				0.242786,
				 4.8861e-03, 5.4787e-03, 6.7533e-03,
				 4.8388e-03, 3.1318e-01, 5.4193e-03,
				0.75664,
				1.47414,
				0.442403,
				 5.7853e-03, 9.4274e-03, 6.7153e-05,
				 5.7115e-03, 3.8601e-01, 9.2332e-03
};

}
