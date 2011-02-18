package snap.util;

import java.util.Vector;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.GammaDistributionImpl;
import org.apache.commons.math.distribution.NormalDistributionImpl;


import beast.util.Randomizer;

/**
 * reads a tree and generates coalescent rates according to CIR process
 *
 */
public class CIRSimulator {

	static double m_fAlpha = 200.0;
	static double m_fBeta = 2.0;
	static double m_fKappa = 138;

//	Generate a normal X  (mean 0, sd 1)
//	Generate Y with gamma(df-1,2)
//	return (X + sqrt(nc))^2 + Y
	static double generateFromNonCentralChiSquare(double df, double nc) throws MathException {
		NormalDistributionImpl normal = new NormalDistributionImpl(0, 1);
		double p = Randomizer.nextDouble();
		double X = normal.inverseCumulativeProbability(p);
		GammaDistributionImpl gamma = new GammaDistributionImpl(df-1.0, 2.0);
		p = Randomizer.nextDouble();
		double Y = gamma.inverseCumulativeProbability(p);
		double tmp = (X + Math.sqrt(nc));
		return tmp * tmp + Y;
	}
	
	static void generateThetasThroughCIR(Node node, double fParentRate) throws MathException {
		double fRate = 0.0;
		if (node.isRoot()) {
			// The stationary distribution (and the distribution of the root theta), is gamma with parameters \alpha, \beta.
			GammaDistributionImpl d = new GammaDistributionImpl(m_fAlpha, m_fBeta);
			double p = Randomizer.nextDouble();
			fRate = d.inverseCumulativeProbability(p);
			node.m_sMetaData = "theta=" +2.0/fRate +"";
		} else {
    		double t = node.m_fLength;
    		double r0 = fParentRate;
    		double df = 2 * m_fAlpha;
    		double nc = 2*(2.0/r0) * m_fBeta * Math.exp(-m_fKappa * t) / (1 - Math.exp(-m_fKappa*t));
			fRate = generateFromNonCentralChiSquare(df, nc);
			node.m_sMetaData = "theta=" +2.0/fRate +"";
		}
		
		if (!node.isLeaf()) {
			generateThetasThroughCIR(node.m_left, fRate);
			generateThetasThroughCIR(node.m_right, fRate);
		}
	}

	public static void main(String[] args) {
		try {
			Vector<String> sLabels = new Vector<String>();
			TreeFileParser parser = new TreeFileParser(sLabels, null, null, 0);
			Node tree = parser.parseNewick("(0:0.005,(1:0.003,(2:0.002,3:0.002):0.001):0.002):0.0;");
			//Node tree = parser.parseNewick("((0:0.004,(1:0.003,(2:0.002,3:0.002):0.001):0.001):0.001,(4:0.003,5:0.003):0.002):0.0;");
			generateThetasThroughCIR(tree, -1);
			System.out.println(tree.toNewick());
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
