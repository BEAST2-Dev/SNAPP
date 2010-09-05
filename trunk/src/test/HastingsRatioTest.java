package test;

import java.util.Random;

/** clean room test for Hastings ratio for tree scaler **/
public class HastingsRatioTest {
	
	class Node {
		Node(double fHeight) {m_fHeight = fHeight;}
		Node m_left;
		Node m_right;
		double m_fHeight;
		Node copy() {
			Node copy = new Node(m_fHeight);
			if (m_left != null) {
				copy.m_left = m_left.copy();
			}
			if (m_right != null) {
				copy.m_right = m_right.copy();
			}
			return copy;
		}
	}
	
	void scale(Node node, double fScale) {
		node.m_fHeight *= fScale;
		if (node.m_left != null) {
			scale(node.m_left, fScale);
		}
		if (node.m_right != null) {
			scale(node.m_right, fScale);
		}
	}
	
	double sumOfHeights(Node node) {
		double fSumOfHeights = node.m_fHeight;
		if (node.m_left != null) {
			fSumOfHeights += sumOfHeights(node.m_left);
		}
		if (node.m_right != null) {
			fSumOfHeights += sumOfHeights(node.m_right);
		}
		return fSumOfHeights;
	}
	
	double logYuleOf(Node tree, double fLambda) {
		double fSumOfHeights = sumOfHeights(tree);
		return -fLambda * fSumOfHeights; // ignore constant
	}
	
	
	void doMCMC(double fLambda, double fScaleFactor, Node tree) {
		Random rand = new Random();
		double fLogProb = logYuleOf(tree, fLambda);
    	System.out.println("sample\tposterior\tsum.heights");
		for (int i = 0; i < 10000000; i++) {
			Node proposedTree = tree.copy();
	        //double fScale = (fScaleFactor + (rand.nextDouble() * ((1.0 / fScaleFactor) - fScaleFactor)));
	        double fScale = Math.exp((2*rand.nextDouble()-1.0)*fScaleFactor);
	        scale(proposedTree, fScale);
	        
	        
	        
	        // HERE THE HASTINGS RATIO IS CALCULATED!!!
	        
	        //double fLogHastingsRatio = -Math.log(fScale) * 1;
	        double fLogHastingsRatio = Math.log(fScale);
	        
       
	        
	        double fNewLogProb =  logYuleOf(proposedTree, fLambda);
            if (fNewLogProb -fLogProb + fLogHastingsRatio >=0 || 
            		rand.nextDouble() < Math.exp(fNewLogProb -fLogProb + fLogHastingsRatio)) {
            	// accept
            	tree = proposedTree;
            	fLogProb = fNewLogProb;
            }
	        if (i % 10000 == 0) {
	        	System.out.println(i + "\t" + fLogProb + "\t" + sumOfHeights(tree));
	        }
		}
	}
	

	void doTest() {
		Node tree = new Node(1.0);
		tree.m_left = new Node(0.5);
		tree.m_left.m_left = new Node(0.25); 
		tree.m_right = new Node(0.75);
		tree.m_right.m_left = new Node(0.5);
		tree.m_right.m_right = new Node(0.5);
		
		
		doMCMC(10, 0.75, tree);
		
	}
	
	public static void main(String [] args) {
		
		HastingsRatioTest test = new HastingsRatioTest();
		test.doTest();
	}
} // class HastingsRatioTest
