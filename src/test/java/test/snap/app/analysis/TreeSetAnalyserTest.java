package test.snap.app.analysis;

import java.io.File;

import org.junit.Test;

import snap.util.TreeSetAnalyser;
import snap.util.TreeSetAnalyser3;

import beast.base.inference.Logger;
import beast.base.inference.MCMC;
import beast.base.util.Randomizer;
import beast.base.parser.XMLParser;

import junit.framework.TestCase;

public class TreeSetAnalyserTest  extends TestCase {

	@Test
	// rudimentary test to see that the analyser does not crash
	public void testTreeSetAnalyser() throws Exception {
		
		// create tree log file
		try {
			Logger.FILE_MODE = Logger.LogFileMode.overwrite;
	
			int nSeed = 127;
			String sFileName = System.getProperty("user.dir") + "/examples2/test1.xml";
			Randomizer.setSeed(nSeed);
			System.out.println("Processing " + sFileName);
			XMLParser parser = new XMLParser();
			try {
				beast.base.inference.Runnable runable = parser.parseFile(new File(sFileName));
				if (runable instanceof MCMC) {
					MCMC mcmc = (MCMC) runable;
					mcmc.setInputValue("preBurnin", 0);
					mcmc.setInputValue("chainLength", 10000);
					mcmc.run();
				}
			} catch (Exception e) {
				System.out.println("ExampleXmlParsing::Failed for " + sFileName
						+ ": " + e.getMessage());
			}
			System.out.println("Done " + sFileName);
		} catch (Exception e) {
			System.out.println("exception thrown ");
			System.out.println(e.getMessage());;
		}
		
		// run TreeSetAnalyser
		try {
			String [] args = new String[]{"-tree", "((A,B),(C,D))", System.getProperty("user.dir") +"/test.127.trees"};
			TreeSetAnalyser3 analyser = new TreeSetAnalyser3();
        	if (analyser.parseArgs(args)) {
        		analyser.run();
        		// new TreeSetAnalyser(analyser, null, 1);
        	}
		} catch (Exception e) {
			e.printStackTrace();
			throw e;
		}
	}
	
	
}
