package test.snap.app.analysis;

import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.junit.Test;

import snap.app.analysis.TreeSetAnalyser;

import beast.core.Logger;
import beast.core.MCMC;
import beast.util.Randomizer;
import beast.util.XMLParser;

import junit.framework.TestCase;

public class TreeSetAnalyserTest  extends TestCase {

	@Test
	// rudimentary test to see that the analyser does not crash
	public void testTreeSetAnalyser() throws IOException {
		
		// create tree log file
		try {
			Logger.FILE_MODE = Logger.LogFileMode.overwrite;
	
			int nSeed = 127;
			String sFileName = System.getProperty("user.dir") + "/examples2/test1.xml";
			Randomizer.setSeed(nSeed);
			System.out.println("Processing " + sFileName);
			XMLParser parser = new XMLParser();
			try {
				beast.core.Runnable runable = parser.parseFile(new File(sFileName));
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
		TreeSetAnalyser.main(new String[]{"-file", "test.127.trees", "-tree", "((A,B),(C,D))"});
	}
	
	
}
