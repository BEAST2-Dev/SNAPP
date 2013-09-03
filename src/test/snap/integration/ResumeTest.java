package test.snap.integration;

import java.io.File;

import junit.framework.TestCase;

import org.junit.Test;

import beast.core.Logger;
import beast.core.MCMC;
import beast.util.Randomizer;
import beast.util.XMLParser;



/** check that a chain can be resumed after termination **/
public class ResumeTest  extends TestCase {
	
	final static String XML_FILE = "test1.xml";
	
	@Test
	public void test_ThatXmlExamplesRun() throws Exception {
		Randomizer.setSeed(127);
		Logger.FILE_MODE = Logger.LogFileMode.overwrite;
		String sDir = System.getProperty("user.dir") + "/examples";
		String sFileName = sDir + "/" + XML_FILE;

		System.out.println("Processing " + sFileName);
		XMLParser parser = new XMLParser();
		beast.core.Runnable runable = parser.parseFile(new File(sFileName));
		runable.setStateFile("tmp.state", false);
		if (runable instanceof MCMC) {
			MCMC mcmc = (MCMC) runable;
			mcmc.setInputValue("preBurnin", 0);
			mcmc.setInputValue("chainLength", 1000);
			mcmc.run();
		}
		System.out.println("Done " + sFileName);

		System.out.println("Resuming " + sFileName);
		Logger.FILE_MODE = Logger.LogFileMode.resume;
		parser = new XMLParser();
		runable = parser.parseFile(new File(sFileName));
		runable.setStateFile("tmp.state", true);
		if (runable instanceof MCMC) {
			MCMC mcmc = (MCMC) runable;
			mcmc.setInputValue("preBurnin", 0);
			mcmc.setInputValue("chainLength", 1000);
			mcmc.run();
		}
		System.out.println("Done " + sFileName);
	} // test_ThatXmlExamplesRun

} // class ResumeTest
