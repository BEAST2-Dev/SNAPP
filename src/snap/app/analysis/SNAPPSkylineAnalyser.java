package snap.app.analysis;

import snap.util.SkylineAnalyser;
import beast.app.beauti.BeautiConfig;
import beast.app.beauti.BeautiDoc;
import beast.app.draw.BEASTObjectDialog;
import beast.app.draw.BEASTObjectPanel;
import beast.app.util.Application;
import beast.app.util.ConsoleApp;

public class SNAPPSkylineAnalyser {
		
	public static void main(final String[] args) throws Exception {
		Application main = null;
		try {
			// create the runnable class with application that we want to launch
			SkylineAnalyser analyser = new SkylineAnalyser();
			
			if (args.length == 0) {
				// try the GUI version

				// need to set the ID of the BEAST-object
				analyser.setID("SkylineAnalyser");
				
				// then initialise
				analyser.initAndValidate();
				
				// create BeautiDoc and beauti configuration
				BeautiDoc doc = new BeautiDoc();
				doc.beautiConfig = new BeautiConfig();
				doc.beautiConfig.initAndValidate();
			
				// create panel with entries for the application
				BEASTObjectPanel panel = new BEASTObjectPanel(analyser, analyser.getClass(), doc);
				
				// wrap panel in a dialog
				BEASTObjectDialog dialog = new BEASTObjectDialog(panel, null);
				if (dialog.showDialog()) {
					dialog.accept(analyser, doc);
					analyser.initAndValidate();

					// create a console to show standard error and standard output
					//analyser.consoleApp = 
					new ConsoleApp("SkylineAnalyser", // name 
							"Skyline Analyser" // console title
							);

					analyser.run();
				}
				return;
			}

			// continue with the command line version
			main = new Application(analyser);
			main.parseArgs(args, false);
			analyser.initAndValidate();
			analyser.run();
		} catch (Exception e) {
			System.out.println(e.getMessage());
			e.printStackTrace();
			if (main != null) {
				System.out.println(main.getUsage());
			}
		}
	}

}
