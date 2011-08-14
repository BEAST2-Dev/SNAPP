package snap.app.mcmc;

import beast.app.BeastMCMC;
import beast.app.util.Arguments;
import beast.app.util.ErrorLogHandler;
import beast.app.util.MessageLogHandler;
import beast.app.util.Version;
import beast.util.Randomizer;
import beast.util.XMLParserException;
import jam.util.IconUtils;

import javax.swing.*;
import javax.swing.filechooser.FileFilter;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.logging.*;

public class SNAPPMCMC {

    private final static Version version = new SNAPPVersion();

    static class BeastConsoleApp extends jam.console.ConsoleApplication {

        public BeastConsoleApp(String nameString, String aboutString, javax.swing.Icon icon) throws IOException {
            super(nameString, aboutString, icon, false);
            getDefaultFrame().setDefaultCloseOperation(WindowConstants.DO_NOTHING_ON_CLOSE);
        }

        public void doStop() {
            // thread.stop is deprecated so need to send a message to running threads...
//            Iterator iter = parser.getThreads();
//            while (iter.hasNext()) {
//                Thread thread = (Thread) iter.next();
//                thread.stop(); // http://java.sun.com/j2se/1.5.0/docs/guide/misc/threadPrimitiveDeprecation.html
//            }
        }

        public void setTitle(String title) {
            getDefaultFrame().setTitle(title);
        }

        BeastMCMC beastMCMC;
    }

    public SNAPPMCMC(BeastMCMC beastMCMC, BeastConsoleApp consoleApp, int maxErrorCount) {

        final Logger infoLogger = Logger.getLogger("beast.app");
        try {

            if (consoleApp != null) {
                consoleApp.beastMCMC = beastMCMC;
            }

            // Add a handler to handle warnings and errors. This is a ConsoleHandler
            // so the messages will go to StdOut..
            Logger logger = Logger.getLogger("beast");

            Handler handler = new MessageLogHandler();
            handler.setFilter(new Filter() {
                public boolean isLoggable(LogRecord record) {
                    return record.getLevel().intValue() < Level.WARNING.intValue();
                }
            });
            logger.addHandler(handler);

//            // Add a handler to handle warnings and errors. This is a ConsoleHandler
//            // so the messages will go to StdErr..
//            handler = new ConsoleHandler();
//            handler.setFilter(new Filter() {
//                public boolean isLoggable(LogRecord record) {
//                    if (verbose) {
//                        return record.getLevel().intValue() >= Level.WARNING.intValue();
//                    } else {
//                        return record.getLevel().intValue() >= Level.SEVERE.intValue();
//                    }
//                }
//            });
//            logger.addHandler(handler);

            logger.setUseParentHandlers(false);

//            infoLogger.info("Parsing XML file: " + fileName);
//            infoLogger.info("  File encoding: " + fileReader.getEncoding());

            // This is a special logger that is for logging numerical and statistical errors
            // during the MCMC run. It will tolerate up to maxErrorCount before throwing a
            // RuntimeException to shut down the run.
            //Logger errorLogger = Logger.getLogger("error");
            handler = new ErrorLogHandler(maxErrorCount);
            handler.setLevel(Level.WARNING);
            logger.addHandler(handler);

            beastMCMC.run();

        } catch (java.io.IOException ioe) {
            infoLogger.severe("File error: " + ioe.getMessage());
        /* Catch exceptions and report useful information

        }  catch (org.xml.sax.SAXParseException spe) {
            if (spe.getMessage() != null && spe.getMessage().equals("Content is not allowed in prolog")) {
                infoLogger.severe("Parsing error - the input file, " + fileName + ", is not a valid XML file.");
            } else {
                infoLogger.severe("Error running file: " + fileName);
                infoLogger.severe("Parsing error - poorly formed XML (possibly not an XML file):\n" +
                        spe.getMessage());
            }
        } catch (org.w3c.dom.DOMException dome) {
            infoLogger.severe("Error running file: " + fileName);
            infoLogger.severe("Parsing error - poorly formed XML:\n" +
                    dome.getMessage());
        } catch (dr.xml.XMLParseException pxe) {
            if (pxe.getMessage() != null && pxe.getMessage().equals("Unknown root document element, beauti")) {
                infoLogger.severe("Error running file: " + fileName);
                infoLogger.severe(
                        "The file you just tried to run in BEAST is actually a BEAUti document.\n" +
                                "Although this uses XML, it is not a format that BEAST understands.\n" +
                                "These files are used by BEAUti to save and load your settings so that\n" +
                                "you can go back and alter them. To generate a BEAST file you must\n" +
                                "select the 'Generate BEAST File' option, either from the File menu or\n" +
                                "the button at the bottom right of the window.");

            } else {
                infoLogger.severe("Parsing error - poorly formed BEAST file, " + fileName + ":\n" +
                        pxe.getMessage());
            }

        } catch (RuntimeException rex) {
            if (rex.getMessage() != null && rex.getMessage().startsWith("The initial posterior is zero")) {
                infoLogger.warning("Error running file: " + fileName);
                infoLogger.severe(
                        "The initial model is invalid because state has a zero probability.\n\n" +
                                "If the log likelihood of the tree is -Inf, his may be because the\n" +
                                "initial, random tree is so large that it has an extremely bad\n" +
                                "likelihood which is being rounded to zero.\n\n" +
                                "Alternatively, it may be that the product of starting mutation rate\n" +
                                "and tree height is extremely small or extremely large. \n\n" +
                                "Finally, it may be that the initial state is incompatible with\n" +
                                "one or more 'hard' constraints (on monophyly or bounds on parameter\n" +
                                "values. This will result in Priors with zero probability.\n\n" +
                                "The individual components of the posterior are as follows:\n" +
                                rex.getMessage() + "\n" +
                                "For more information go to <http://beast.bio.ed.ac.uk/>.");
            } else {
                // This call never returns as another RuntimeException exception is raised by
                // the error log handler???
                infoLogger.warning("Error running file: " + fileName);
                System.err.println("Fatal exception: " + rex.getMessage());
                rex.printStackTrace(System.err);
            }
            */
        } catch (Exception ex) {
//            infoLogger.warning("Error running file: " + fileName);
            infoLogger.severe("Fatal exception: " + ex.getMessage());
            System.err.println("Fatal exception: " + ex.getMessage());
            ex.printStackTrace(System.err);
        }
    }

    static String getFileNameByDialog(String title) {
        JFileChooser fc = new JFileChooser(System.getProperty("user.dir"));
        fc.addChoosableFileFilter(new FileFilter() {
            public boolean accept(File f) {
                if (f.isDirectory()) {
                    return true;
                }
                String name = f.getName().toLowerCase();
                if (name.endsWith(".xml")) {
                    return true;
                }
                return false;
            }

            // The description of this filter
            public String getDescription() {
                return "xml files";
            }
        });

        fc.setDialogTitle(title);
        int rval = fc.showOpenDialog(null);

        if (rval == JFileChooser.APPROVE_OPTION) {
            return fc.getSelectedFile().toString();
        }
        return null;
    } // getFileNameByDialog


    public static void centreLine(String line, int pageWidth) {
        int n = pageWidth - line.length();
        int n1 = n / 2;
        for (int i = 0; i < n1; i++) {
            System.out.print(" ");
        }
        System.out.println(line);
    }

    public static void printTitle() {
        System.out.println();
        centreLine("SNAPP " + version.getVersionString() + ", " + version.getDateString(), 60);
        centreLine("SNP and AFLP Phylogenetic analysis", 60);
        for (String creditLine : version.getCredits()) {
            centreLine(creditLine, 60);
        }
        System.out.println();

    }

    public static void printUsage(Arguments arguments) {

        arguments.printUsage("snapp", "[<input-file-name>]");
        System.out.println();
        System.out.println("  Example: snapp test.xml");
        System.out.println("  Example: snapp -window test.xml");
        System.out.println("  Example: snapp -help");
        System.out.println();
        System.exit(0);
    }

    //Main method
    public static void main(String[] args) throws java.io.IOException {
    	List<String> MCMCargs = new ArrayList<String>();

        Arguments arguments = new Arguments(
                new Arguments.Option[]{

                        new Arguments.Option("window", "Provide a console window"),
                        new Arguments.Option("options", "Display an options dialog"),
                        new Arguments.Option("working", "Change working directory to input file's directory"),
                        new Arguments.LongOption("seed", "Specify a random number generator seed"),
                        new Arguments.StringOption("prefix", "PREFIX", "Specify a prefix for all output log filenames"),
                        new Arguments.Option("overwrite", "Allow overwriting of log files"),
                        new Arguments.Option("resume", "Allow appending of log files"),
                        // RRB: not sure what effect this option has
                        new Arguments.IntegerOption("errors", "Specify maximum number of numerical errors before stopping"),
                        new Arguments.IntegerOption("threads", "The number of computational threads to use (default auto)"),
                        new Arguments.Option("help", "Print this information and stop"),
                });

        try {
            arguments.parseArguments(args);
        } catch (Arguments.ArgumentException ae) {
            System.out.println();
            System.out.println(ae.getMessage());
            System.out.println();
            printUsage(arguments);
            System.exit(1);
        }

        if (arguments.hasOption("help")) {
            printUsage(arguments);
            System.exit(0);
        }

        final boolean window = arguments.hasOption("window");
        final boolean options = arguments.hasOption("options");
        final boolean working = arguments.hasOption("working");
        String fileNamePrefix = null;

        long seed = Randomizer.getSeed();

        int threadCount = 1;


        if (arguments.hasOption("prefix")) {
            fileNamePrefix = arguments.getStringOption("prefix");
        }


        if (arguments.hasOption("threads")) {
            threadCount = arguments.getIntegerOption("threads");
            if (threadCount < 0) {
                printTitle();
                System.err.println("The the number of threads should be >= 0");
                System.exit(1);
            }
        }

        if (arguments.hasOption("seed")) {
            seed = arguments.getLongOption("seed");
            if (seed <= 0) {
                printTitle();
                System.err.println("The random number seed should be > 0");
                System.exit(1);
            }
        }

        int maxErrorCount = 0;
        if (arguments.hasOption("errors")) {
            maxErrorCount = arguments.getIntegerOption("errors");
            if (maxErrorCount < 0) {
                maxErrorCount = 0;
            }
        }

        BeastConsoleApp consoleApp = null;

        String nameString = "SNAPP " + version.getVersionString();

        if (window) {
            System.setProperty("com.apple.macos.useScreenMenuBar", "true");
            System.setProperty("apple.laf.useScreenMenuBar", "true");
            System.setProperty("apple.awt.showGrowBox", "true");
            System.setProperty("beast.useWindow", "true");

            javax.swing.Icon icon = IconUtils.getIcon(SNAPPMCMC.class, "snapp.png");

            String aboutString = "<html><div style=\"font-family:sans-serif;\"><center>" +
                    "<div style=\"font-size:12;\"><p>Bayesian Evolutionary Analysis Sampling Trees<br>" +
                    "Version " + version.getVersionString() + ", " + version.getDateString() + "</p>" +
                    version.getHTMLCredits() +
                    "</div></center></div></html>";

            consoleApp = new BeastConsoleApp(nameString, aboutString, icon);
        }

        printTitle();

        File inputFile = null;

        if (options) {

            String titleString = "<html><center><p>SNAPP<br>" +
                    "Version " + version.getVersionString() + ", " + version.getDateString() + "</p></center></html>";
            javax.swing.Icon icon = IconUtils.getIcon(SNAPPMCMC.class, "snapp.png");

            SNAPPDialog dialog = new SNAPPDialog(new JFrame(), titleString, icon);

            if (!dialog.showDialog(nameString, seed)) {
               	System.exit(0);
            }

//            if (dialog.allowOverwrite()) {
//                allowOverwrite = true;
//            }
            switch (dialog.getLogginMode()) {
            case 0:/* do not ovewrite */ break;
            case 1: MCMCargs.add("-overwrite"); break;
            case 2: MCMCargs.add("-resume"); break;
            }

            seed = dialog.getSeed();
            threadCount = dialog.getThreadPoolSize();


            inputFile = dialog.getInputFile();

        } else {
        	if (arguments.hasOption("overwrite")) {
        		MCMCargs.add("-overwrite");	
        	} else if (arguments.hasOption("resume")) {
        		MCMCargs.add("-resume");	
        	}
        }

        if (inputFile == null) {

            String[] args2 = arguments.getLeftoverArguments();

            if (args2.length > 1) {
                System.err.println("Unknown option: " + args2[1]);
                System.err.println();
                printUsage(arguments);
                return;
            }

            String inputFileName = null;


            if (args2.length > 0) {
                inputFileName = args2[0];
                inputFile = new File(inputFileName);
            }

            if (inputFileName == null) {
                // No input file name was given so throw up a dialog box...
                inputFile = new File(getFileNameByDialog("SNAPP " + version.getVersionString() + " - Select XML input file"));
            }
        }

        if (inputFile != null && inputFile.getParent() != null && working) {
            System.setProperty("user.dir", inputFile.getParent());
        }

        if (window) {
            if (inputFile == null) {
                consoleApp.setTitle("null");
            } else {
                consoleApp.setTitle(inputFile.getName());
            }
        }

        if (fileNamePrefix != null && fileNamePrefix.trim().length() > 0) {
            System.setProperty("file.name.prefix", fileNamePrefix.trim());
        }

//        if (allowOverwrite) {
//            System.setProperty("log.allow.overwrite", "true");
//        }


        if (threadCount > 0) {
            System.setProperty("thread.count", String.valueOf(threadCount));
            MCMCargs.add("-threads");
            MCMCargs.add(threadCount + "");
        } else {
            MCMCargs.add("-threads");
            MCMCargs.add(Runtime.getRuntime().availableProcessors() + "");        	
        }

        MCMCargs.add("-seed");
        MCMCargs.add(seed + "");
        Randomizer.setSeed(seed);

        System.out.println();
        System.out.println("Random number seed: " + seed);
        System.out.println();

        // Construct the beast object
        BeastMCMC beastMCMC = new BeastMCMC();

        try {
            // set all the settings...
        	MCMCargs.add(inputFile.getAbsolutePath());
            beastMCMC.parseArgs(MCMCargs.toArray(new String[0]));
            
            new SNAPPMCMC(beastMCMC, consoleApp, maxErrorCount);
       } catch (RuntimeException rte) {
            if (window) {
                // This sleep for 2 seconds is to ensure that the final message
                // appears at the end of the console.
                try {
                    Thread.sleep(2000);
                } catch (InterruptedException e) {
                    e.printStackTrace();
                }
                System.out.println();
                System.out.println("SNAPP has terminated with an error. Please select QUIT from the menu.");
            }
            // logger.severe will throw a RTE but we want to keep the console visible
		} catch (XMLParserException e) {
			System.out.println(e.getMessage());
        } catch (Exception e) {
        	e.printStackTrace();
        }

        if (!window) {
            System.exit(0);
        }
    }
}


