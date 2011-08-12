                    SNAPP v1.0.a 2011
                 SNAPP 2 development team 2011

Last updated: 12 Aug 2011

Contents:
1) INTRODUCTION
2) INSTALLING SNAPP
3) CONVERTING SEQUENCES
4) RUNNING SNAPP
5) ANALYZING RESULTS
6) SUPPORT & LINKS
7) ACKNOWLEDGMENTS 

___________________________________________________________________________
1) INTRODUCTION

SNAPP (SNP and AFLP Package for Phylogenetic analysis) is package for
evolutionary inference from molecular sequences. 

SNAPP uses a complex and powerful input format (specified in XML) to
describe the evolutionary model. This has advantages in terms of
flexibility in that the developers of SNAPP do not have to try and predict
every analysis that researchers may wish to perform and explicitly provide
an option for doing it. However, this flexibility means it is possible to
construct models that don't perform well under the Markov chain Monte Carlo
(MCMC) inference framework used. We cannot test every possible model that
can be used in SNAPP. There are two solutions to this: Firstly, we  supply
a range of recipes for commonly performed analyses that we know should work
in SNAPP and provide example input files for these (although, the actual
data can also produce unexpected behavour). Secondly, we provide advice and
tools for the diagnosis of problems and suggestions on how to fix them:

<http://snapp?????????????.cs.auckland.ac.nz/>

SNAPP is not a black-box into which you can put your data and expect an
easily interpretable answer. It requires careful inspection of the output
to check that it has performed correctly and usually will need tweaking,
adjustment and a number of runs to get a valid answer. Sorry.

SNAPP is build on top of Beast 2. More information on Beast 2 can be obtained
at <http://beast2.cs.auckland.ac.nz/>.
___________________________________________________________________________
2) INSTALLING SNAPP

SNAPP requires a Java Virtual Machine to run. Many systems will already
have this installed. It requires at least version 1.6 of Java to run. The
latest versions of Java can be downloaded from:

<http://java.sun.com/>

If in doubt type "java -version" to see what version of java is installed
(or if it is installed at all).

Mac OS X will already have a suitable version of Java installed.

Within the SNAPP package will be the following directories:
Directory       Contents
doc/            Documentation of SNAPP
examples/       Some NEXUS and XML files
lib/            Java & native libraries used by SNAPP 
bin/            Scripts of the corresponding OS
templates/      Templates to initiate BEAUti

___________________________________________________________________________
3) CONVERTING SEQUENCES

A program called "BEAUti" will import data in NEXUS format, allow you to
select various models and options and generate an XML file ready for use in
SNAPP.

To run BEAUti simply double-click the "BEAUti.exe" file in the SNAPP
folder. If this doesn't work then you may not have Java installed correctly. 
Try opening an MS-DOS window and typing:

	java -cp lib/SNAPP.jar beast.app.beauti.Beauti

__________________________________________________________________________
4) RUNNING SNAPP

To run SNAPP simply double-click the "SNAPP.exe" file in the SNAPP
folder. You will be asked to select a SNAPP XML input file.

Alternatively open a Command window and type:
	
	java -jar lib/SNAPP.jar input.xml

Where "input.xml" is the name of a SNAPP XML format file. This file can
either be created from scratch using a text editor or be created by the
BEAUti program from a NEXUS format file. 

For documentation on creating and tuning the input files look at the
documentation and tutorials on-line at:

Help -      <http://SNAPP2.cs.auckland.ac.nz/>
FAQ -       <http://SNAPP2.cs.auckland.ac.nz/index.php/FAQ>
Tutorials - <http://SNAPP2.cs.auckland.ac.nz/index.php/Main_Page#SNAPP_2_Tutorials>

  Usage: snapp [-window] [-options] [-working] [-seed] [-prefix <PREFIX>] [-overwrite] [-resume] [-errors <i>] [-threads <i>] [-help] [<input-file-name>]
    -window Provide a console window
    -options Display an options dialog
    -working Change working directory to input file's directory
    -seed Specify a random number generator seed
    -prefix Specify a prefix for all output log filenames
    -overwrite Allow overwriting of log files
    -resume Allow appending of log files
    -errors Specify maximum number of numerical errors before stopping
    -threads The number of computational threads to use (default auto)
    -help Print this information and stop

  Example: snapp test.xml
  Example: snapp -window test.xml
  Example: snapp -help
     
For example:

     java -jar lib/SNAPP.jar -seed 123456 -overwrite input.xml

___________________________________________________________________________
5) ANALYZING RESULTS

We have produced a powerful graphical program for analysing MCMC log files
(it can also analyse output from MrBayes and other MCMCs). This is called
'Tracer' and is available from the Tracer web site:

<http://tree.bio.ed.ac.uk/software/tracer>

Additionally, DensiTree is provided to analyse tree sets.

Alternatively, LogCombiner & TreeAnnotator distributed with Beast can be used. 
LogCombiner can combine log or tree files from multiple runs of SNAPP into a 
single combined results file (after removing appropriate burn-ins). TreeAnnotator 
can summarize a sample of trees from SNAPP using a single target tree, annotating 
it with posterior probabilities, HPD node heights and rates. This tree can then be
viewed in a new program called 'FigTree' which is available from:

<http://tree.bio.ed.ac.uk/software/figtree>

___________________________________________________________________________
6) SUPPORT & LINKS

SNAPP is an extremely complex program and as such will inevitably have
bugs. Please email us to discuss any problems:

<david.bryant@otago.ac.nz>
<remco@cs.auckland.ac.nz>

The SNAPP users' mailing-list is coming soon.

The website for SNAPP is here:

<http://SNAPP2.cs.auckland.ac.nz/>

Source code distributed under the GNU Lesser General Public License:

<http://code.google.com/p/snap-mcmc/>

___________________________________________________________________________
7) ACKNOWLEDGMENTS

Thanks for supplying code or assisting with the creation
or testing of SNAPP 2 development team.



