package snap.app.mcmc;

import beast.app.util.Version;

/**
 * This class provides a mechanism for returning the version number of the
 * dr software. It relies on the administrator of the dr source using the
 * module tagging system in CVS. The method getVersionString() will return
 * the version of dr under the following condition: <BR>
 * 1. the dr source has been checked out *by tag* before being packaged for
 * distribution.
 * <p/>
 * Version last changed 2009/08/1 by AER
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 */
public class SNAPPVersion extends Version {

    /**
     * Version string: assumed to be in format x.x.x
     */
    private static final String VERSION = "1.0.a";

    private static final String DATE_STRING = "2009-2011";

    private static final boolean IS_PRERELEASE = true;

    private static final String REVISION = "$Rev$";

    public String getVersion() {
        return VERSION;
    }

    public String getVersionString() {
        return "v" + VERSION + (IS_PRERELEASE ? " Prerelease " + getBuildString() : "");
    }

    public String getDateString() {
        return DATE_STRING;
    }

    public String[] getCredits() {
        return new String[]{
                "Designed and developed by",
                "David Bryant, Remco Bouckaert",
                "",
                "Dept. Mathematics and Statistics",
                "University of Otago",
                "david.bryant@otago.ac.nz",
                "",
                "Department of Computer Science",
                "University of Auckland",
                "remco@cs.auckland.ac.nz",
                "",
                "Downloads, Help & Resources:",
                "\thttp://snapp.cs.auckland.ac.nz",
                "",
                "Source code distributed under the GNU Lesser General Public License:",
                "\thttp://code.google.com/p/snap-mcmc",
                ""
        };    
    }
    
    public String getBuildString() {
        return "r" + REVISION.split(" ")[1];
    }
}
