package nexus;

import jloda.util.Alert;
import jloda.util.Basic;
import jloda.util.parse.NexusStreamParser;

//import splits.core.Document;
import core.TaxaSet;

import java.io.*;
import java.net.URLClassLoader;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;


/**
 * The nexus assumptions block. This is where we keep track of the
 * assumptions under which the data is processed.
 */
public class Assumptions extends NexusBlock {
    /**
     * Identification string
     */
    public final static String NAME = "snapp_Assumptions";
    private TaxaSet extaxa;
    private List useTaxaSets;
    private List exchar;
    private List useCharSets;

    public final static String USE_ALL = "all"; // use all taxa or characters

    private boolean uptodate;



    private boolean excludeGaps;
    private boolean excludeMissing;
    private boolean excludeNonParsimony;
    private boolean excludeCodon1;
    private boolean excludeCodon2;
    private boolean excludeCodon3;
    private int excludeConstant;



    /**
     * Constructor
     */
    public Assumptions() {
        extaxa = null;
        useTaxaSets = null;
        exchar = null;
        useCharSets = null;


        uptodate = false;


        excludeGaps = false;
        excludeMissing = false;
        excludeNonParsimony = false;
        excludeCodon1 = false;
        excludeCodon2 = false;
        excludeCodon3 = false;
        excludeConstant = 0;
    }

    /**
     * has data been marked uptodate i.e. in a complete input file?
     *
     * @return true, if data doesn't need immediate updating
     */
    public boolean isUptodate() {
        return uptodate;
    }

    /**
     * data is update-to-date and next call to update will be ignored
     *
     * @param uptodate
     */
    public void setUptodate(boolean uptodate) {
        this.uptodate = uptodate;
    }



    /**
     * Gets the excludeGaps
     *
     * @return the excludeGaps
     */
    public boolean getExcludeGaps() {
        return excludeGaps;
    }

    /**
     * Sets the excludeGaps
     *
     * @param excGaps is excludeGaps
     */
    public void setExcludeGaps(boolean excGaps) {
        this.excludeGaps = excGaps;
    }

    /**
     * Gets the excludeMissing
     *
     * @return the excludeMissing
     */
    public boolean getExcludeMissing() {
        return excludeMissing;
    }

    /**
     * Sets the excludeMissing
     *
     * @param excMissing is excludeMissing
     */
    public void setExcludeMissing(boolean excMissing) {
        this.excludeMissing = excMissing;
    }

    /**
     * Gets the excludeNonParsimony
     *
     * @return the excludeNonParsimony
     */
    public boolean getExcludeNonParsimony() {
        return excludeNonParsimony;
    }

    /**
     * Sets the excludeNonParsimony
     *
     * @param excNonParsi is excludeNonParsimony
     */
    public void setExcludeNonParsimony(boolean excNonParsi) {
        this.excludeNonParsimony = excNonParsi;
    }

    /**
     * Gets the excludeCodon1
     *
     * @return the excludeCodon1
     */
    public boolean getExcludeCodon1() {
        return excludeCodon1;
    }

    /**
     * Sets the excludeCodon1
     *
     * @param excCodon1 is excludeCodon1
     */
    public void setExcludeCodon1(boolean excCodon1) {
        this.excludeCodon1 = excCodon1;
    }

    /**
     * Gets the excludeCodon2
     *
     * @return the excludeCodon2
     */
    public boolean getExcludeCodon2() {
        return excludeCodon2;
    }

    /**
     * Sets the excludeCodon2
     *
     * @param excCodon2 is excludeCodon2
     */
    public void setExcludeCodon2(boolean excCodon2) {
        this.excludeCodon2 = excCodon2;
    }

    /**
     * Gets the excludeCodon3
     *
     * @return the excludeCodon3
     */
    public boolean getExcludeCodon3() {
        return excludeCodon3;
    }

    /**
     * Sets the excludeCodon3
     *
     * @param excCodon3 is excludeCodon3
     */
    public void setExcludeCodon3(boolean excCodon3) {
        this.excludeCodon3 = excCodon3;
    }

    /**
     * Gets the excludeConstant
     *
     * @return the excludeConstant
     */
    public int getExcludeConstant() {
        return excludeConstant;
    }

    /**
     * Sets the excludeConstant
     *
     * @param excConstant is excludeConstant
     */
    public void setExcludeConstant(int excConstant) {
        this.excludeConstant = excConstant;
    }

    /**
     * Gets the set of taxa that are excluded from all computations
     *
     * @return the list of excluded taxa
     */
    public TaxaSet getExTaxa() {
        return extaxa;
    }

    /**
     * Gets the set of character positions that are excluded from all
     * computations
     *
     * @return the list of excluded characters
     */
    public List getExChar() {
        return exchar;
    }

    /**
     * Sets the set of taxa that are excluded from all computations
     *
     * @param extaxa the set of excluded taxa
     */
    public void setExTaxa(TaxaSet extaxa) {
        this.extaxa = extaxa;
    }

    public List getUseTaxaSets() {
        return useTaxaSets;
    }

    public void setUseTaxaSets(List useTaxaSets) {
        if (useTaxaSets != null && useTaxaSets.contains(USE_ALL))
            this.useTaxaSets = null;
        else
            this.useTaxaSets = useTaxaSets;
    }

    /**
     * Sets the list of character positions that are excluded from all
     * computations
     *
     * @param exchar the list of excluded characters
     */
    public void setExChar(List exchar) {
        this.exchar = exchar;
    }


    public List getUseCharSets() {
        return useCharSets;
    }

    public void setUseCharSets(List useCharSets) {
        if (useCharSets != null && useCharSets.contains(USE_ALL))
            this.useCharSets = null;
        else
            this.useCharSets = useCharSets;
    }



    /**
     * Read the assumptions block.
     *
     * @param np the nexus parser
     */
    public void read(NexusStreamParser np, Taxa taxa) throws IOException {
        setUptodate(false);

        if (taxa.getMustDetectLabels() == true)
            throw new IOException("line " + np.lineno() +
                    ": Can't read ASSUMPTIONS block because no taxlabels given in TAXA block");

        np.matchBeginBlock(NAME);

        try {
            np.pushPunctuationCharacters(NexusStreamParser.ASSIGNMENT_PUNCTUATION);
            // catch any expections below and reset punctuation there

            while (!np.peekMatchIgnoreCase("END;")) {
                if (np.peekMatchIgnoreCase("CHARSET"))
                    new Alert("Found CHARSET in ASSUMPTIONS block,please put into a SET-block");
                if (np.peekMatchIgnoreCase("TAXSET"))
                    new Alert("Found TAXSET in ASSUMPTIONS block,please put into a SET-block");
                if (np.peekMatchIgnoreCase("CHARPARTITION"))
                    new Alert("Found CHARPARTITION in ASSUMPTIONS block,please put into a SET-block");
                else if (np.peekMatchIgnoreCase("extaxa")) {
                    np.matchIgnoreCase("extaxa");
                    extaxa = new TaxaSet();

                    if (np.peekMatchIgnoreCase("=none;")) // restore all
                    {
                        np.matchIgnoreCase("=none;");
                    } else // hide the listed ones
                    {
                        List tokens = np.getTokensRespectCase("=", ";");

                        Iterator it = tokens.listIterator();
                        while (it.hasNext()) {
                            boolean ok = false; // found taxon yet?
                            String label = (String) it.next();

                            int i = taxa.getOriginalTaxa().indexOf(label);
                            if (i != -1) // label found
                                ok = true;
                            else // label not found, perhaps find its id?
                            {
                                try {
                                    i = Integer.parseInt(label);
                                    if (i > 0 && i <= taxa.getOriginalTaxa().getNtax())
                                        ok = true;
                                } catch (Exception ex) {
                                }
                            }
                            if (ok)
                                extaxa.set(i);
                        }
                    }
                } else if (np.peekMatchIgnoreCase("exchar")) {
                    np.matchIgnoreCase("exchar");

                    exchar = np.getIntegerList("=", ";");

                } else if (np.peekMatchIgnoreCase("usetaxset")) {
                    np.matchIgnoreCase("usetaxset");
                    List sets = np.getTokensRespectCase("=", ";");
                    setUseTaxaSets(sets);
                } else if (np.peekMatchIgnoreCase("usecharset")) {
                    np.matchIgnoreCase("usecharset");
                    List sets = np.getTokensRespectCase("=", ";");
                    setUseCharSets(sets);
                }  else if (np.peekMatchIgnoreCase("exclude")) {
                    List tokens = np.getTokensLowerCase("exclude", ";");
                    excludeGaps = np.findIgnoreCase(tokens, "no gaps", false, excludeGaps);
                    excludeGaps = np.findIgnoreCase(tokens, "gaps", true, excludeGaps);
                    excludeNonParsimony = np.findIgnoreCase(tokens, "no nonparsimony", false,
                            excludeNonParsimony);
                    excludeNonParsimony = np.findIgnoreCase(tokens, "nonparsimony", true,
                            excludeNonParsimony);
                    excludeMissing = np.findIgnoreCase(tokens, "no missing", false,
                            excludeMissing);
                    excludeMissing = np.findIgnoreCase(tokens, "missing", true, excludeMissing);
                    if (np.findIgnoreCase(tokens, "no constant", false, true) == false)
                        excludeConstant = 0;
                    excludeConstant = (int) np.findIgnoreCase(tokens, "constant=", excludeConstant);
                    if (np.findIgnoreCase(tokens, "constant", true, false) == true)
                        excludeConstant = -1;
                    excludeCodon1 = np.findIgnoreCase(tokens, "no codon1", false, excludeCodon1);
                    excludeCodon1 = np.findIgnoreCase(tokens, "codon1", true, excludeCodon1);
                    excludeCodon2 = np.findIgnoreCase(tokens, "no codon2", false, excludeCodon2);
                    excludeCodon2 = np.findIgnoreCase(tokens, "codon2", true, excludeCodon2);
                    excludeCodon3 = np.findIgnoreCase(tokens, "no codon3", false, excludeCodon3);
                    excludeCodon3 = np.findIgnoreCase(tokens, "codon3", true, excludeCodon3);

                    if (tokens.size() > 0)
                        throw new IOException("line " + np.lineno() + ": `"
                                + tokens.get(0) + "' unexpected in EXCLUDE");
                }  else
                    throw new IOException("line " + np.lineno() + ": unexpected: "
                            + np.getWordRespectCase());
            }
        }catch (IOException ex) {
            np.popPunctuationCharacters(); // restore punctation
            throw ex;
        }
        np.matchEndBlock();
    }



    public void write(Writer w, Taxa taxa) throws IOException {
        w.write("\nBEGIN " + Assumptions.NAME + ";\n");
        if (isUptodate())
            w.write("\tuptodate;\n");
        if (extaxa != null && extaxa.cardinality() > 0) {
            w.write("\textaxa=");
            for (int t = extaxa.getBits().nextSetBit(1); t > 0; t = extaxa.getBits().nextSetBit(t + 1)) {
                try {
                    if (taxa.getOriginalTaxa() != null) {
                        String name = taxa.getOriginalTaxa().getLabel(t);
                        w.write(" '" + name + "'");
                    }
                } catch (Exception ex) {
                    Basic.caught(ex);
                }
            }
            w.write(";\n");
        }
        /*
        if (useTaxaSets != null && useTaxaSets.size() > 0) {
            w.write("\tusetaxset=");
            Iterator it = useTaxaSets.iterator();
            while (it.hasNext()) {
                String label = (String) it.next();
                w.write(" '" + label + "'");
            }
            w.write(";\n");
        }
        */
        if (exchar != null && exchar.size() > 0) {
            w.write("\texchar=");
            Iterator it = exchar.listIterator();
            int first = 0, prev = 0;
            while (it.hasNext()) {
                int c = ((Integer) (it.next())).intValue();
                if (first == 0)
                    first = prev = c;
                else if (c == prev + 1)
                    prev = c;
                else // end of interval
                {
                    if (prev == first)
                        w.write(" " + first);
                    else
                        w.write(" " + first + "-" + prev);
                    first = prev = c;
                }
            }
            if (first > 0) {
                if (prev == first)
                    w.write(" " + first);
                else
                    w.write(" " + first + "-" + prev);
            }
            w.write(";\n");
        }
        if (useCharSets != null && useCharSets.size() > 0) {
            w.write("\tusecharset=");
            Iterator it = useCharSets.iterator();
            while (it.hasNext()) {
                String label = (String) it.next();
                w.write(" '" + label + "'");
            }
            w.write(";\n");
        }




        {
            StringWriter sw = new StringWriter();
            if (excludeGaps)
                sw.write(" gaps");
            if (excludeMissing)
                sw.write(" missing");
            if (excludeNonParsimony)
                sw.write(" nonparsimony");
            if (excludeConstant == -1)
                sw.write(" constant");
            else if (excludeConstant > 0)
                sw.write(" constant " + excludeConstant);
            if (excludeCodon1)
                sw.write(" codon1");
            if (excludeCodon2)
                sw.write(" codon2");
            if (excludeCodon3)
                sw.write(" codon3");
            if (sw.toString().length() > 0)
                w.write("\texclude " + sw.toString() + ";\n");
        }

        w.write("END; [" + Assumptions.NAME + "]\n");
    }




    /**
     * Shows assumptions. Shows all assumptions, whether blocks are defined or not
     *
     * @return object in nexus format
     */
    public String toString(Taxa taxa) {
        StringWriter sw = new StringWriter();
        try {
            write(sw, taxa);
        } catch (java.io.IOException ex) {
            return "";
        }
        return sw.toString();
    }


    /**
     * Show the usage of this block
     *
     * @param ps the PrintStream
     */
    public static void showUsage(PrintStream ps) {
        ps.println("BEGIN ST_ASSUMPTIONS;");
        ps.println("\t[EXTAXA={NONE|list-of-original-taxa-labels};]");
        // ps.println("\t[USETAXSET={ALL|list-of-taxset-labels};]");
        ps.println("\t[EXCHAR={NONE|list-of-original-char-positions};]");
        ps.println("\t[EXCLUDE [[NO] GAPS] [[NO] NONPARSIMONY]");
        ps.println("\t\t[{NO CONSTANT|CONSTANT [number]}]");
        ps.println("\t\t[[NO] CODON1] [[NO] CODON2] [[NO] CODON3];]");
        ps.println("\t[USECHARSET={ALL|list-of-charset-labels};]");
        ps.println("END;");
    }




    /**
     * clone the assumptions block
     *
     * @param taxa
     * @return a clone
     */
    public Assumptions clone(Taxa taxa) {
        Assumptions assumptions = new Assumptions();
        StringWriter w = new StringWriter();
        try {
            write(w, taxa);
            assumptions.read(new NexusStreamParser(new StringReader(w.toString())), taxa);
        } catch (Exception ex) {
            Basic.caught(ex);
        }
        return assumptions;

    }


}

//EOF
