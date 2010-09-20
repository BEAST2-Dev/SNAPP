/**
 * @version $Id: Trees.java,v 1.62 2008-07-10 15:01:00 bryant Exp $
 *
 * @author Daniel Huson
 *
 * This is a greatly reduced tree class, designed only to read and write tree
 * files, and the trees are stored as strings, rather than parsed.
 *
 */

package nexus;

import jloda.util.Basic;
import jloda.util.NotOwnerException;
import jloda.util.parse.NexusStreamParser;
import core.GenericException;

import java.io.*;
import java.util.*;

/**
 * NexusBlock trees class
 */
public class Trees extends NexusBlock {
    /**
     * Identification string
     */
    final public static String NAME = "Trees";
    private int ntrees = 0; // number of trees
    final private Vector<String> names = new Vector<String>(); // list of names
    final private Vector<String> trees = new Vector<String>(); // list of tree strings
    final private Vector<Boolean> isRooted = new Vector<Boolean>(); //Flag indicating if a tree is rooted.
    final private Map<String,String> translate = new HashMap<String,String>(); // maps node labels to taxon labels

    /**
     * Construct a new Trees object.
     */
    public Trees() {
        super();
    }

    /**
     * Constructs a new Trees object and adds the given phylogenetic tree to it
     *
     * @param name the name of the tree
     * @param tree the phylogenetic tree
     * @param taxa the taxa
     */
    public Trees(String name, boolean treeIsRooted, String tree) {
        super();
        try {
            this.addTree(name, treeIsRooted, tree);
        } catch (Exception ex) {
            Basic.caught(ex);
        }
    }

    /**
     * clears all the data associated with this trees block
     */
    public void clear() {
        translate.clear();
        trees.clear();
        names.clear();
        isRooted.clear();
        ntrees = 0;
    }



    /**
     * Get the number of trees
     *
     * @return number of trees
     */
    public int getNtrees() {
        return ntrees;
    }

    /**
     * Returns the i-th tree name.
     * Trees are numbered 1 to ntrees
     *
     * @param i the number of the tree
     * @return the i-th tree name
     */
    public String getName(int i) {
        return  names.elementAt(i - 1);
    }

    /**
     * sets the i-th tree name
     *
     * @param i
     * @param name
     */
    public void setName(int i, String name) {
        names.setElementAt(name, i - 1);
    }


    /**
     * returns the index of the named tree.
     * Trees are numbered 1 to ntrees
     *
     * @param name
     * @return index of named tree
     */
    public int indexOf(String name) {
        return names.indexOf(name) + 1;
    }

    /**
     * Returns the i-th tree
     *
     * @param i the number of the tree
     * @return the i-th tree
     */
    public String getTree(int i) {
        return trees.elementAt(i - 1);
    }

    /**
     * Returns the nexus flag [&R] indicating whether the tree should be considered
     * as rooted
     *
     * @param i
     * @return String  Returns [&R] if rooted, and "" otherwise.
     */
    public String getFlags(int i) {
        if (isRooted.get(i))
            return "[&R]";
        else
            return "";
    }

    /**
     * remove the tree at index i.
     * The index must be between 1 and ntrees
     *
     * @param i
     */
    public void removeTree(int i) {

        if (names.remove(i - 1) != null && trees.remove(i - 1) != null && isRooted.remove(i-1))
            ntrees--;
    }

    /**
     * Adds a tree. Does no checking of the translate table, just assumes that this
     *  is a valid string for the taxa file.
     *
     * @param name the name of the tree
     * @param tree the phylogenetic tree
     * @param taxa the taxa block
     */
    public void addTree(String name, boolean treeIsRooted, String tree)
            throws GenericException, NotOwnerException {


        ntrees++;
        trees.add(tree);
        names.add(name); //TODO: Check that names are unique.
        isRooted.add(treeIsRooted);
    }


    /**
     * Writes trees taxa object in nexus format
     *
     * @param w a writer
     */
    public void write(Writer w, Taxa taxa) throws IOException {
        w.write("\nBEGIN " + Trees.NAME + ";\n");

        if (translate.size() > 0 &&
                (translate.size() > taxa.getNtax()  || !translateIsOneToOne(translate))) {
            w.write("TRANSLATE\n");
            Set keys = translate.keySet();

            Iterator it = keys.iterator();
            while (it.hasNext()) {
                String nodeLabel = (String) it.next();
                w.write("\t '" + nodeLabel + "'\t'" + translate.get(nodeLabel) + "',\n");
            }
            w.write(";\n");
        }

        for (int t = 1; t <= getNtrees(); t++) {
            w.write("[" + t + "] tree '" + getName(t) + "'=" + getFlags(t) + " " + getTree(t) + ";\n");
        }
        w.write("END; [" + Trees.NAME + "]\n");
    }

    /**
     * is translate one-to-one?
     *
     * @param translate
     * @return true, if one-to-one
     */
    public boolean translateIsOneToOne(Map translate) {
        Set keys = translate.keySet();

        Iterator it = keys.iterator();
        while (it.hasNext()) {
            String nodeLabel = (String) it.next();
            if (!nodeLabel.equals(translate.get(nodeLabel)))
                return false;
        }
        return true;

    }


    /**
     * Reads a tree object in NexusBlock format
     *
     * @param np   nexus stream parser
     * @param taxa the taxa block
     */
    public void read(NexusStreamParser np, Taxa taxa) throws GenericException, IOException {
        np.matchBeginBlock(NAME);
        clear();

        /*if (np.peekMatchIgnoreCase("properties")) {
            List tokens = np.getTokensLowerCase("properties", ";");
            if (np.findIgnoreCase(tokens, "no partialtrees"))
                partial = false;
            if (np.findIgnoreCase(tokens, "partialtrees=no"))
                partial = false;
            if (np.findIgnoreCase(tokens, "partialtrees=yes"))
                partial = true;
            if (np.findIgnoreCase(tokens, "partialtrees"))
                partial = true;
            if (np.findIgnoreCase(tokens, "rooted=yes")) {
                rooted = true;
                rootedGloballySet = true;
            }
            if (np.findIgnoreCase(tokens, "rooted=no")) {
                rooted = false;
                rootedGloballySet = true;
            }
            if (tokens.size() != 0)
                throw new IOException("line " + np.lineno() + ": `" + tokens + "' unexpected in PROPERTIES");
        }*/

        if (np.peekMatchIgnoreCase("translate")) {
            List<String> taxlabels = new ArrayList<String>();
            np.matchIgnoreCase("translate");
            while (!np.peekMatchIgnoreCase(";")) {
                String nodelabel = np.getWordRespectCase();
                String taxlabel = np.getWordRespectCase();
// if we have a translate and have to detect the Tasa use the taxlabels
                taxlabels.add(taxlabel);
                translate.put(nodelabel, taxlabel);

                if (!np.peekMatchIgnoreCase(";"))
                    np.matchIgnoreCase(",");
            }
            np.matchIgnoreCase(";");           
        } else if (taxa.getMustDetectLabels()) {
            throw new GenericException("line " + np.lineno() +
                    ": Taxon labels not given in taxa block, thus TRANSLATE-statement required");
        } else {
            // set the translation table from the taxa:
            translate.clear();
            for (int t = 1; t <= taxa.getNtax(); t++)
                translate.put(taxa.getLabel(t), taxa.getLabel(t));
        }
        while (np.peekMatchIgnoreCase("tree")) {
            np.matchIgnoreCase("tree");
            if (np.peekMatchRespectCase("*"))
                np.matchRespectCase("*"); // don't know why PAUP puts this star in the file....

            String name = np.getWordRespectCase();
            name = name.replaceAll("[ \t\b]+", "_");
            name = name.replaceAll("[:;,]+", ".");
            name = name.replaceAll("\\[", "(");
            name = name.replaceAll("\\]", ")");

            np.matchIgnoreCase("=");
            np.getComment(); // clears comments


            StringBuffer buf = new StringBuffer();

            LinkedList tokensToCome = (LinkedList) np.getTokensRespectCase(null, ";");
            for (Iterator p = tokensToCome.iterator(); p.hasNext();) {
                buf.append(p.next());
            }

            /*
            while (!np.peekMatchIgnoreCase(";"))
                buf.append(np.getWordRespectCase());
            np.matchIgnoreCase(";");
              */
            String comment = np.getComment();
            boolean isRooted = (comment != null && comment.equalsIgnoreCase("&R"));
            addTree(name,isRooted,buf.toString());
        }
        np.matchEndBlock();
    }



    /**
     * Returns true if tree i is rooted, else false
     *
     * @param i number of the tree
     * @return true if tree i is rooted
     */
    public boolean isRooted(int i) {
        return this.isRooted.get(i);
    }

    /**
     * Produces a string representation of a NexusBlock object
     *
     * @return object in nexus format
     */
    public String toString(Taxa taxa) {
        StringWriter sw = new StringWriter();
        try {
            write(sw, taxa);
        } catch (Exception ex) {
            return "()";
        }
        return sw.toString();
    }

    /**
     * show the usage of this block
     *
     * @param ps the print stream
     */
    public static void showUsage(PrintStream ps) {
        ps.println("BEGIN " + Trees.NAME);
       // ps.println("[PROPERTIES PARTIALTREES={YES|NO} ROOTED={YES|NO};]");
        ps.println("[TRANSLATE");
        ps.println("    nodeLabel1 taxon1,");
        ps.println("    nodeLabel2 taxon2,");
        ps.println("    ...");
        ps.println("    nodeLabelN taxonN");
        ps.println(";]");
        ps.println("[TREE name1 = tree1-in-Newick-format;]");
        ps.println("[TREE name2 = tree2-in-Newick-format;]");
        ps.println("...");
        ps.println("[TREE nameM = treeM-in-Newick-format;]");
        ps.println("END;");
    }





    /**
     * clones a trees object
     *
     * @param taxa
     * @return a clone
     */
    public Trees clone(Taxa taxa) {
        Trees trees = new Trees();
        try {
            StringWriter sw = new StringWriter();
            this.write(sw, taxa);
            StringReader sr = new StringReader(sw.toString());
            trees.read(new NexusStreamParser(sr), taxa);
        } catch (Exception ex) {
            Basic.caught(ex);
        }
        return trees;
    }


    /**
     * copies a trees object
     *
     * @param taxa
     * @param src  source tree
     */
    public void copy(Taxa taxa, Trees src) {
        try {
            StringWriter sw = new StringWriter();
            src.write(sw, taxa);
            StringReader sr = new StringReader(sw.toString());
            this.read(new NexusStreamParser(sr), taxa);
        } catch (Exception ex) {
            Basic.caught(ex);
        }
    }

    /**
     * sets the translate map to the identity mapping taxa->taxa
     *
     * @param taxa
     */
    public void setIdentityTranslate(Taxa taxa) {
        translate.clear();
        for (int t = 1; t <= taxa.getNtax(); t++)
            translate.put(taxa.getLabel(t), taxa.getLabel(t));
    }

    /**
     * Sets the translate map to taxaid->taxa
     *
     * @param taxa
     */
    public void setNumberedIdentityTranslate(Taxa taxa) {
        translate.clear();
        for (int t = 1; t <= taxa.getNtax(); t++)
            translate.put("" + t, taxa.getLabel(t));
    }
}

// EOF
