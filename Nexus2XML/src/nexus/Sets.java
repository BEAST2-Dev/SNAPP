/*
 * Created on Feb 6, 2005
 *
 * The NEXUS Sets block allows the user to define sets of characters
 * and taxa, as well as partitions of taxo or splits.
 *
 * The syntax for this block is taken directly from PAUP manuel, 4.0b6
 *TODO: Current version only implements taxa sets.
 */
package nexus;

import jloda.phylo.PhyloTree;
import jloda.util.Basic;
import jloda.util.Partition;
import jloda.util.parse.NexusStreamParser;
import core.TaxaSet;

import java.io.IOException;
import java.io.PrintStream;
import java.io.StringReader;
import java.io.Writer;
import java.util.*;


/**
 * Sets of taxa and chars
 *
 * @author bryant, taxonomy added by Daniel Huson, 12.05#
 */
public class Sets extends NexusBlock implements Cloneable {
    /**
     * Identification string
     */
    final public static String NAME = "Sets";
    private Map taxSets; // Assigns taxa sets to names of taxa sets
    private Map charSets; //Assign char set names to char sets.
    private Map charParts;

    /**
     * Construct a new Sets object
     */
    public Sets() {
        taxSets = new TreeMap();
        charSets = new TreeMap();
        charParts = new TreeMap();
    }

    /**
     * get number of taxa sets
     *
     * @return int number of taxa sets
     */
    public int getNumTaxSets() {
        return taxSets.size();
    }


    /**
     * get number of character sets
     *
     * @return int number of character sets
     */
    public int getNumCharSets() {
        return charSets.size();
    }

    /**
     * get number of character partitions
     *
     * @return int number of character partitions
     */
    public int getNumCharPartitions() {
        return charParts.size();
    }



    /**
     * Returns the taxa set associated with the given name
     *
     * @param name
     * @return set of names of taxa. Will return null if name does not appear
     *         or if name is present, but there is no set for that name.
     */
    public Set getTaxSet(String name) {
        return (Set) taxSets.get(name);
    }

    /**
     * Returns the taxa set associated with the given name
     *
     * @param name
     * @param taxa
     * @return taxa set. Will return null if name does not appear
     *         or if name is present, but there is no set for that name.
     */
    public TaxaSet getTaxSet(String name, Taxa taxa) {
        Set taxaSet = (Set) taxSets.get(name);
        if (taxaSet != null && taxaSet.size() > 0) {
            TaxaSet result = new TaxaSet();
            for (Iterator it = taxaSet.iterator(); it.hasNext();) {
                String tax = (String) it.next();
                int pos = taxa.indexOf(tax);
                if (pos > 0)
                    result.set(pos);
            }
            return result;
        }
        return null;
    }

    /**
     * gets the set of character set names
     *
     * @return character set names
     */
    public Set getCharSetNames() {
        return charSets.keySet();
    }

    /**
     * gets the set of character partition names
     *
     * @return character partition names
     */
    public Set getCharPartitionNames() {
        return charParts.keySet();
    }


    /**
     * gets the set of taxa set names
     *
     * @return taxa set names
     */
    public Set getTaxaSetNames() {
        return taxSets.keySet();
    }

   

    /**
     * gets the named char set
     *
     * @param charSetName
     * @return Set
     */
    public Set getCharSet(String charSetName) {
        return (Set) charSets.get(charSetName);
    }

    /**
     * gets the named char partition
     *
     * @param charPartitionName
     * @return Partition
     */
    public Partition getCharPartition(String charPartitionName) {
        return (Partition) charParts.get(charPartitionName);
    }

    /**
     * Clear all  sets and partitions
     */
    public void clear() {
        taxSets.clear();
        charSets.clear();
        charParts.clear();
    }

    /**
     * Clear all taxa sets
     */
    public void clearTaxaSets() {
        taxSets.clear();
    }


    /**
     * adds a new tax set
     *
     * @param name   Name of taxa set being added
     * @param taxSet Taxa set being added
     * @param taxa   Taxa block
     * @return true if a set with this name already exists (
     *         in which case it is replaced)
     */
    public boolean addTaxSet(String name, TaxaSet taxSet, Taxa taxa) {
        Set set = new HashSet();
        for (int t = 1; t <= taxa.getNtax(); t++)
            if (taxSet.get(t))
                set.add(taxa.getLabel(t));
        return addTaxSet(name, set);
    }

    /**
     * adds a new tax set
     *
     * @param name Name of taxa set being added
     * @param set  names of taxa
     * @return true if a set with this name already exists (
     *         in which case it is replaced)
     */
    public boolean addTaxSet(String name, Set set) {
        boolean alreadyThere = taxSets.containsKey(name);
        taxSets.put(name, set);
        return alreadyThere;
    }



    /**
     * add a character
     *
     * @param name
     * @param set
     * @return true, if name already used
     */
    public boolean addCharSet(String name, Set set) {
        boolean alreadyThere = charSets.containsKey(name);
        charSets.put(name, set);
        return alreadyThere;
    }

    /**
     * add a character
     *
     * @param name
     * @param first
     * @param last
     * @return true, if name already used
     */
    public boolean addCharSet(String name, int first, int last) {
        boolean alreadyThere;
        alreadyThere = charSets.containsKey(name);
        SortedSet set = new TreeSet();
        for (int i = first; i <= last; i++)
            set.add(new Integer(i));
        charSets.put(name, set);
        return alreadyThere;
    }

    /**
     * add character partition
     *
     * @param name
     * @param partition
     * @return true, if name already used
     */
    public boolean addCharPartition(String name, Partition partition) {
        boolean alreadyThere = charParts.containsKey(name);
        charParts.put(name, partition);
        return alreadyThere;
    }

    /**
     * removes taxa set
     *
     * @param name
     * @return true iff set with that name was present
     */
    public boolean removeTaxSet(String name) {
        boolean keyPresent = taxSets.containsKey(name);
        taxSets.remove(name);
        return keyPresent;
    }

   

    /**
     * removes char set
     *
     * @param charSetName
     * @return true iff set with that name was present
     */
    public boolean removeCharSet(String charSetName) {
        boolean keyPresent = charSets.containsKey(charSetName);
        charSets.remove(charSetName);
        return keyPresent;
    }

    /**
     * removes char partition
     *
     * @param charPartName
     * @return true iff set with that name was present
     */
    public boolean removeCharPart(String charPartName) {
        boolean keyPresent = charParts.containsKey(charPartName);
        charSets.remove(charPartName);
        return keyPresent;
    }


    /**
     * usage
     *
     * @param ps
     */
    public static void showUsage(PrintStream ps) {
        ps.println("BEGIN " + NAME + ";");
        ps.println("\t[TAXSET taxset-name    = taxon-list;]");
        ps.println("\t[TAXONOMY taxonomy-name = Newick-tree;]");
        ps.println("\t[CHARSET charset-name  = character-list;]");
        ps.println("\t[CHARPARTITION charpart-name = 1:charset-name|character-list [,...]; ");
        ps.println("\t...");
        ps.println("END;");
    }

    /**
     * print an int set to a string
     *
     * @param intSet
     * @return a string
     */
    private String toString(Set intSet) {
        int runstart = 0;
        int runlength = 0;
        String output = "";

        for (Iterator p = intSet.iterator(); p.hasNext();) {
            int x = ((Integer) p.next()).intValue();
            if (runlength == 0) { //Just started
                runstart = x;
                runlength = 1;
            } else if (x == runstart + runlength) { //extend current run
                runlength += 1;
            } else {
                //end of current run. print last run and reset.
                if (runlength <= 2) {
                    for (int y = runstart; y < runstart + runlength; y++) {
                        output += y;
                        output += ' ';
                    }
                } else {
                    //print the first, followed by a dash, followed by the last
                    output += runstart;
                    output += '-';
                    output += (runstart + runlength - 1);
                    output += ' ';
                }

                runstart = x;
                runlength = 1;
            }
        }
        //Now print out the remaining entries.
        if (runlength <= 2) {
            for (int y = runstart; y < runstart + runlength; y++) {
                output += y;
            }
        } else {
            //print the first, followed by a dash, followed by a period (denoting last element)
            output += runstart;
            output += '-';
            output += runstart + runlength - 1;
        }
        return output;
    }

    /**
     * return a int partition as a string
     *
     * @param partition
     * @return string
     */
    private String toString(Partition partition) {
        String s = "";
        for (int i = 1; i <= partition.getNumBlocks(); i++) {
            s += i;
            s += ":";
            s += toString(partition.getBlock(i));
            if (i < partition.getNumBlocks())
                s += ",";
        }
        return s;
    }

    /**
     * write a sets block
     *
     * @param w
     * @param taxa
     * @throws java.io.IOException
     */
    public void write(Writer w, Taxa taxa) throws java.io.IOException {
        w.write("\nBEGIN " + Sets.NAME + ";\n");
        for (Iterator p = taxSets.keySet().iterator(); p.hasNext();) {
            String name = (String) (p.next());
            w.write("  TAXSET '" + name + "' =");
            Set set = getTaxSet(name);
            for (Iterator it = set.iterator(); it.hasNext();)
                w.write(" '" + it.next() + "'");
            w.write(";\n");
        }


        for (Iterator p = charSets.keySet().iterator(); p.hasNext();) {
            String name = (String) (p.next());
            w.write("  CHARSET '" + name + "' = ");
            w.write(toString(getCharSet(name)));
            w.write(";\n");
        }
        for (Iterator p = charParts.keySet().iterator(); p.hasNext();) {
            String name = (String) (p.next());
            w.write("  CHARPARTITION '" + name + "' = ");
            w.write(toString(getCharPartition(name)));
            w.write(";\n");
        }
        w.write("END; [Sets]\n");
    }


    /**
     * Reads a string specifying a set of taxon.
     * The set is specified using taxon ids or taxon names. Ranges can be specified using a hyphen,
     * and a range ending with '.' is shorthand for a range ending with the last taxon.
     *
     * @param np
     * @param taxa
     * @return Set of sets read in.
     * @throws IOException
     */
    private Set readTaxSet(NexusStreamParser np, Taxa taxa) throws IOException {
        Set result = new HashSet();

        int previous = -1;
        boolean inRun = false;

        while (!np.peekMatchIgnoreCase(";")) {
            if (np.peekMatchIgnoreCase("-")) {
                if (inRun || previous == -1)
                    throw new IOException("line " + np.lineno() + ": '-' unexpected");
                np.matchIgnoreCase("-");
                inRun = true;
            } else {
                String next = np.getWordRespectCase();
                int current = -1;
                if (Basic.isInteger(next)) {
                    np.pushBack();
                    current = np.getInt(1, taxa.getNtax());
                }
                if (current == -1 && next.equals("."))
                    current = taxa.getNtax();
                if (current == -1) // not a number, maybe a name?
                {
                    result.add(next);
                    current = taxa.indexOf(next);
                }
                if (inRun) {
                    if (current < previous)
                        throw new IOException("line " + np.lineno() + ": illegal range: "
                                + previous + " - " + current);
                    for (int i = previous + 1; i <= current; i++)
                        result.add(taxa.getLabel(i));
                    inRun = false;
                    previous = -1;
                } else {
                    result.add(taxa.getLabel(current));
                    previous = current;
                }
            }
        }
        if (inRun)
            throw new IOException("line " + np.lineno() + ": '-' unexpected");
        np.matchIgnoreCase(";");
        return result;
    }

    /**
     * reads a characters set
     *
     * @param np
     * @param chars
     * @return Set of sets read in.
     * @throws IOException
     */

    private Set readCharSet(NexusStreamParser np, Characters chars) throws IOException {
        Set set = new TreeSet();

        int previous = -1;
        boolean inRun = false;

        while (!np.peekMatchIgnoreCase(";")) {
            if (np.peekMatchIgnoreCase("-")) {
                if (inRun || previous == -1)
                    throw new IOException("line " + np.lineno() + ": '-' unexpected");
                np.matchIgnoreCase("-");
                inRun = true;
            } else {
                int current;
                if (np.peekMatchRespectCase(".")) {
                    np.matchRespectCase(".");
                    current = chars.getNchar();
                } else
                    current = np.getInt(1, chars.getNchar());
                if (inRun) {
                    if (current < previous)
                        throw new IOException("line " + np.lineno() + ": illegal range: "
                                + previous + " - " + current);
                    for (int i = previous + 1; i <= current; i++)
                        set.add(new Integer(i));
                    inRun = false;
                    previous = -1;
                } else {
                    set.add(new Integer(current));
                    previous = current;
                }
            }
        }
        np.matchIgnoreCase(";");
        if (inRun)
            throw new IOException("line " + np.lineno() + ": '-' unexpected");

        //System.err.println("charset: " + set);
        return set;
    }

    /**
     * read a character parititon
     *
     * @param np
     * @param chars
     * @return Set of set read in.
     * @throws IOException
     */

    private Partition readCharPartition(NexusStreamParser np, Characters chars) throws IOException {
        Partition result = new Partition();

        int previous = -1;
        boolean inRun = false;

        boolean expectNextBlock = true;
        boolean expectComma = false;
        Set block = null;
        int blockNumber = 1;
        String blockName = "";

        while (!np.peekMatchIgnoreCase(";")) {
            if (expectComma) {
                if (!np.peekMatchIgnoreCase(",")) {
                    throw new IOException("line " + np.lineno() + ": ',' expected");
                }
                expectComma = false;
            }
            if (expectNextBlock) {
                if (inRun)
                    throw new IOException("line " + np.lineno() + ": '-' unexpected");

                //Blocks in the partition can be numbered (in which case they have
                //to correspond to block number) or names.
                if (np.peekMatchIgnoreCase("" + blockNumber)) {
                    np.matchIgnoreCase(blockNumber + " :");
                    blockName = Integer.toString(blockNumber);
                } else {
                    String next = np.getWordRespectCase();
                    if (Basic.isInteger(next))
                        throw new IOException("line " + np.lineno() + ": Incorrect partition block number: " + next);
                    np.matchIgnoreCase(":");
                    blockName = next;
                }

                blockNumber++;
                block = new TreeSet();
                expectNextBlock = false;
            } else {
                if (np.peekMatchIgnoreCase("-")) {
                    if (inRun || previous == -1)
                        throw new IOException("line " + np.lineno() + ": '-' unexpected");
                    np.matchIgnoreCase("-");
                    inRun = true;
                } else {
                    if (inRun) {
                        int current = np.getInt(1, chars.getNchar());
                        if (current < previous)
                            throw new IOException("line " + np.lineno() + ": illegal range: "
                                    + previous + " - " + current);
                        for (int i = previous + 1; i <= current; i++)
                            block.add(new Integer(i));
                        inRun = false;
                        previous = -1;
                    } else if (np.peekMatchIgnoreCase(",")) {
                        np.matchIgnoreCase(",");
                        try {
                            result.addBlock(block, blockName);
                        }
                        catch (Exception ex) {
                            throw new IOException(ex.toString());
                        }
                        expectNextBlock = true;
                    } else {
                        String next = np.getWordRespectCase();
                        if (Basic.isInteger(next)) {
                            np.pushBack();
                            int current = np.getInt(1, chars.getNchar());
                            block.add(new Integer(current));
                            previous = current;
                        } else // should be a defined character set
                        {
                            Set cSet = this.getCharSet(next);
                            if (cSet == null)
                                throw new IOException("line " + np.lineno() + ": Illegal character or character set: "
                                        + previous + " - " + next);
                            block.addAll(cSet);
                            previous = -1;
                            expectComma = true;
                        }
                    }
                }
            }
        }

        //Problem here. Last block is not followed by a comma. Have commented this out.
        // if (expectNextBlock == false)
        //    throw new IOException("line " + np.lineno() + ": expected ','");
        np.matchIgnoreCase(";");
        return result;
    }

    /**
     * reads a sets block
     *
     * @param np
     * @param taxa
     * @param chars
     * @throws IOException
     */
    public void read(NexusStreamParser np, Taxa taxa, Characters chars) throws IOException {
        if (taxa.getMustDetectLabels())
            throw new IOException("line " + np.lineno() +
                    ": Can't read SETS block because no taxlabels given in TAXA block");
        np.matchBeginBlock(NAME);

        np.pushPunctuationCharacters(NexusStreamParser.STRICT_PUNCTUATION);
        try {
            /* Read in the taxa sets */
            while (np.peekMatchIgnoreCase("taxset")) {
                np.matchIgnoreCase("taxset");
                String taxSetName = np.getWordRespectCase();
                np.matchIgnoreCase("=");
                /* Check that set does not have the same name as a taxon */
                //int taxindex = taxa.indexOf(taxSetName);
                if (taxa.indexOf(taxSetName) > 0)
                    throw new IOException("line " + np.lineno() +
                            ": Can't read SETS block because the label for a TAXSET is the same as a taxon name");
                Set thisSet = readTaxSet(np, taxa);
                this.addTaxSet(taxSetName, thisSet);
            }


            while (np.peekMatchIgnoreCase("charset")) {
                if (chars == null) {
                    throw new IOException("line " + np.lineno() +
                            ": Can't read SETS block because there has been no CHARACTERS block");
                }
                np.matchIgnoreCase("charset");
                String charSetName = np.getWordRespectCase();
                np.matchIgnoreCase("=");
                Set thisSet = readCharSet(np, chars);
                this.addCharSet(charSetName, thisSet);
            }
            while (np.peekMatchIgnoreCase("charpartition")) {
                if (chars == null) {
                    throw new IOException("line " + np.lineno() +
                            ": Can't read SETS block because there has been no CHARACTERS block");
                }
                np.matchIgnoreCase("charpartition");
                String charPartName = np.getWordRespectCase();
                np.matchIgnoreCase("=");
                Partition thisPartition = readCharPartition(np, chars);
                this.addCharPartition(charPartName, thisPartition);
            }
            np.matchEndBlock();
        } catch (IOException ex) {
            throw ex;
        } catch (RuntimeException ex) {
            throw ex;
        } finally {
            np.popPunctuationCharacters();
        }
    }
}
