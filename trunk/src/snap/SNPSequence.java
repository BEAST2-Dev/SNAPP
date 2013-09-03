package snap;

import java.util.ArrayList;
import java.util.List;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.evolution.alignment.Sequence;
import beast.evolution.datatype.DataType;



@Description("A SNPSequence holds a collection of binary sequences that are summarized in a new sequence " +
		"by counting the number of sites with a 1.")
public class SNPSequence extends Sequence {
	public Input<List<Sequence>> m_sequences = new Input<List<Sequence>>("sequence","binary sequence", new ArrayList<Sequence>());

	public SNPSequence() {
		dataInput.setRule(Validate.OPTIONAL);
	}
	
	@Override
	public void initAndValidate() throws Exception {
		int nNrOfStates = 0;
		
		boolean warned = false;
		for (Sequence s : m_sequences.get()) {
			// remove one for zero state
			nNrOfStates += s.totalCountInput.get() - 1;
			
			if (!warned && s.getClass().equals(Sequence.class) && s.totalCountInput.get() != null && s.totalCountInput.get() > 3) {
				warned = true;
				System.out.println("WARNING: there is a sequence with a state count larger than 2.\n" +
					"If this XML file was generated with BEAUti, this is most likely wrong.\n" +
					"To fix this, change totalCount in the XML to 2 for binary sequences, or 3 for diploid data.");
			}
			if (!warned && s.getClass().equals(Sequence.class) && s.totalCountInput.get() != null && s.totalCountInput.get() > 3) {
				System.out.println("WARNING: there is a sequence with a state count of 1 or less, which will be ignored.");
			}
			
		}
		// add one for zero state
		totalCountInput.setValue(nNrOfStates + 1, this);
	}

	@Override
    public List<Integer> getSequence(DataType dataType) throws Exception {
        Integer [] sequences = null;
        
        // grab info from sub sequences and add them up
        for (Sequence sequence: m_sequences.get()) {
        	if (sequences == null) {
        		// start with first sequence
        		sequences = sequence.getSequence(dataType).toArray(new Integer[0]);
                for (int i = 0; i < sequences.length; i++) {
                	sequences[i] = Math.max(0, sequences[i]);
                }
        	} else {
                List<Integer> sequence2 = sequence.getSequence(dataType);
                // sanity check: make sure sequence2 is of same length as rest
                if (sequence2.size() != sequences.length) {
                	throw new Exception("sequence " + sequence.taxonInput.get() + " is of length " + sequence2.size() + 
                			" but it was expected to be of length " + sequences.length);
                }
                // add instances from sequence2 to sequences
                for (int i = 0; i < sequences.length; i++) {
                	sequences[i] += Math.max(0, sequence2.get(i));
                }
        	}
        }
        
        // convert array to list
        List<Integer> list = new ArrayList<Integer>();
        for (Integer i : sequences) {
        	list.add(i);
        }
        return list;
	}
	
	/** For every site in the sequence, determine the total nr of reds + greens 
	 *  Ambiguous sites are disregarded **/
    public List<Integer> getLineageCounts(DataType dataType) throws Exception {
    	int [] statecounts = null;
        
        // grab info from sub sequences and add them up
        for (Sequence sequence: m_sequences.get()) {
        	System.out.println("Seq: " + sequence.getID());
        	if (statecounts == null) {
        		int maxStateCount = sequence.totalCountInput.get();
        		// start with first sequence
        		Integer[] sequence0 = sequence.getSequence(dataType).toArray(new Integer[0]);
        		statecounts = new int[sequence0.length];
                for (int i = 0; i < statecounts.length; i++) {
                	statecounts[i] = (dataType.isAmbiguousState(sequence0[i])? 0 : maxStateCount - 1);
                }
        	} else {
        		int maxStateCount = sequence.totalCountInput.get();
                List<Integer> sequence2 = sequence.getSequence(dataType);
                // sanity check: make sure sequence2 is of same length as rest
                if (sequence2.size() != statecounts.length) {
                	throw new Exception("sequence " + sequence.taxonInput.get() + " is of length " + sequence2.size() + 
                			" but it was expected to be of length " + statecounts.length);
                }
                // add instances from sequence2 to sequences
                for (int i = 0; i < statecounts.length; i++) {
                	statecounts[i] += (dataType.isAmbiguousState(sequence2.get(i))? 0 : maxStateCount - 1);
                }
        	}
        }
        // convert array to list
        List<Integer> list = new ArrayList<Integer>();
        for (Integer i : statecounts) {
        	list.add(i);
        }
        return list;
	}
	
} // class SNPSequence
