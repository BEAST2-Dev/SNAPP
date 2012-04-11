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
		m_sData.setRule(Validate.OPTIONAL);
	}
	
	@Override
	public void initAndValidate() throws Exception {
		m_nTotalCount.setValue(m_sequences.get().size(), this);
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
                	throw new Exception("sequence " + sequence.m_sTaxon.get() + " is of length " + sequence2.size() + 
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
	 *  Ambiguous sites are disgarded **/
    public List<Integer> getStateCounts(DataType dataType) throws Exception {
    	int [] statecounts = null;
        
        // grab info from sub sequences and add them up
        for (Sequence sequence: m_sequences.get()) {
        	if (statecounts == null) {
        		// start with first sequence
        		Integer[] sequence0 = sequence.getSequence(dataType).toArray(new Integer[0]);
        		statecounts = new int[sequence0.length];
                for (int i = 0; i < statecounts.length; i++) {
                	statecounts[i] = (dataType.isAmbiguousState(sequence0[i])? 0 : 1);
                }
        	} else {
                List<Integer> sequence2 = sequence.getSequence(dataType);
                // sanity check: make sure sequence2 is of same length as rest
                if (sequence2.size() != statecounts.length) {
                	throw new Exception("sequence " + sequence.m_sTaxon.get() + " is of length " + sequence2.size() + 
                			" but it was expected to be of length " + statecounts.length);
                }
                // add instances from sequence2 to sequences
                for (int i = 0; i < statecounts.length; i++) {
                	statecounts[i] += (dataType.isAmbiguousState(sequence2.get(i))? 0 : 1);
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
