package snap;

import java.util.ArrayList;
import java.util.List;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.evolution.alignment.Sequence;

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
	public List<Integer> getSequence(String sDataMap) throws Exception {
        Integer [] sequences = null;
        
        // grab info from sub sequences and add them up
        for (Sequence sequence: m_sequences.get()) {
        	if (sequences == null) {
        		sequences = sequence.getSequence(sDataMap).toArray(new Integer[0]);
        	} else {
                List<Integer> sequence2 = sequence.getSequence(sDataMap);
                // sanity check: make sure sequence2 is of same length as rest
                if (sequence2.size() != sequences.length) {
                	throw new Exception("sequence " + sequence.m_sTaxon.get() + " is of length " + sequence2.size() + 
                			" but it was expected to be of length " + sequences.length);
                }
                for (int i = 0; i < sequences.length; i++) {
                	sequences[i] += sequence2.get(i);
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
	
	
} // class SNPSequence
