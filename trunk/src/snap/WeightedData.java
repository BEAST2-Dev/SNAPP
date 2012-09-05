package snap;

import java.util.Arrays;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.util.Randomizer;

@Description("Alignment with weights attached to patterns")
public class WeightedData extends Data {
	public Input<String> weightInput = new Input<String>("weights", "comma separated list of weights, one for each site in the alignment", Validate.REQUIRED);
	// balance
	public Input<Boolean> balance = new Input<Boolean>("balanced", "whether to attempt to reassign labels so #reds approx equal to #greens. " +
			"This assumes there is no significance in colour of the labellings. " +
			"Patterns and complementary patterns are found and weights swapped", false);

	@Override
	public void initAndValidate() throws Exception {
		super.initAndValidate();
		
		String [] sStr = weightInput.get().split(",");
		
		int siteCount = getSiteCount() + m_nWeight[m_nWeight.length - 2] + m_nWeight[m_nWeight.length - 1];
		if (sStr.length != siteCount) {
			throw new Exception("Number of weights (" + sStr.length + ") does not match number of sites (" + siteCount + ") in alignment");
		}
		Arrays.fill(m_nWeight, 0);
		for (int i = 0; i < sStr.length; i++) {
			int userWeight = Integer.parseInt(sStr[i]);
			int pattern = getPatternIndex(i);
			if (pattern >= m_nWeight.length - 2) {
				System.err.println("WARNING: Constant pattern detected. This will be ignored.");
			} else {
				m_nWeight[pattern] += userWeight;
			}
		}
		if (balance.get()) {
			attemptToBalance();
		}
	} // initAndValidate

	private void attemptToBalance() {
		int nPatterns = m_nPatterns.length-2;
		int [] complement = new int[nPatterns];
		Arrays.fill(complement, -1);

		// find complement of each pattern
		for (int i = 0; i < nPatterns; i++) {
			int [] pattern0 = getPattern(i);
			for (int j = i+1; j < nPatterns; j++) {
				int [] pattern1 = getPattern(j);
				boolean match = true;
				for (int k = 0; k < pattern0.length; k++) {
					if (pattern0[k] + pattern1[k] != m_nStateCounts.get(k)) {
						match = false;
						break;
					}
				}
				if (match) {
					complement[i] = j;
					complement[j] = i;
				}
			}
		}

		// determine contribution of each pattern
		int [] patternscore = new int[nPatterns];
		int [] score = new int[nPatterns];
		int currentScore = 0;
		for (int i = 0; i < nPatterns; i++) {
			int [] pattern0 = getPattern(i);
			for (int k = 0; k < pattern0.length; k++) {
				patternscore[i] += pattern0[k];
			}
			score[i] = patternscore[i] * m_nWeight[i];
			currentScore += score[i]; 
		}	

		// calc target score
		int targetScore = 0;
		int totalLineages = 0;
		for (int k = 0; k < getPattern(0).length; k++) {
			totalLineages += m_nStateCounts.get(k);
		}
		for (int i = 0; i < nPatterns; i++) {
			targetScore += totalLineages * m_nWeight[i];
		}
		targetScore = targetScore / 2;
		
		// swap weights
		int delta = targetScore - currentScore;
		int attempts = 0;
		while (delta != 0 && attempts < 10 * nPatterns) {
			int i = Randomizer.nextInt(nPatterns);
			if (complement[i] >= 0) {
				int j = complement[i];
				int newScore = currentScore - m_nWeight[i] * patternscore[i] - m_nWeight[j] * patternscore[j] +
						m_nWeight[i] * patternscore[j] + m_nWeight[j] * patternscore[i];
				int newDelta = targetScore - newScore;
				if (Math.abs(newDelta) < Math.abs(delta)) {
					int tmp = m_nWeight[i];
					m_nWeight[i] = m_nWeight[j];
					m_nWeight[j] = tmp;
					delta = newDelta;
					currentScore = newScore;
				}
			}
			attempts++;
		}
		
		System.out.println(Arrays.toString(m_nWeight));
	} // attemptToBalance
	

	@Override
	public String toString() {
		StringBuffer buf = new StringBuffer();
		if (m_nPatterns != null) {
			int nTaxa = m_nPatterns[0].length;
			int nPatterns = m_nPatterns.length;
			for (int i = 0; i < nTaxa; i++) {
				for (int j = 0; j < nPatterns; j++) {
					buf.append(m_nPatterns[j][i]);
				}
				buf.append("\n");
			}
			if (m_nWeight != null) {
				for (int j = 0; j < nPatterns; j++) {
					buf.append(m_nWeight[j] + ",");
				}
				buf.append("\n");
			}
		}
		return buf.toString();
	}
}
