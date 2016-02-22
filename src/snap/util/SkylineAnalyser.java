package snap.util;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.Arrays;

import beast.core.Description;
import beast.core.Input;
import beast.core.Runnable;
import beast.core.util.Log;
import beast.util.LogAnalyser;

@Description("Creates skyline plot from trace log")
public class SkylineAnalyser extends Runnable {
	public Input<File> traceLogInput = new Input<>("tracelog", "file containing the trace log to be analysed", new File("trace.log"));
	public Input<File> outputInput = new Input<>("out", "file where results are being written to", new File("/tmp/skyline.log"));
	public Input<String> nodeNameInput = new Input<>("nodeName", "name of the node heights in the trace log (minus the number)", "node");
	public Input<String> thetaNameInput = new Input<>("thetaName", "name of the thetas in the trace log (minus the number)", "theta");
	public Input<Integer> groupSizeInput = new Input<>("groupSize", "number of groups = number of branches in tree", 2);
	public Input<Double> totalHeightInput = new Input<>("totalHeight", "total height of the tree", 0.1);
	public Input<Integer> numberOfBinsInput = new Input<>("numberOfBins", "number of bins used for reconstructing the skyline", 100);
	public Input<Integer> burnInPercentageInput = new Input<>("burnInPercentage", "percentage of samples in trace log to be considered burnin", 10);


	@Override
	public void initAndValidate() {

	}

	@Override
	public void run() throws Exception {
		LogAnalyser trace = new LogAnalyser(traceLogInput.get().getPath(), burnInPercentageInput.get());
		int traceLength = trace.getTrace(0).length;
		
		Double [][] nodeHeights0 = new Double[groupSizeInput.get() + 1][];
		Double [][] thetas0 = new Double[groupSizeInput.get() + 1][];
		for (int i = 0; i < groupSizeInput.get() + 1; i++) {
			nodeHeights0[i] = trace.getTrace(nodeNameInput.get() + i);
			thetas0[i] = trace.getTrace(thetaNameInput.get() + i);
		}
		// transpose
		double [][] nodeHeights = new double[traceLength][groupSizeInput.get() + 1];
		double [][] thetas = new double[traceLength][groupSizeInput.get() + 1];
		for (int i = 0; i < groupSizeInput.get() + 1; i++) {
			for (int j = 0; j < traceLength; j++) {
				nodeHeights[j][i] = nodeHeights0[i][j];
				thetas[j][i] = thetas0[i][j];
			}		
		}
		// save memory
		nodeHeights0 = null;
		thetas0 = null;
		trace = null;
		
		
		int bins = numberOfBinsInput.get();
		
		double[] heights = new double[bins];
		for (int i = 0; i < bins; i++) {
			heights[i] = i * totalHeightInput.get() / (bins -1);
		}
		
		double[][] plot = new double [bins][traceLength];
		for (int i = 0; i < traceLength; i++) {
			populate(plot, i,  nodeHeights[i], thetas[i], heights);
		}
		for (int i = 0; i < bins; i++) {
			Arrays.sort(plot[i]);
		}
		createOutput(plot, heights, traceLength);
		
		Log.info.println("Done!");
	}

	private void createOutput(double[][] plot, double[] heights, int traceLength) throws FileNotFoundException {
        PrintStream outfile = new PrintStream(outputInput.get());
        
		int medianIndex = traceLength / 2;
		int lowerIndex = (int) (traceLength * 0.025);
		int upperIndex = (int) (traceLength * 0.975);
		outfile.println("height\tlower\tmedian\tupper");
		for (int i = 0; i < heights.length; i++) {
			outfile.print(heights[i] + "\t" + plot[i][lowerIndex] + "\t");
			outfile.print(plot[i][medianIndex] + "\t");
			outfile.print(plot[i][upperIndex] + "\n");
		}
        outfile.close();
	}

	private void populate(double[][] plot, int i, double[] nodeHeights, double[] thetas, double[] heights) {
		int h = 0;
		double theta = thetas[h];
		int k = 0;
		int bins = plot.length;
		while (k < bins) {
			plot[k][i] = theta;
			if (h < nodeHeights.length && heights[k] > nodeHeights[h]) {
				theta = thetas[h];
				h++;
			}
			k++;
		}
	}

}
