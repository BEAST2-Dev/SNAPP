package snap.app.analysis;

import jam.mac.Utils;
import jam.panels.OptionsPanel;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import javax.swing.border.TitledBorder;
import javax.swing.filechooser.FileFilter;

import snap.app.mcmc.WholeNumberField;

import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;

public class TSADialog {
	private final JFrame frame;

	private final OptionsPanel optionPanel;

	private File inputFile = null;

	final WholeNumberField burninText = new WholeNumberField((long) 10, Long.MAX_VALUE);
	final JTextField newickText = new JTextField();
	final JTextField inputFileNameText = new JTextField("not selected", 16);

	public TSADialog(final JFrame frame, final String titleString, final Icon icon) {
		this.frame = frame;

		optionPanel = new OptionsPanel(12, 12);

		// this.frame = frame;

		JPanel panel = new JPanel(new BorderLayout());
		panel.setOpaque(false);

		final JLabel titleText = new JLabel(titleString);
		titleText.setIcon(icon);
		optionPanel.addSpanningComponent(titleText);
		titleText.setFont(new Font("sans-serif", 0, 12));

		final JButton inputFileButton = new JButton("Choose File...");

		inputFileButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent ae) {
				if (!Utils.isMacOSX()) {
					JFileChooser fc = new JFileChooser(System.getProperty("user.dir"));
					fc.addChoosableFileFilter(new FileFilter() {
						public boolean accept(File f) {
							if (f.isDirectory()) {
								return true;
							}
							String name = f.getName().toLowerCase();
							if (name.endsWith(".trees")) {
								return true;
							}
							return false;
						}

						// The description of this filter
						public String getDescription() {
							return "tree set files";
						}
					});

					fc.setDialogTitle("Load tree set file");
					int rval = fc.showOpenDialog(null);
					if (rval == JFileChooser.APPROVE_OPTION) {
						inputFile = fc.getSelectedFile();
						inputFileNameText.setText(inputFile.getName());
					}
				} else {
					FileDialog dialog = new FileDialog(frame, "Select tree set file...", FileDialog.LOAD);

					dialog.setVisible(true);
					if (dialog.getFile() == null) {
						// the dialog was cancelled...
						return;
					}

					inputFile = new File(dialog.getDirectory(), dialog.getFile());
					inputFileNameText.setText(inputFile.getName());

				}

			}
		});
		inputFileNameText.setEditable(true);

		JPanel panel1 = new JPanel(new BorderLayout(0, 0));
		panel1.add(inputFileNameText, BorderLayout.CENTER);
		panel1.add(inputFileButton, BorderLayout.EAST);
		optionPanel.addComponentWithLabel("Tree set file: ", panel1);

		// optionPanel.addComponent(overwriteCheckBox);

		optionPanel.addSeparator();

		// burninText.setValue(10);
		burninText.setColumns(12);
		optionPanel.addComponentWithLabel("Burn in percentage: ", burninText);

		optionPanel.addSeparator();
		newickText.setColumns(60);
		newickText.setToolTipText("check whether this tree is in the 95% credible set");
		optionPanel.addComponentWithLabel("Tree (optional): ", newickText);

		final OptionsPanel optionPanel1 = new OptionsPanel(0, 12);
		// optionPanel1.setBorder(BorderFactory.createEmptyBorder());
		optionPanel1.setBorder(new TitledBorder(""));

		OptionsPanel optionPanel2 = new OptionsPanel(0, 12);
		optionPanel2.setBorder(BorderFactory.createEmptyBorder());

		optionPanel1.addComponent(optionPanel2);

		// optionPanel.addSpanningComponent(optionPanel1);

	}

	public boolean showDialog(String title, long seed) {

		JOptionPane optionPane = new JOptionPane(optionPanel, JOptionPane.PLAIN_MESSAGE, JOptionPane.OK_CANCEL_OPTION,
				null, new String[] { "Run", "Quit" }, "Run");
		optionPane.setBorder(new EmptyBorder(12, 12, 12, 12));

		final JDialog dialog = optionPane.createDialog(frame, title);
		// dialog.setResizable(true);
		dialog.pack();

		dialog.setVisible(true);

		if (optionPane.getValue() == null) {
			System.exit(0);
		}

		return optionPane.getValue().equals("Run");
	}

	public File getInputFile() {
		return inputFile;
	}
}