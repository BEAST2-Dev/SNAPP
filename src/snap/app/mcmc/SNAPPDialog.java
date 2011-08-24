package snap.app.mcmc;

import jam.mac.Utils;
import jam.panels.OptionsPanel;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import javax.swing.border.TitledBorder;
import javax.swing.filechooser.FileFilter;

import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;

public class SNAPPDialog {
	private final JFrame frame;

	private final OptionsPanel optionPanel;

	private final WholeNumberField seedText = new WholeNumberField((long) 1, Long.MAX_VALUE);
	// private final JCheckBox overwriteCheckBox = new
	// JCheckBox("Allow overwriting of log files");
	private final JComboBox logginMode = new JComboBox(new String[] { "default: only write new log files",
			"overwrite: overwrite log files", "resume: appends log to existing files (if any)" });

	private final JComboBox threadsCombo = new JComboBox(new Object[] { "Automatic", 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11,
			12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24 });

	private File inputFile = null;

	public SNAPPDialog(final JFrame frame, final String titleString, final Icon icon) {
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
		final JTextField inputFileNameText = new JTextField("not selected", 16);

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
							if (name.endsWith(".xml")) {
								return true;
							}
							return false;
						}

						// The description of this filter
						public String getDescription() {
							return "xml files";
						}
					});

					fc.setDialogTitle("Load xml file");
					int rval = fc.showOpenDialog(null);
					if (rval == JFileChooser.APPROVE_OPTION) {
						inputFile = fc.getSelectedFile();
						inputFileNameText.setText(inputFile.getName());
					}
				} else {
					FileDialog dialog = new FileDialog(frame, "Select SNAPP xml file...", FileDialog.LOAD);

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
		inputFileNameText.setEditable(false);

		JPanel panel1 = new JPanel(new BorderLayout(0, 0));
		panel1.add(inputFileNameText, BorderLayout.CENTER);
		panel1.add(inputFileButton, BorderLayout.EAST);
		optionPanel.addComponentWithLabel("SNAPP XML File: ", panel1);

		optionPanel.addComponent(logginMode);
		// optionPanel.addComponent(overwriteCheckBox);

		optionPanel.addSeparator();

		seedText.setColumns(12);
		optionPanel.addComponentWithLabel("Random number seed: ", seedText);

		optionPanel.addComponentWithLabel("Thread pool size: ", threadsCombo);
		threadsCombo.setSelectedIndex(1);

		optionPanel.addSeparator();

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

		seedText.setValue(seed);

		final JDialog dialog = optionPane.createDialog(frame, title);
		// dialog.setResizable(true);
		dialog.pack();

		dialog.setVisible(true);

		if (optionPane.getValue() == null) {
			System.exit(0);
		}

		return optionPane.getValue().equals("Run");
	}

	public long getSeed() {
		return seedText.getLongValue();
	}

	// public boolean allowOverwrite() {
	// return overwriteCheckBox.isSelected();
	// }

	public int getLogginMode() {
		return logginMode.getSelectedIndex();
	}

	public int getThreadPoolSize() {
		if (threadsCombo.getSelectedIndex() == 0) {
			// Automatic
			return -1;
		}
		return (Integer) threadsCombo.getSelectedItem();
	}

	public File getInputFile() {
		return inputFile;
	}
}