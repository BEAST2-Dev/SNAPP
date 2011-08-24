package beast.app.draw;

import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Rectangle;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.util.*;
import java.util.regex.PatternSyntaxException;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import javax.swing.event.CellEditorListener;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;
import javax.swing.table.DefaultTableModel;
import javax.swing.table.JTableHeader;
import javax.swing.table.TableCellEditor;
import javax.swing.table.TableCellRenderer;
import javax.swing.table.TableColumnModel;

import beast.core.Input;
import beast.core.Plugin;
import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;

public class DataInputEditor extends InputEditor {
	private static final long serialVersionUID = 1L;
	List<TaxonSet> m_taxonset;
	List<Taxon> m_lineageset;
	Map<String,String> m_taxonMap;
	JTable m_table;
	DefaultTableModel m_model = new DefaultTableModel();
	
	JTextField filterEntry;
	String m_sFilter = ".*";
	int m_sortByColumn = 0;
	boolean m_bIsAscending = true;
	
	@Override
	public Class<?> type() {
		return snap.Data.class;
	}

	
	@Override
	public void init(Input<?> input, Plugin plugin, EXPAND bExpand, boolean bAddButtons) {
		m_input = input;
		m_plugin = plugin;
		List<TaxonSet> taxonset = ((snap.Data)input.get()).m_taxonsets.get(); 
		add(getContent(taxonset));
	}
	
	
	private Component getContent(List<TaxonSet> taxonset) {
		m_taxonset = taxonset;
		m_taxonMap = new HashMap<String, String>();
		m_lineageset = new ArrayList<Taxon>();
		for (Taxon taxonset2 : m_taxonset) {
			for (Taxon taxon : ((TaxonSet)taxonset2).m_taxonset.get()) {
				m_lineageset.add(taxon);
				m_taxonMap.put(taxon.getID(), taxonset2.getID());
			}
		}
		
		// set up table.
		// special features: background shading of rows
		// custom editor allowing only Date column to be edited.		
		m_model = new DefaultTableModel();
		m_model.addColumn("Lineage");
		m_model.addColumn("Taxon");
		taxonSetToModel();

		m_table = new JTable(m_model) {
			private static final long serialVersionUID = 1L;

			// method that induces table row shading 
			@Override
			public Component prepareRenderer (TableCellRenderer renderer,int Index_row, int Index_col) {
				Component comp = super.prepareRenderer(renderer, Index_row, Index_col);
				//even index, selected or not selected
				if (isCellSelected(Index_row, Index_col)) {
					comp.setBackground(Color.gray);
				} else 	if (Index_row % 2 == 0) {
					comp.setBackground(new Color(237,243,255));
				} else {
					comp.setBackground(Color.white);
				}
				return comp;
			}
		};
		
		// set up editor that makes sure only doubles are accepted as entry
		// and only the Date column is editable.
		m_table.setDefaultEditor(Object.class, new TableCellEditor() {
			JTextField m_textField = new JTextField();
			int m_iRow, m_iCol;
			@Override
			public boolean stopCellEditing() {
				m_table.removeEditor();
				String sText = m_textField.getText();
				System.err.println(sText);
				m_model.setValueAt(sText, m_iRow, m_iCol);
//				try {
//					Double.parseDouble(sText);
//				} catch (Exception e) {
//					return false;
//				}
				modelToTaxonset();
				return true;
			}
		
			@Override
			public boolean isCellEditable(EventObject anEvent) {
				return m_table.getSelectedColumn() == 1;
			}
			
			
			@Override
			public Component getTableCellEditorComponent(JTable table, Object value, boolean isSelected, int iRow, int iCol) {
				if (!isSelected) {
					return null;
				}
				m_iRow = iRow;
				m_iCol = iCol;
				m_textField.setText((String)value);
				return m_textField; 			
			}

			@Override
			public boolean shouldSelectCell(EventObject anEvent) {return false;}
			@Override
			public void removeCellEditorListener(CellEditorListener l) {}
			@Override
			public Object getCellEditorValue() {return null;}
			@Override
			public void cancelCellEditing() {}
			@Override
			public void addCellEditorListener(CellEditorListener l) {}
		
		});				

		JTableHeader header = m_table.getTableHeader();
		header.addMouseListener(new ColumnHeaderListener());

		JScrollPane pane = new JScrollPane(m_table);

		Box box = Box.createVerticalBox();
		box.add(createFilterBox());
		box.add(pane);
		box.add(createButtonBox());
		return box;
	}
	
	private Component createButtonBox() {
		Box buttonBox = Box.createHorizontalBox();
		
		JButton fillDownButton = new JButton("Fill down");
		fillDownButton.setToolTipText("replaces all taxons in selection with the one that is selected at the top");
		fillDownButton.addActionListener(new ActionListener() {
			
			@Override
			public void actionPerformed(ActionEvent e) {
				int [] rows = m_table.getSelectedRows();
				if (rows.length < 2) {
					return;
				}
				String sTaxon = (String) ((Vector) m_model.getDataVector().elementAt(rows[0])).elementAt(1);
				for (int i = 1; i < rows.length; i++) {
					m_model.setValueAt(sTaxon, rows[i], 1);
				}
				modelToTaxonset();
			}
		});
		
		JButton guessButton = new JButton("Guess");
		guessButton.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				guess();
			}
		});
		
		buttonBox.add(Box.createHorizontalGlue());
		buttonBox.add(fillDownButton);
		buttonBox.add(Box.createHorizontalGlue());
		buttonBox.add(guessButton);
		buttonBox.add(Box.createHorizontalGlue());
		return buttonBox;
	}
	
	

	public class ColumnHeaderListener extends MouseAdapter {
	    public void mouseClicked(MouseEvent evt) {
	        // The index of the column whose header was clicked
	        int vColIndex = m_table.getColumnModel().getColumnIndexAtX(evt.getX());
	        if (vColIndex == -1) {
	            return;
	        }
	        if (vColIndex != m_sortByColumn)
	        	m_sortByColumn = vColIndex;
	        else 
	        	m_bIsAscending = !m_bIsAscending;
	        taxonSetToModel();
	    }
	}

	private void guess() {
		GuessDlg dlg = new GuessDlg(this);
		String sPattern = dlg.showDialog("Guess taxon sets");
		if (sPattern != null) {
			try {
				((snap.Data)m_input.get()).guessTaxonSets(sPattern, 0);
				for (Taxon taxonset2 : m_taxonset) {
					for (Taxon taxon : ((TaxonSet)taxonset2).m_taxonset.get()) {
						m_lineageset.add(taxon);
						m_taxonMap.put(taxon.getID(), taxonset2.getID());
					}
				}
				taxonSetToModel();
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}
	
	class GuessDlg extends JDialog {
		String m_sPattern = "^(.+)[-_\\. ](.*)$";
		Component m_parent;
		Box guessPanel;
		ButtonGroup group;
        JRadioButton b1 = new JRadioButton("use everything");
        JRadioButton b2 = new JRadioButton("use regular expression");
		
		int m_location = 0;
		String m_sDelimiter = ".";
		JTextField regexpEntry;
		
		GuessDlg(Component parent) {
			m_parent = parent;
	        guessPanel = Box.createVerticalBox();

	        group = new ButtonGroup();
	        group.add(b1);
	        group.add(b2);
	        group.setSelected(b1.getModel(), true);
	        
	        guessPanel.add(createDelimiterBox(b1));
	        guessPanel.add(Box.createVerticalStrut(20));
	        guessPanel.add(createRegExtpBox(b2));
	        guessPanel.add(Box.createVerticalStrut(20));
		}

		
		private Component createDelimiterBox(JRadioButton b) {
			Box box = Box.createHorizontalBox();
			box.add(b);
			
			JComboBox combo = new JComboBox(new String[]{"after first","after last","before first","before last"});
			box.add(Box.createHorizontalGlue());
			box.add(combo);
			combo.addActionListener(new ActionListener() {
				@Override
				public void actionPerformed(ActionEvent e) {
					JComboBox combo = (JComboBox) e.getSource();
					m_location = combo.getSelectedIndex();
				}
			});
			
			JComboBox combo2 = new JComboBox(new String[]{".",",","_","-"," ","/",":",";"});
			box.add(Box.createHorizontalGlue());
			box.add(combo2);
			combo2.addActionListener(new ActionListener() {
				@Override
				public void actionPerformed(ActionEvent e) {
					JComboBox combo = (JComboBox) e.getSource();
					m_sDelimiter = (String) combo.getSelectedItem();
				}
			});
			box.add(Box.createHorizontalGlue());
			return box;
		}

		
		
		public Component createRegExtpBox(JRadioButton b) {
			Box box = Box.createHorizontalBox();
			box.add(b);
			regexpEntry = new JTextField();
			regexpEntry.setText(m_sPattern);
			regexpEntry.setColumns(30);
			regexpEntry.setToolTipText("Enter regular expression to match taxa");
			regexpEntry.setMaximumSize(new Dimension(1024, 20));
			box.add(Box.createHorizontalGlue());
			box.add(regexpEntry);
			box.add(Box.createHorizontalGlue());
			return box;
		}

		public String showDialog(String title) {

	        JOptionPane optionPane = new JOptionPane(guessPanel,
	                JOptionPane.PLAIN_MESSAGE,
	                JOptionPane.OK_CANCEL_OPTION,
	                null,
	                new String[] { "Cancel", "OK" },
	                "OK");
	        optionPane.setBorder(new EmptyBorder(12, 12, 12, 12));

	        final JDialog dialog = optionPane.createDialog(m_parent, title);
	        //dialog.setResizable(true);
	        dialog.pack();

	        dialog.setVisible(true);

//	        if (optionPane.getValue() == null) {
//	        	System.exit(0);
//	        }
	        if (b1.getModel() == group.getSelection()) {
	        	String sDelimiter = m_sDelimiter;
	        	if (sDelimiter.equals(".") || sDelimiter.equals("/")) {
	        		sDelimiter = "\\" + sDelimiter;
	        	}
	        	switch (m_location) {
	        	case 0:
	        		m_sPattern = "^[^"+sDelimiter+"]+" + sDelimiter + "(.*)$";
	        		break;
	        	case 1:
	        		m_sPattern = "^.*" + sDelimiter + "(.*)$";
	        		break;
	        	case 2:
	        		m_sPattern = "^([^"+sDelimiter+"]+)" + sDelimiter + ".*$";
	        		break;
	        	case 3:
	        		m_sPattern = "^(.*)" + sDelimiter + ".*$";
	        		break;
	        	}
	        }	        
	        if (b2.getModel() == group.getSelection()) {
	        	m_sPattern = regexpEntry.getText();
	        }	        
	        
	        // sanity check
	        try {
				m_sPattern.matches(m_sPattern);
			} catch (PatternSyntaxException e) {
				JOptionPane.showMessageDialog(this, "This is not a valid regular expression");
				return null;
			}
	        
	        if (optionPane.getValue().equals("OK")) {
	        	System.err.println("Pattern = " + m_sPattern);
	        	return m_sPattern;
	        } else {
	        	return null;
	        }
	    }
	}

	private Component createFilterBox() {
		Box filterBox = Box.createHorizontalBox();
		filterBox.add(new JLabel("filter: "));
		//Dimension size = new Dimension(100,20);
		filterEntry = new JTextField();
		filterEntry.setColumns(80);
//		filterEntry.setMinimumSize(size);
//		filterEntry.setPreferredSize(size);
//		filterEntry.setSize(size);
		filterEntry.setToolTipText("Enter regular expression to match taxa");
		filterEntry.setMaximumSize(new Dimension(1024, 20));
		filterBox.add(filterEntry);
		filterBox.add(Box.createHorizontalGlue());
		filterEntry.getDocument().addDocumentListener(new DocumentListener() {
			@Override
			public void removeUpdate(DocumentEvent e) {
				processFilter();
			}
			@Override
			public void insertUpdate(DocumentEvent e) {
				processFilter();
			}
			@Override
			public void changedUpdate(DocumentEvent e) {
				processFilter();
			}
			private void processFilter() {
				String sFilter = ".*" + filterEntry.getText() + ".*";
				try {
					// sanity check: make sure the filter is legit
					sFilter.matches(sFilter);
					m_sFilter = sFilter;
					taxonSetToModel();
					m_table.repaint();
				} catch (PatternSyntaxException e) {
					// ignore
				}
			}
		});
		return filterBox;
	}




	/** for convert taxon sets to table model **/
	private void taxonSetToModel() {
		// count number if lineages that match the filter
		int i = 0;
		for (String sLineageID : m_taxonMap.keySet()) {
			if (sLineageID.matches(m_sFilter)) {
				i++;
			}
		}
		
		// clear table model
		while (m_model.getRowCount() > 0) {
			m_model.removeRow(0);
		}	
		
		// fill table model with lineages matching the filter
		for (String sLineageID : m_taxonMap.keySet()) {
			if (sLineageID.matches(m_sFilter)) {
				Object [] rowData = new Object[2];
				rowData[0] = sLineageID;
				rowData[1] = m_taxonMap.get(sLineageID);
				m_model.addRow(rowData);
			}
		}
		
	    Vector data = m_model.getDataVector();
	    Collections.sort(data, new Comparator<Vector>() {
			@Override
			public int compare(Vector v1, Vector v2) {
		        String o1 = (String) v1.get(m_sortByColumn);
		        String o2 = (String) v2.get(m_sortByColumn);
		        if (o1.equals(o2)) {
			        o1 = (String) v1.get(1 - m_sortByColumn);
			        o2 = (String) v2.get(1 - m_sortByColumn);
		        }
		        if (m_bIsAscending) {
		        	return o1.compareTo(o2);
		        } else {
		        	return o2.compareTo(o1);
		        }
			}
	    	
		});
//	    m_model.fireTableStructureChanged();

	    m_model.fireTableRowsInserted(0, m_model.getRowCount());
	}

	/** for convert table model to taxon sets **/
	private void modelToTaxonset() {
		// update map
		for (int i = 0; i < m_model.getRowCount();i++) {
			String sLineageID = (String) ((Vector)m_model.getDataVector().elementAt(i)).elementAt(0);
			String sTaxonSetID = (String) ((Vector)m_model.getDataVector().elementAt(i)).elementAt(1);
			
			// new taxon set?
			if (!m_taxonMap.containsValue(sTaxonSetID)) {
				// create new taxon set
				TaxonSet taxonset = new TaxonSet();
				taxonset.setID(sTaxonSetID);
				m_taxonset.add(taxonset);
			}
			m_taxonMap.put(sLineageID, sTaxonSetID);
		}	
		
		// clear old taxon sets
		for (TaxonSet set : m_taxonset) {
			set.m_taxonset.get().clear();
		}
		
		// group lineages with their taxon sets
		for (String sLineageID : m_taxonMap.keySet()) {
			for (Taxon taxon : m_lineageset) {
				if (taxon.getID().equals(sLineageID)) {
					String sTaxonSet = m_taxonMap.get(sLineageID);
					for (TaxonSet set : m_taxonset) {
						if (set.getID().equals(sTaxonSet)) {
							try {
								set.m_taxonset.setValue(taxon, set);
							} catch (Exception e) {
								e.printStackTrace();
							}
						}
					}
				}
			}
		}
		
		// remove unused taxon sets
		for (int i = m_taxonset.size()-1; i >= 0; i--) {
			if (m_taxonset.get(i).m_taxonset.get().size() == 0) {
				m_taxonset.remove(i);
			}
		}
	}

}
