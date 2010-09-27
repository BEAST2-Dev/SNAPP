package snap.util;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.Box;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JSlider;
import javax.swing.JTextField;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import snap.FMatrix;
import snap.likelihood.MatrixExponentiator;

/** panel for showing distribution of # reds for given start # reds for 
 * values of U, V, Gamma over a time range
 */
public class QViewer extends JPanel implements ChangeListener {
	private static final long serialVersionUID = 1L;
	Box box;
	JTextField m_entryNBottom;
	int m_nBottom = 2;
	JTextField m_entryNTop;
	int m_nTop = 2;
	JTextField m_entryU;
	double m_fU = 1;
	JTextField m_entryV;
	double m_fV = 1;
	JTextField m_entryGamma;
	double m_fMaxGamma = 0.01;
	double m_fGamma = m_fMaxGamma;
	JTextField m_entryT;
	double m_fT = 1;
	JTextField m_entryScale;
	double m_fScale = 1;
	JSlider m_slider;

	Color [] m_color;
	
	public QViewer() {
		m_entryNBottom = new JTextField(m_nBottom+"");
		m_entryNBottom.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				m_nBottom = Integer.parseInt(m_entryNBottom.getText());
				recalc();
			}
		});		
		m_entryNTop = new JTextField(m_nTop+"");
		m_entryNTop.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				m_nTop = Integer.parseInt(m_entryNTop.getText());
				recalc();
			}
		});		
		m_entryU = new JTextField(m_fU+"");
		m_entryU.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				m_fU = Double.parseDouble(m_entryU.getText());
				recalc();
			}
		});		
		m_entryV = new JTextField(m_fV+"");
		m_entryV.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				m_fV = Double.parseDouble(m_entryV.getText());
				recalc();
			}
		});		
		m_entryGamma = new JTextField(m_fGamma+"");
		m_entryGamma.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				m_fMaxGamma = Double.parseDouble(m_entryGamma.getText());
				m_fGamma = m_fMaxGamma; 
				recalc();
			}
		});		
		m_entryT = new JTextField(m_fT+"");
		m_entryT.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				m_fT = Double.parseDouble(m_entryT.getText());
				recalc();
			}
		});		
		m_entryScale = new JTextField(m_fScale+"");
		m_entryScale.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				m_fScale = Double.parseDouble(m_entryScale.getText());
				recalc();
			}
		});		
		
		m_slider = new JSlider(JSlider.HORIZONTAL, 0, 1000, 1000);
		m_slider.setPreferredSize(new Dimension(100, 5));
		m_slider.addChangeListener(this);
		
		box = Box.createHorizontalBox();
		box.add(new JLabel("#lineages at bottom"));
		box.add(m_entryNBottom);
		box.add(new JLabel("#lineages at top"));
		box.add(m_entryNTop);
		box.add(new JLabel("gamma"));
		box.add(m_entryGamma);
		box.add(new JLabel("u"));
		box.add(m_entryU);
		box.add(new JLabel("v"));
		box.add(m_entryV);
		box.add(new JLabel("t"));
		box.add(m_entryT);
		box.add(new JLabel("Scale:"));
		box.add(m_entryScale);
        box.add(m_slider);

		m_color = new Color[12];
		int k = 0;
		m_color[k++] = Color.blue;
		m_color[k++] = Color.green;
		m_color[k++] = Color.red;
		m_color[k++] = Color.gray;
		m_color[k++] = Color.orange;
		m_color[k++] = Color.yellow;
		m_color[k++] = Color.pink;
		m_color[k++] = Color.black;
		m_color[k++] = Color.cyan;
		m_color[k++] = Color.darkGray;
		m_color[k++] = Color.magenta;
	} // c'tor
	
	@Override
	public void stateChanged(ChangeEvent e) {
		int nPercentage = m_slider.getValue();
		m_fGamma = m_fMaxGamma * nPercentage / 100;
		recalc();
	}

	
	double [][][]m_fP; // # bottom lineages x # top lineages x # time stamps
	static int NR_TIMES = 50;
	void recalc() {
		m_fP = new double[m_nBottom+1][NR_TIMES][m_nTop+1];
		for (int nReds = 0; nReds <= m_nBottom; nReds++) {
			FMatrix bottomDistribution = new FMatrix(m_nBottom, nReds);
			for (int iTime = 1; iTime <= NR_TIMES; iTime++) {
				try {
					FMatrix tmp = MatrixExponentiator.expQTtx(m_nBottom, m_fU, m_fV, m_fGamma, m_fT * iTime / NR_TIMES, bottomDistribution);
					for (int iTop = 0; iTop <= m_nTop; iTop++) {
						m_fP[nReds][iTime-1][iTop] = tmp.get(m_nTop, iTop); 
					}
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
		}
		repaint();
	}
	

	protected void paintComponent(Graphics g) {
		Graphics2D g2 = (Graphics2D) g;
		g2.setBackground(Color.WHITE);
		int nWidth = getWidth();
		int nHeight = getHeight();
		g2.clearRect(0, 0, nWidth, nHeight);
		g2.drawString("gamma=" + m_fGamma, 10, 10);
		if (m_fP != null) {
			for (int iBottom = 0; iBottom <= m_nBottom; iBottom++) {
				g2.setColor(m_color[iBottom % m_color.length]);
				g2.drawString(iBottom + " ", 10, 10*iBottom+20);
				for (int iTop = 0; iTop <= m_nTop; iTop++) {
					for (int iTime = 1; iTime < NR_TIMES; iTime++) {
						int x1 = nWidth * (iTime - 1) / NR_TIMES;
						int y1 = nHeight-(int)(m_fP[iBottom][iTime-1][iTop] * nHeight * m_fScale);
						int x2 = nWidth * iTime / NR_TIMES;
						int y2 = nHeight-(int)(m_fP[iBottom][iTime][iTop] * nHeight * m_fScale);
						g2.drawLine(x1, y1, x2, y2);
					}
				}
			}
		}
		
	} // paintComponent
	
	public static void main(String [] args) {
		JFrame frame = new JFrame("RateMatrixBySampling");		
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		QViewer viewer =  new QViewer();
		viewer.setVisible(true);
		frame.add(viewer,BorderLayout.CENTER);
		frame.add(viewer.box, BorderLayout.SOUTH);
		frame.setSize(600, 800);
		frame.setVisible(true);
	} // main
	
} // class QViewer
