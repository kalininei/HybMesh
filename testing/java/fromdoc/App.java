import java.util.ArrayList;
import javax.swing.*;
import java.beans.*;
import java.awt.*;
import java.awt.event.*;

public class App extends JFrame implements ActionListener{
	public static final long serialVersionUID = 1L;

	// Entry point
	public static void main(String[] args){
		try{ 
			UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName()); 
			// set path to currently used hybmesh executable
			Hybmesh.hybmesh_exec_path = "../../../src/py/hybmesh.py";
			GridData gd = new GridData();
			javax.swing.SwingUtilities.invokeLater(new Runnable() {
				public void run() { createAndShowGUI(gd); }
			});
		} catch(Exception e){ 
			e.printStackTrace(); 
		} 
	}

	Drawer tab1, tab2;
	GridData gd;
	JToggleButton togAdd, togMove, togRem;

	private static void createAndShowGUI(GridData gd){
		JFrame app = new App(gd);
		app.pack();
		app.setLocationRelativeTo(null);
		app.setVisible(true);
	}

	App(GridData gd){
		super("Hybmesh/java mapping demo");
		this.gd = gd;
		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		//Create the toolbar.
		JToolBar toolbar = new JToolBar();
		toolbar.setFloatable(false);
		addButtons(toolbar);
	
		//Create split pane
		tab1 = new Drawer(gd.from, this);
		tab2 = new Drawer(gd.to, this);
		JSplitPane splitPane = new JSplitPane(
			JSplitPane.HORIZONTAL_SPLIT, tab1, tab2);
		splitPane.setDividerLocation(400);

		//main panel
		add(new JPanel(new BorderLayout()));
		getContentPane().setPreferredSize(new Dimension(800, 400));
		getContentPane().add(toolbar, BorderLayout.PAGE_START);
		getContentPane().add(splitPane, BorderLayout.CENTER);
	}

	// checks for what toggle button is pressed now. Used in Drawer mouse click handle.
	public boolean statusAddRef(){ return togAdd.isSelected(); }
	public boolean statusMoveRef(){ return togMove.isSelected(); }
	public boolean statusRemoveRef(){ return togRem.isSelected(); }

	void addButtons(JToolBar jtb){
		setButton(togAdd = new JToggleButton("Add refpoints"), "add", jtb);
		setButton(togMove = new JToggleButton("Move refpoints"), "move", jtb);
		setButton(togRem = new JToggleButton("Remove refpoints"), "remove", jtb);
		jtb.addSeparator();
		setButton(new JButton("Basic grid"), "basic", jtb);
		setButton(new JButton("Run mapping"), "run", jtb);
	}

	void setButton(AbstractButton b, String act, JToolBar jtb){
		jtb.add(b);
		b.setActionCommand(act);
		b.addActionListener(this);
	}

	public void actionPerformed(ActionEvent e) {
		String cmd = e.getActionCommand();
		if ("add".equals(cmd)){
			togMove.setSelected(false); togRem.setSelected(false);
		} else if ("move".equals(cmd)){
			togAdd.setSelected(false); togRem.setSelected(false);
		} else if ("remove".equals(cmd)){
			togAdd.setSelected(false); togMove.setSelected(false);
		} else if ("basic".equals(cmd)){
			//Build a dialog to ask new dimensions of the basic grid
			JSpinner nxField = new JSpinner(new SpinnerNumberModel(10, 1, 99, 1));
			JSpinner nyField = new JSpinner(new SpinnerNumberModel(10, 1, 99, 1));
			final JComponent[] inputs = new JComponent[] {
				new JLabel("Nx"), nxField,
				new JLabel("Ny"), nyField,
			};
			int result = JOptionPane.showConfirmDialog(this, inputs,
				"Enter grid segmentation", JOptionPane.PLAIN_MESSAGE);
			if (result != JOptionPane.OK_OPTION) return;
			// call builder function. It could throw.
			try{
				int nx = (int)nxField.getValue();
				int ny = (int)nyField.getValue();
				gd.setBasicSegmentation(nx, ny);
				tab1.validate();
				tab1.repaint();
			} catch (Exception ee){
				handleException(ee);
			}
		} else if ("run".equals(cmd)){
			//build a dialog to set mapping method
			String[] algos = {"Direct Laplace", "Inverse Laplace"};
			JComboBox<String> algoField = new JComboBox<String>(algos);
			final JComponent[] inputs = new JComponent[] {
				new JLabel("Algorithm"), algoField
			};
			if (JOptionPane.showConfirmDialog(null, inputs, "Map grid options",
				JOptionPane.PLAIN_MESSAGE) != JOptionPane.OK_OPTION) return;

			String ta = (String)algoField.getSelectedItem();
			final String algo = (algos[0].equals(ta)) ? "direct_laplace"
			                                          : "inverse_laplace";
			//call building method with progress bar popup window
			try{
				ProgressBarExecutor.execute(new ProgressBarExecutor.IExec(){
					public Exception run(Hybmesh.ICallback cb){
						return gd.doMapping(cb, algo);
					}
				});
				tab2.validate();
				tab2.repaint();
			} catch (Exception ee){
				handleException(ee);
			}
		}
	}

	// shows message boxes for errors which have occured in grid building routines
	void handleException(Exception e){
		try{
			throw e;
		} catch (Hybmesh.EUserInterrupt ee){
			JOptionPane.showConfirmDialog(this, "Interrupted",
					"Interrupted", JOptionPane.PLAIN_MESSAGE);
		} catch (Hybmesh.ERuntimeError ee){
			JOptionPane.showMessageDialog(this, ee.getMessage(),
					"Hybmesh runtime error", JOptionPane.ERROR_MESSAGE);
		} catch (Exception ee){
			JOptionPane.showMessageDialog(this, ee.getMessage(),
					"Error", JOptionPane.ERROR_MESSAGE);
		}
	}
}

// Transforms coordingate from screen to physical
class CoordinateTransform{
	int margin = 20;
	Rectangle canvas;
	double xmin, xmax, ymin, ymax;

	public CoordinateTransform(Rectangle canvas,
			double xmin, double ymin, double xmax, double ymax){
		this.canvas = canvas;
		this.xmin = xmin; this.xmax = xmax;
		this.ymin = ymin; this.ymax = ymax;
		//adjust bounding box to keep aspect ratio.
		double rel1 = (double)canvas.height/canvas.width;
		double Ly = rel1 * (xmax - xmin) - (ymax - ymin);
		double Lx = (ymax-ymin)/rel1 - (xmax - xmin);
		if (Ly > 0){
			this.ymin -= Ly/2; this.ymax += Ly/2;
		} else {
			this.xmin -= Lx/2; this.xmax += Lx/2;
		}
	}
	public int x(double rx){
		return (int)((canvas.width-2*margin)*(rx-xmin)/(xmax-xmin)) + margin;
	}
	public int y(double ry){
		return (int)((canvas.height-2*margin)*(ymax-ry)/(ymax-ymin)) + margin;
	}
	public double realx(int x){
		return (double)(x-margin)/(canvas.width-2*margin)*(xmax-xmin)+xmin;
	}
	public double realy(int y){
		return -((double)(y-margin)/(canvas.height-2*margin)*(ymax-ymin)-ymax);
	}
}

// Component for drawing grid, contour and vertices
class Drawer extends JPanel implements MouseListener, MouseMotionListener{
	public static final long serialVersionUID = 1L;
	GridData.Grid target;
	CoordinateTransform transf;
	App app;
	int moved_point = -1;
	int point_radius = 6;

	public Drawer(GridData.Grid target, App app){
		this.target = target;
		this.app = app;
		addMouseListener(this);
		addMouseMotionListener(this);
		setBackground(Color.WHITE);
	}

	void refreshTransf(){
		Rectangle rect = getBounds();
		transf = new CoordinateTransform(rect,
				target.minX(), target.minY(), target.maxX(), target.maxY());
	}

	void paintGrid(Graphics2D g){
		g.setColor(Color.BLACK);
		g.setStroke(new BasicStroke(1));
		for (int i=0; i<target.edge_vert.length; i+=2){
			g.drawLine(transf.x(target.vertices[2*target.edge_vert[i]]),
				transf.y(target.vertices[2*target.edge_vert[i]+1]),
				transf.x(target.vertices[2*target.edge_vert[i+1]]),
				transf.y(target.vertices[2*target.edge_vert[i+1]+1]));
		}
	}

	void paintContour(Graphics2D g){
		g.setColor(Color.BLUE);
		g.setStroke(new BasicStroke(3));
		int[] x = new int[4], y = new int[4];
		for (int i=0; i<4; ++i){
			x[i] = transf.x(target.corner_vert.get(i).x);
			y[i] = transf.y(target.corner_vert.get(i).y);
		}
		g.drawPolygon(x, y, 4);
	}

	void drawPoint(Graphics2D g, Hybmesh.Point2 p, String s){
		g.fillOval(transf.x(p.x)-point_radius, transf.y(p.y)-point_radius,
			2*point_radius, 2*point_radius);
		g.drawString(s, transf.x(p.x)+point_radius, transf.y(p.y));
	}
	void paintPoints(Graphics2D g){
		g.setColor(Color.BLACK);
		g.setStroke(new BasicStroke(3));
		ArrayList<Hybmesh.Point2> pts = target.corner_vert;
		for (int i=0; i<pts.size(); ++i){
			drawPoint(g, pts.get(i), Integer.toString(i));
		}
	}
	void paintUPoints(Graphics2D g){
		g.setColor(Color.BLUE);
		g.setStroke(new BasicStroke(3));
		for (int i=0; i<target.user_defined_vert.size(); ++i){
			drawPoint(g, target.udefToPoint(i), Integer.toString(4+i));
		}
	}

	@Override
	public void paintComponent(Graphics g){
		super.paintComponent(g);
		refreshTransf();
		//1) grid
		if (target.hm_grid != null) paintGrid((Graphics2D)g);
		//2) contour
		paintContour((Graphics2D)g);
		//3) corner points
		paintPoints((Graphics2D)g);
		//4) user defined points
		paintUPoints((Graphics2D)g);
	}
	@Override
	public void mouseClicked(MouseEvent e){ }
	@Override
	public void mouseExited(MouseEvent e){ }
	@Override
	public void mouseEntered(MouseEvent e){ }
	@Override
	public void mousePressed(MouseEvent e){
		double x = transf.realx(e.getX()), y = transf.realy(e.getY());
		if (app.statusAddRef()){
			target.addUserDefinedPoint(x, y);
			validate();
			repaint();
		} else if (app.statusMoveRef()){
			moved_point = target.closestVertex(x, y);
		} else if (app.statusRemoveRef()){
			int rem_point = target.closestVertex(x, y);
			target.removeVertex(rem_point);
			validate();
			repaint();
		}
	}
	@Override
	public void mouseReleased(MouseEvent e){
		if (app.statusMoveRef()){ moved_point = -1; }
	}
	@Override
	public void mouseMoved(MouseEvent e) { }
	@Override
	public void mouseDragged(MouseEvent e) {
		if (moved_point != -1){
			double x = transf.realx(e.getX()), y = transf.realy(e.getY());
			target.moveVertex(moved_point, x, y);
			validate();
			repaint();
		}
	}
}

// Provides static method to call grid building routing with popup progress bar
class ProgressBarExecutor{
	public static final long serialVersionUID = 1L;

	//Popup window form.
	static class Popup extends JDialog implements Hybmesh.ICallback,
	                                              PropertyChangeListener{
		public static final long serialVersionUID = 1L;
		JLabel lab1, lab2;
		JProgressBar pb1, pb2;
		boolean isCancelled = false;

		public Popup(){
			super((JFrame)null, "Hybmesh operation", true);
			setDefaultCloseOperation(DO_NOTHING_ON_CLOSE);
			lab1 = new JLabel("n1");
			lab2 = new JLabel("n2");
			pb1 = new JProgressBar(0, 100);
			pb2 = new JProgressBar(0, 100);
			JButton cancel_button = new JButton("Cancel");
			cancel_button.addActionListener(new ActionListener(){
				public void actionPerformed(ActionEvent e){
					isCancelled = true;
				}
			});
			JOptionPane pane = new JOptionPane(
				new Object[]{lab1, pb1, lab2, pb2},
				JOptionPane.PLAIN_MESSAGE,
				JOptionPane.DEFAULT_OPTION, null,
				new Object[]{cancel_button});
			setContentPane(pane);
			setSize(new Dimension(400, 200));
			setResizable(false);
			setLocationRelativeTo(null);
			addPropertyChangeListener(this);
		}
		@Override
		public void propertyChange(PropertyChangeEvent e){
			String nm = e.getPropertyName();
			if ("n1".equals(nm)) lab1.setText((String)e.getNewValue());
			else if ("n2".equals(nm)) lab2.setText((String)e.getNewValue());
			else if ("p1".equals(nm)) pb1.setValue((Integer)e.getNewValue());
			else if ("p2".equals(nm)){
				int v = (Integer)e.getNewValue();
				if (v < 0) pb2.setIndeterminate(true);
				else{
					pb2.setIndeterminate(false);
					pb2.setValue(v);
				}
			}
		}
		@Override
		public int callback(String n1, String n2, double p1, double p2){
			try{
				if (isCancelled) throw new InterruptedException();
				//Using events to change window components because
				//this function is called from a thread separated from gui.
				firePropertyChange("n1", null, n1);
				firePropertyChange("n2", null, n2);
				firePropertyChange("p1", null, new Integer((int)(100*p1)));
				firePropertyChange("p2", null, new Integer((int)(100*p2)));
				//!!!! Intentinal delay for testing reasons.
				Thread.sleep(300);
			} catch (InterruptedException e){
				return 1;
			}
			return 0;
		}
	};

	// Backgroung worker which runs building routine in a separate thread.
	static class Worker extends SwingWorker<Exception, Void>{
		IExec func;
		Popup cb;
		Worker(IExec func, Popup cb){
			this.func = func;
			this.cb = cb;
			execute();
			cb.setVisible(true);
		}
		@Override
		public Exception doInBackground(){ return func.run(cb); }
		@Override
		public void done(){ cb.setVisible(false); }
	}

	// Interfaces for grid building routine
	public static interface IExec{
		public Exception run(Hybmesh.ICallback cb);
	}

	// Main static function. Rethrows all grid building exception.
	// If process was cancelled throws Hybmesh.EUserInterrupt.
	public static void execute(IExec func) throws Exception{
		Worker worker = new Worker(func, new Popup());
		if(worker.get() != null) throw worker.get();
	}
};
