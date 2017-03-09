using System;
using System.Windows.Forms;
using System.Drawing;
using System.Drawing.Drawing2D;
using System.ComponentModel;
using CBFunc = System.Func<System.String, System.String,
	System.Double, System.Double, System.Int32>;
using ExecutorFunc = System.Func<
	System.Func<System.String, System.String, System.Double, System.Double, System.Int32>,
	System.Exception>;
using P2 = Hybmesh.Point2;

class App{
	/// <summary> Entry point </summary>
	public static void Main(){
		//implicitly set needed pathes
		Hybmesh.hybmesh_exec_path = "../../../src/py/hybmesh.py";
		Hybmesh.hybmesh_lib_path = "../../../build/bin";
		//create grid storage
		Builder builder = new Builder();
		//gui loop
		Application.Run(new AppGui.MainForm(builder));
	}
}

/// <summary> Holds properties of regular rectangular grid </summary>
class GridProps{
	double x0, y0, x1, y1;
	int nx, ny;
	public GridProps(int nx, int ny, double x0, double y0, double x1, double y1){
		this.nx = nx; this.ny = ny;
		this.x0 = x0; this.y0 = y0; this.x1 = x1; this.y1 = y1;
	}
	public double XMin{ get{ return x0; } set{ x0 = value;} }
	public double XMax{ get{ return x1; } set{ x1 = value;} }
	public double YMin{ get{ return y0; } set{ y0 = value;} }
	public double YMax{ get{ return y1; } set{ y1 = value;} }
	public int Nx{ get{ return nx; } set{ nx = value;} }
	public int Ny{ get{ return ny; } set{ ny = value;} }
}

/// <summary> Holds properties of unite operation </summary>
class UniteProps{
	double buffer;
	bool fix_bnd;
	public UniteProps(double buffer, bool fix_bnd){
		this.buffer = buffer; this.fix_bnd = fix_bnd;
	}
	public double Buffer{ get{ return buffer; } set{ buffer=value;} }
	public bool FixBoundary{ get{ return fix_bnd; } set{ fix_bnd=value;} }
}

/// <summary> Storage of grid data used for drawing </summary>
class Grid{
	public double[] xy;                //Coordinates in a raw array
	public int[] edge_vert;            //Edge->Vertex connectivity in a raw array
	public int[] bnd;                  //(bnd edge, bnd type) pairs in a raw array
	public Hybmesh.Grid2D hm;          //handle used in UniteGrids operation
	double _xmin, _xmax, _ymin, _ymax; //Bounding Box

	Grid(Hybmesh.Grid2D descriptor){
		hm = descriptor;
		//get data from hybmesh descriptor
		xy = descriptor.RawVertices();
		edge_vert = descriptor.RawTab("edge_vert");
		bnd = descriptor.RawTab("bnd_bt");
		//boundary box
		_xmin = xy[0]; _xmax = xy[0]; _ymin = xy[1]; _ymax = xy[1];
		for (int i=0; i<xy.Length; ++i){
			if (i % 2 == 0 && xy[i]<_xmin) _xmin = xy[i];
			if (i % 2 == 0 && xy[i]>_xmax) _xmax = xy[i];
			if (i % 2 == 1 && xy[i]<_ymin) _ymin = xy[i];
			if (i % 2 == 1 && xy[i]>_ymax) _ymax = xy[i];
		}
	}
	/// <summary> Builds a grid from descriptor or returns null </summary>
	public static Grid Build(Hybmesh.Grid2D descriptor){
		return (descriptor != null) ? new Grid(descriptor) : null;
	}

	/// <summary> fills bounding box of the grid </summary>
	public void BBox(out double xmin, out double ymin, out double xmax, out double ymax){
		xmin = _xmin; ymin = _ymin; xmax = _xmax; ymax = _ymax;
	}
}

/// <summary> Holds hybmesh handle, grids and provides methods for grid manipulations </summary>
class Builder{
	public Hybmesh hm = new Hybmesh();
	public Grid g1, g2, result;
	public GridProps props_grid1 = new GridProps(10, 10, 0, 0, 1, 1);
	public GridProps props_grid2 = new GridProps(10, 10, 0.3, 0.3, 0.5, 0.5);
	public UniteProps props_unite = new UniteProps(0.1, false);

	Grid BuildGridB(GridProps props, int bnd){
		return Grid.Build(hm.AddUnfRectGrid(
			new P2(props.XMin, props.YMin), new P2(props.XMax, props.YMax),
			props.Nx, props.Ny, new int[]{bnd, bnd, bnd, bnd}));
	}

	// Fills this.g1/this.g2 grids or throws;
	public void BuildGrid1Operation(){
		g1 = BuildGridB(props_grid1, 1);
	}
	public void BuildGrid2Operation(){
		g2 = BuildGridB(props_grid2, 2);
	}

	// Fills this.result by union of g1 and g2. Does not throw but returns exception if any.
	// Signature fits ProgressBarExecutor.Exec requirements for visual callback.
	public Exception UniteGridsOperation(CBFunc cb){
		try{
			if (g1 == null || g2 == null) throw new Exception(
					"Source grids have not been built");
			hm.AssignCallback(cb.Invoke);
			result = Grid.Build(hm.UniteGrids1(
				g1.hm, g2.hm,
				props_unite.Buffer, false, props_unite.FixBoundary));
			return null;
		} catch (Exception e){
			return e;
		} finally {
			hm.ResetCallback();
		}
	}

	/// <summary> bounding box with respect to all existing grids </summary>
	public void BBox(out double xmin, out double ymin, out double xmax, out double ymax){
		xmin=ymin=xmax=ymax=0;
		//calculate
		if (g1 != null){
			g1.BBox(out xmin, out ymin, out xmax, out ymax);
			if (g2 != null){
				double x1, y1, x2, y2;
				g2.BBox(out x1, out y1, out x2, out y2);
				if (xmin > x1) xmin = x1; if (ymin > y1) ymin = y1;
				if (xmax < x2) xmax = x2; if (ymax < y2) ymax = x2;
			}
		} else if (g2 != null){
			g2.BBox(out xmin, out ymin, out xmax, out ymax);
		}
	}
};

/// <summary>
/// Provides method to execute Action with popup cancel dialog with a pair of progress bars.
/// </summary>
class ProgressBarExecutor{

	/// <summary> Needed to execute worker in background and do not hang gui </summary>
	class BGWorker: BackgroundWorker{
		ExecutorFunc worker;
		Action finalizer;
		CBFunc reporter;
		public Exception result;

		public BGWorker(ExecutorFunc worker, CBFunc reporter, Action finalizer){
			this.worker = worker;
			this.reporter = reporter;
			this.finalizer = finalizer;
			WorkerSupportsCancellation = true;
			DoWork += EventHandler;
		}
		void EventHandler(object sender, DoWorkEventArgs e){
			e.Result = worker(reporter);
		}
		protected override void OnRunWorkerCompleted(RunWorkerCompletedEventArgs e){
			result = e.Result as Exception;
			finalizer();
		}
	};

	///<summary> Popup dialog form: progress bar + cancel button </summary>
	class CBForm: Form{
		ProgressBar pb1 = new ProgressBar(), pb2 = new ProgressBar();
		TableLayoutPanel layout = new TableLayoutPanel();
		Label lab1 = new Label(), lab2 = new Label();
		Button cancel = new Button();
		BGWorker bg;
		bool stopped = false;
		public Exception ExecError { get { return bg.result; } }

		public CBForm(ExecutorFunc worker){
			Text = "Hybmesh operation";
			cancel.Anchor = lab1.Anchor = lab2.Anchor = 0;
			lab1.TextAlign = lab2.TextAlign = ContentAlignment.MiddleCenter;
			lab1.Width = lab2.Width = pb1.Width = pb2.Width = 400;
			pb1.Maximum = pb2.Maximum = 100;
			cancel.Text = "Cancel";
			cancel.Click += OnCancel;
			layout.AutoSize = true;
			layout.ColumnCount = 1;
			layout.RowCount = 5;
			layout.Controls.Add(lab1, 0, 0);
			layout.Controls.Add(pb1, 0, 1);
			layout.Controls.Add(lab2, 0, 2);
			layout.Controls.Add(pb2, 0, 3);
			layout.Controls.Add(cancel, 0, 4);
			layout.Dock = System.Windows.Forms.DockStyle.Fill;
			AutoSize = true;
			AutoSizeMode = AutoSizeMode.GrowAndShrink;
			Controls.Add(layout);
			FormBorderStyle = FormBorderStyle.FixedDialog;
			CenterToScreen();
			bg = new BGWorker(worker, Report, Close);
			bg.RunWorkerAsync();  //Start operation in background
		}
		protected override void OnFormClosing(FormClosingEventArgs e){
			stopped = true; //to stop on forced form close.
		}
		void OnCancel(object sender, EventArgs e){
			stopped = true; //to stop on cancel click.
		}

		//function passed as CBFunc to worker.
		int Report(string n1, string n2, double p1, double p2){
			lab1.Text = n1; lab2.Text = n2;
			pb1.Value = (p1 < 0) ? 0 : (int)(p1*100);
			pb2.Value = (p2 < 0) ? 0 : (int)(p2*100);
			//!!!!! Intentional delay for testing reasons.
			System.Threading.Thread.Sleep(500);
			return stopped ? 1 : 0;
		}
	};

	/// <summary>
	/// Executes worker with popup progress bar form.
	/// Rethrows all exceptions. Throws Hybmesh.EUserInterrupt on cancellation request.
	/// </summary>
	static public void Exec(ExecutorFunc worker){
		using (var dialog = new CBForm(worker)){
			dialog.ShowDialog();
			if (dialog.ExecError != null) throw dialog.ExecError;
		}
	}
}

namespace AppGui{
	
/// <summary> Panel with automatic size and position adjustment. </summary>
class AppPanel: Panel{
	double xpos, ypos, xsz, ysz;
	/// <param name = xpos> Horizontal location normalized to [0, 1]</param> 
	/// <param name = ypos> Vertical location normalized to [0, 1]</param> 
	/// <param name = xsz>  Panel width normalized to [0, 1] </param> 
	/// <param name = ysz>  Panel height normalized to [0, 1] </param> 
	public AppPanel(Control parent,
			double xpos, double ypos, double xsz, double ysz): base(){
		this.xpos = xpos; this.ypos = ypos; this.xsz = xsz; this.ysz = ysz;
		this.BorderStyle = BorderStyle.Fixed3D;
		parent.Controls.Add(this);
		parent.Resize += AutoResize;
	}

	//protected override void OnResize(EventArgs e){
	void AutoResize(object sender, EventArgs e){
		int x = Parent.ClientSize.Width;
		int y = Parent.ClientSize.Height;
		Location = new Point((int)(x*xpos), (int)(y*ypos));
		Size = new Size((int)(x*xsz), (int)(y*ysz));
	}
}

/// <summary> Options section: caption + properties + Apply button </summay>
class OptDisplay: AppPanel{
	Button OkButton = new Button();
	PropertyGrid props = new PropertyGrid();
	public CheckBox cap = new CheckBox();

	/// <param name="Props"> Property structure for PropertyGrid widget</param>
	public OptDisplay(string caption, object Props,
			MainForm parent, double xpos, double ypos, double xsz, double ysz)
				: base(parent, xpos, ypos, xsz, ysz){
		cap.Text = caption;
		cap.Checked = true;
		cap.CheckedChanged += parent.RedrawCanvasHandler;
		OkButton.Text = "Apply";
		OkButton.Click += parent.ApplyHandler;
		OkButton.Click += parent.RedrawCanvasHandler;
		props.ToolbarVisible =  props.HelpVisible = false;
		props.PropertySort=PropertySort.NoSort;
		props.Location = new Point(0, cap.Size.Height);
		props.SelectedObject = Props;
		Controls.Add(props);
		Controls.Add(OkButton);
		Controls.Add(cap);
	}

	protected override void OnResize(EventArgs e){
		var s = ClientSize;
		props.Size = new Size(s.Width, s.Height-OkButton.Size.Height-cap.Size.Height);
		OkButton.Location = new Point(s.Width-OkButton.Size.Width,
					      props.Location.Y+props.Size.Height);
	}

	///<summary> checks whether sender is this.ApplyButton </summary>
	public bool IsSender(object sender){ return Object.ReferenceEquals(sender, OkButton); }
}

/// <summary> Grid draw canvas </summary>
class Drawer: AppPanel{
	//phisical bounding box of g1 and g2
	double xmin, ymin, xmax, ymax;
	int margin = 5;
	int CWidth{ get{ return ClientSize.Width - 2*margin; }}
	int CHeight{ get{ return ClientSize.Height - 2*margin; }}
	void AdjustBBox(){
		((MainForm)Parent).builder.BBox(out xmin, out ymin, out xmax, out ymax);
		//adjust bounding box to form sizes so that Lx/Ly == Width/Heigh
		//to keep aspect ratio.
		double rel1 = (double)Size.Height/Size.Width;
		double Ly = rel1 * (xmax - xmin) - (ymax - ymin);
		double Lx = (ymax-ymin)/rel1 - (xmax - xmin);
		if (Ly > 0){
			ymin -= Ly/2; ymax += Ly/2;
		} else {
			xmin -= Lx/2; xmax += Lx/2;
		}
	}
	//physical coordinates => canvas coordinates
	int XToCanvas(double x){ return (int)(CWidth*(x-xmin)/(xmax-xmin)) + margin; }
	int YToCanvas(double y){ return (int)(CHeight*(ymax-y)/(ymax-ymin)) + margin; }

	public Drawer(Control parent, double xpos, double ypos, double xsz, double ysz)
			: base(parent, xpos, ypos, xsz, ysz){
		BackColor = Color.White;
	}

	//Draw grid edges. This is called from OnPaint(e).
	void DrawGrid(Grid g, PaintEventArgs e){
		Graphics X = e.Graphics;
		Pen pen0 = new Pen(Color.Black, 1); //internal edges
		Pen pen1 = new Pen(Color.Red, 2);   //bnd=1 edge pen
		Pen pen2 = new Pen(Color.Blue, 2);  //bnd=2 edge pen
		//internal edges
		for (int i=0; i<g.edge_vert.Length/2; ++i){
			int v1 = g.edge_vert[2*i], v2 = g.edge_vert[2*i+1];
			X.DrawLine(pen0, XToCanvas(g.xy[2*v1]), YToCanvas(g.xy[2*v1+1]),
				XToCanvas(g.xy[2*v2]), YToCanvas(g.xy[2*v2+1]));
		}
		//boundary edges
		for (int i=0; i<g.bnd.Length/2; ++i){
			int v1 = g.edge_vert[2*g.bnd[2*i]], v2 = g.edge_vert[2*g.bnd[2*i]+1];
			Pen peni = (g.bnd[2*i+1] == 1 ? pen1 : pen2);
			X.DrawLine(peni, XToCanvas(g.xy[2*v1]), YToCanvas(g.xy[2*v1+1]),
				XToCanvas(g.xy[2*v2]), YToCanvas(g.xy[2*v2+1]));
		}
	}

	protected override void OnPaint(PaintEventArgs e){
		AppGui.MainForm b = (AppGui.MainForm)Parent;
		AdjustBBox();
		if (b.VisibleGrid1) DrawGrid(b.builder.g1, e);
		if (b.VisibleGrid2) DrawGrid(b.builder.g2, e);
		if (b.VisibleResult) DrawGrid(b.builder.result, e);
	}
};

class MainForm: Form{
	public Builder builder;
	Drawer canvas;
	OptDisplay gui_grid1, gui_grid2, gui_unite;
	public bool VisibleGrid1{get{return gui_grid1.cap.Checked && builder.g1 != null;}}
	public bool VisibleGrid2{get{return gui_grid2.cap.Checked && builder.g2 != null;}}
	public bool VisibleResult{get{return gui_unite.cap.Checked && builder.result != null;}}

	public MainForm(Builder builder){
		this.builder = builder;
		canvas = new Drawer(this, 0, 0, 0.7, 1.0);
		gui_grid1 = new OptDisplay("Basic grid", builder.props_grid1,
				this, 0.7, 0, 0.3, 0.4);
		gui_grid2 = new OptDisplay("Secondary grid", builder.props_grid2,
				this, 0.7, 0.4, 0.3, 0.4);
		gui_unite = new OptDisplay("Unite", builder.props_unite,
					this, 0.7, 0.8, 0.3, 0.2);
		Text = "App";
		Size = new Size(700, 500);
		CenterToScreen();
	}

	/// <summary> Called by apply click. Launches builders, shows error messages </summary>
	public void ApplyHandler(object sender, EventArgs e){
		try{
			// Choose operation by checking who called it.
			if (gui_grid1.IsSender(sender)){
				builder.BuildGrid1Operation();
			} else if (gui_grid2.IsSender(sender)){
				builder.BuildGrid2Operation();
			} else if (gui_unite.IsSender(sender)){
				ProgressBarExecutor.Exec(builder.UniteGridsOperation);
			}
		} catch (Hybmesh.EUserInterrupt){
			MessageBox.Show("Interrupted", "", MessageBoxButtons.OK,
					MessageBoxIcon.Information);
		} catch (Hybmesh.ERuntimeError ee){
			MessageBox.Show(ee.Message, "Hybmesh error", MessageBoxButtons.OK,
					MessageBoxIcon.Error);
		} catch (Exception ee){
			MessageBox.Show(ee.Message, "Applicaton error", MessageBoxButtons.OK,
					MessageBoxIcon.Error);
		}
	}
	public void RedrawCanvasHandler(object sender, EventArgs e){
		canvas.Invalidate();
	}
}
}
