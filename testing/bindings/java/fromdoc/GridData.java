import java.util.ArrayList;
import static java.lang.Math.*;

// Stores left and right grids and defined points
public class GridData{
	Hybmesh hm;
	Grid from, to;

	public GridData() throws Exception{
		// Establish hybmesh connection
		hm = new Hybmesh();
		// build base grid with default partition
		from = Grid.sqrGrid(hm, 10, 10);
		to = Grid.sqrContour();
		from.allow_contour_move = false;
		to.allow_contour_move = true;
	}

	// Sets basic grid segmentation
	public void setBasicSegmentation(int nx, int ny) throws Exception{
		Grid tmp = Grid.sqrGrid(hm, nx, ny);
		from.hm_grid = tmp.hm_grid;
		from.vertices = tmp.vertices;
		from.edge_vert = tmp.edge_vert;
	}

	// Performs mapping. Signature fits ProgressBarExecutor.IExec interface.
	public void doMapping(Hybmesh.ICallback cb, String algo) throws Exception{
		hm.assignCallback(cb);
		try{
			//1. create target contour
			Hybmesh.Contour2D target_contour = to.mainContour(hm);
			//2. assemble base and target points
			Hybmesh.Point2[] base_points = from.allVertices();
			Hybmesh.Point2[] target_points = to.allVertices();
			//3. do mapping
			Hybmesh.Grid2D ret = hm.mapGrid(
				from.hm_grid, target_contour, base_points, target_points,
				null, null, algo, false, false);
			//4. add to result
			to.fillGrid(ret);
		} finally {
			hm.resetCallback();
		}
	}

	//Grid data storage
	public static class Grid{
		double[] vertices;
		int[] edge_vert;
		Hybmesh.Grid2D hm_grid=null;
		ArrayList<Hybmesh.Point2> corner_vert;
		ArrayList<Double> user_defined_vert;
		boolean allow_contour_move = false;

		public void fillGrid(Hybmesh.Grid2D hm_grid){
			try{
				this.hm_grid = hm_grid;
				vertices = hm_grid.rawVertices();
				edge_vert = hm_grid.rawTab("edge_vert");
			} catch (Exception e){
				this.hm_grid = null; vertices = null; edge_vert = null;
			}
		}
		
		// Bounding box coordinates
		public double minX(){
			return min(min(corner_vert.get(0).x, corner_vert.get(1).x),
				min(corner_vert.get(2).x, corner_vert.get(3).x));
		}
		public double minY(){
			return min(min(corner_vert.get(0).y, corner_vert.get(1).y),
				min(corner_vert.get(2).y, corner_vert.get(3).y));
		}
		public double maxX(){
			return max(max(corner_vert.get(0).x, corner_vert.get(1).x),
				max(corner_vert.get(2).x, corner_vert.get(3).x));
		}
		public double maxY(){
			return max(max(corner_vert.get(0).y, corner_vert.get(1).y),
				max(corner_vert.get(2).y, corner_vert.get(3).y));
		}

		// corner + user-defined vertices
		public Hybmesh.Point2[] allVertices(){
			Hybmesh.Point2[] ret = new Hybmesh.Point2[4 + user_defined_vert.size()];
			for (int i=0; i<4; ++i) ret[i] = corner_vert.get(i);
			for (int i=0; i<user_defined_vert.size(); ++i){
				ret[i+4] = udefToPoint(i);
			}
			return ret;
		}
		// creates a hybmesh contour object from four corner points
		public Hybmesh.Contour2D mainContour(Hybmesh hm){
			try{
				Hybmesh.Point2[] ret = new Hybmesh.Point2[5];
				for (int i=0; i<4; ++i) ret[i] = corner_vert.get(i);
				ret[4] = ret[0];
				return hm.createContour(ret, null);
			} catch (Exception e){
				return null;
			}
		}

		// indexed user-defined point to physical point
		Hybmesh.Point2 udefToPoint(int index){
			ArrayList<Hybmesh.Point2> pts = corner_vert;
			double w = user_defined_vert.get(index);
			int i1, i2;
			if (w > 3){ w-=3; i1 = 3; i2 = 0; }
			else if (w > 2) { w-=2; i1 = 2; i2 = 3;}
			else if (w > 1) { w-=1; i1 = 1; i2 = 2;}
			else { w-=0; i1 = 0; i2 = 1;}
			return new Hybmesh.Point2((1-w)*pts.get(i1).x + (w)*pts.get(i2).x,
					(1-w)*pts.get(i1).y + (w)*pts.get(i2).y);
		}
		
		//squared distance betwen two points
		double meas(Hybmesh.Point2 p, double x, double y){
			return (p.x - x)*(p.x - x) + (p.y - y)*(p.y - y);
		}
		// projects point to contours i-th section.
		// returns {squared distance to section, [0, 1] coordinate of the closest point}
		double[] sect_project(double x, double y, int isec){
			Hybmesh.Point2 p1 = corner_vert.get(isec);
			Hybmesh.Point2 p2 = corner_vert.get(isec == 3 ? 0 : isec + 1);
			double ax = p2.x - p1.x, ay = p2.y - p1.y;
			double bx = x - p1.x, by = y - p1.y;
			double ksi=(ax*bx+ay*by)/(ax*ax + ay*ay);
			if (ksi>=1) { return new double[] {meas(p2, x, y), 1.0}; }
			else if (ksi<=0) { return new double[]{meas(p1, x, y), 0.0}; } 
			else{
				double[] A = new double[]{
					p1.y-p2.y, p2.x-p1.x, p1.x*p2.y-p1.y*p2.x};
				double d0 = A[0]*x+A[1]*y+A[2];
				d0 *= d0; d0 /= (A[0]*A[0] + A[1]*A[1]);
				return new double[]{d0, ksi};
			}
		}
		// gives [0, 4] contour coordinate of given point where
		// [0, 1] are normalized coordinates for the first section,
		// [1, 2]           --- // ---       for the secont section, etc.
		Double calcUserDefinedPoint(double x, double y){
			//find closest contour section
			double[] d1 = sect_project(x, y, 0);
			double[] d2 = sect_project(x, y, 1);
			double[] d3 = sect_project(x, y, 2);
			double[] d4 = sect_project(x, y, 3);
			if (d1[0] <= d2[0] && d1[0] <= d3[0] && d1[0] <= d4[0]){
				return new Double(d1[1]);
			} else if (d2[0] <= d1[0] && d2[0] <= d3[0] && d2[0] <= d4[0]){
				return new Double(d2[1]+1.0);
			} else if (d3[0] <= d1[0] && d3[0] <= d2[0] && d3[0] <= d4[0]){
				return new Double(d3[1]+2.0);
			} else{
				return new Double(d4[1]+3.0);
			}
		}
		// adds a point to the user defined list
		public void addUserDefinedPoint(double x, double y){
			user_defined_vert.add(calcUserDefinedPoint(x, y));
		}

		// finds index of closest given vertex.
		// [0, 1, 2, 3] - are indicies of corner points, all other - user defined ones.
		public int closestVertex(double x, double y){
			int ret = -1;
			double meas_min = 1e100;
			//corner vertices
			for (int i=0; i<4; ++i){
				double m = meas(corner_vert.get(i), x, y);
				if (m <= meas_min){ meas_min = m; ret = i; }
			}
			//user vertices
			for (int i=0; i<user_defined_vert.size(); ++i){
				double m = meas(udefToPoint(i), x, y);
				if (m <= meas_min){ meas_min = m; ret = i+4; }
			}
			return ret;
		}

		// moves indexed vertex to x, y location. 
		// Corner points are moved only if allow_contour_move is true.
		public void moveVertex(int ivert, double x, double y){
			if (ivert < 4){
				//corner point
				if (allow_contour_move){
					corner_vert.get(ivert).x = x;
					corner_vert.get(ivert).y = y;
				}
			} else {
				user_defined_vert.set(ivert-4, calcUserDefinedPoint(x, y));
			}
		}
		//removes indexed points. Corner points (ivert < 4) can not be removed
		public void removeVertex(int ivert){
			if (ivert >= 4){
				user_defined_vert.remove(ivert-4);
			}
		}

		// initilaizes Grid object with zero grid and 4 corner points
		public static Grid sqrContour(){
			Grid ret = new Grid();
			ret.user_defined_vert = new ArrayList<Double>();
			ret.corner_vert = new ArrayList<Hybmesh.Point2>();
			ret.corner_vert.add(new Hybmesh.Point2(0, 0));
			ret.corner_vert.add(new Hybmesh.Point2(1, 0));
			ret.corner_vert.add(new Hybmesh.Point2(1, 1));
			ret.corner_vert.add(new Hybmesh.Point2(0, 1));
			return ret;
		}

		// initializes Grid object with unity grid and 4 corner points
		public static Grid sqrGrid(Hybmesh hm, int nx, int ny)
				throws Hybmesh.EUserInterrupt, Hybmesh.ERuntimeError{
			Grid ret = sqrContour();
			ret.fillGrid(hm.addUnfRectGrid(
				ret.corner_vert.get(0), ret.corner_vert.get(2), nx, ny, null));
			return ret;
		}
	};
};
