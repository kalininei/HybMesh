import java.util.Map;

class GoodCallback implements Hybmesh.ICallback{
	public int callback(String s1, String s2, double p1, double p2){
		System.out.println(s1+" "+Double.toString(p1)+
			" --- "+s2+" "+Double.toString(p2));
		return 0;
	}
}
class BadCallback implements Hybmesh.ICallback{
	public int callback(String s1, String s2, double p1, double p2){
		System.out.println(s1+" "+Double.toString(p1)+
			" --- "+s2+" "+Double.toString(p2));
		return (p1<0.7) ? 0 : 1;
	}
}

class P2 extends Hybmesh.Point2{
	public P2(double x, double y){
		super(x, y);
	}
}
class P3 extends Hybmesh.Point3{
	public P3(double x, double y, double z){
		super(x, y, z);
	}
}

class a{

static void checkCond(boolean cond) throws Exception{
	if (!cond) throw new Exception("check failed");
}

static void checkdims(int[] v1, int[] v2) throws Exception{
	if (v1.length != v2.length) throw new Exception("check failed");
	for (int i=0; i<v2.length; ++i){
		if (v1[i] != v2[i]) throw new Exception("check failed");
	}
}

public static void pings(Hybmesh hm) throws Exception{
	Hybmesh.Grid2D g1 = hm.addUnfRectGrid1(
		new double[]{1, 2, 3}, new double[]{2, 3, 4}, null);
	checkdims(g1.dims(), new int[]{9, 12, 4});
	Hybmesh.Contour2D c1 = hm.gridBndToContour(g1, true);
	Hybmesh.Contour2D c2 = hm.gridBndToContour(g1, false);
	checkdims(c1.dims(), new int[]{4, 4});
	checkdims(c2.dims(), new int[]{8, 8});
	P2[] pv = new P2[4];
	pv[0] = new P2(0, 0);
	pv[1] = new P2(1, 2);
	pv[2] = new P2(3, 1);
	pv[3] = new P2(2, 0);
	Hybmesh.Contour2D c3 = hm.createContour(pv, new int[]{0, 0, 1});
	Hybmesh.Contour2D c4 = hm.createContour(pv, new int[]{3, 2, 1});
	checkdims(c3.dims(), new int[]{4, 3});
	checkdims(c4.dims(), new int[]{4, 3});
	Hybmesh.Contour2D c5 = hm.createSplineContour(pv, new int[]{0, 0, 1}, 100);
	Hybmesh.Contour2D c6 = hm.createSplineContour(pv, new int[]{3, 2, 1}, 100);
	checkdims(c5.dims(), new int[]{101, 100});
	checkdims(c6.dims(), new int[]{101, 100});
	Hybmesh.Contour2D[] c7 = hm.extractSubcontours(c6, pv);
	checkCond(c7.length == 3);
	Hybmesh.Contour2D c8 = hm.addRectContour(new P2(1, 1), new P2(10, 10), new int[]{0, 1, 2, 3});
	checkdims(c8.dims(), new int[]{4, 4});
	Hybmesh.Contour2D c9 = hm.addCircContour(new P2(5.0, 5.0), 10, 32, 0);
	checkdims(c9.dims(), new int[]{32, 32});
	double ar9 = c9.domainArea();
	checkCond(Math.abs(ar9 - 3.1415926*100) < 5);
	Hybmesh.Contour2D c10 = hm.simplifyContour(c2, 30);
	checkdims(c10.dims(), new int[]{4, 4});
	Hybmesh.Contour2D c11 = hm.simplifyContour(c9, 91);
	checkdims(c11.dims(), new int[]{4, 4});
	Hybmesh.Contour2D c12 = hm.uniteContours(new Hybmesh.Contour2D[]{c8, c9});
	checkdims(c12.dims(), new int[]{36, 36});
	Hybmesh.Contour2D[] c13 = hm.decomposeContour(c12);
	checkCond(c13.length == 2);
	Hybmesh.Contour2D c14 = hm.clipDomain(c13[0], c13[1], "difference", true);
	Hybmesh.Contour2D c15 = hm.addRectContour(new P2(0, 0), new P2(1, 1), null);
	Hybmesh.Contour2D c16 = hm.addRectContour(new P2(10, 10), new P2(11, 11), null);
	Hybmesh.Contour2D c17 = hm.clipDomain(c15, c16, "intersection", true);
	Hybmesh.Contour2D c18 = hm.clipDomain(c13[1], c13[1], "difference", true);
	checkdims(c17.dims(), new int[]{0, 0});
	checkdims(c18.dims(), new int[]{0, 0});
	Hybmesh.Contour2D c19 = hm.partitionContourConst(c15, 0.01, 30, false, -1, null, null, null, null);
	checkdims(c19.dims(), new int[]{400, 400});
	Hybmesh.Contour2D c20 = hm.partitionContourConst(c15, 0.01, 30, false, -1, null, null, null, null);
	checkdims(c19.dims(), new int[]{400, 400});
	Hybmesh.Contour2D c21 = hm.partitionContourRefPoints(c15, new double[]{0.01, 0.5},
		new P2[]{new P2(0, 0), new P2(1, 1)},
		30, true, -1, null, null, new P2(0, 0), new P2(1, 1));
	checkdims(c21.dims(), new int[]{18, 18});
	Hybmesh.Contour2D c22 = hm.partitionContourRefLengths(c15, new double[]{0.01, 0.5},
		new double[]{0, 2},
		30, true, -1, null, null, new P2(0, 0), new P2(1, 1));
	checkdims(c21.dims(), c22.dims());
	Hybmesh.Contour2D c23 = hm.matchedPartition(c15, 0.1, 0.5, null,
			new double[]{0.01}, new P2[]{new P2(0.1, 0.1)}, 30, 3);
	checkdims(c23.dims(), new int[]{74, 74});
	double[] c24 = hm.partitionSegment(0.5, 1.5, 0.1, 0.5, new double[]{1, 0.01});
	checkCond(c24.length == 17 && c24[0] == 0.5);

	Hybmesh.Contour2D c25 = hm.createContour(new P2[]{new P2(0, 0), new P2(1, 0), new P2(2, 0)}, null);
	Hybmesh.Contour2D c26 = hm.createContour(new P2[]{new P2(2, 0.1), new P2(3, 0.1), new P2(4, 1)}, null);
	Hybmesh.Contour2D c27 = hm.connectSubcontours(new Hybmesh.Contour2D[]{c25, c26}, new int[]{}, "no", true);
	Hybmesh.Contour2D c28 = hm.connectSubcontours(new Hybmesh.Contour2D[]{c25, c26}, new int[]{}, "yes", true);
	Hybmesh.Contour2D c29 = hm.connectSubcontours(new Hybmesh.Contour2D[]{c25, c26}, new int[]{}, "force", true);
	checkdims(c27.dims(), new int[]{5, 4});
	checkdims(c28.dims(), new int[]{4, 4});
	checkdims(c29.dims(), new int[]{5, 5});

	Hybmesh.Grid2D c30 = hm.addUnfRectGrid1(new double[]{1, 2, 3}, new double[]{4, 5, 6, 7}, null);
	checkdims(c30.dims(), new int[]{12, 17, 6});
	Hybmesh.Grid2D c31 = hm.addUnfCircGrid(new P2(0, 0), 10, 32, 5, 1.2, false, 0);
	checkdims(c31.dims(), new int[]{160, 288, 129});
	Hybmesh.Grid2D c32 = hm.addUnfRingGrid(new P2(1, 1), 3, 7, 16, 2, 1.0, null);
	checkdims(c32.dims(), new int[]{48, 80, 32});
	Hybmesh.Grid2D c33 = hm.addUnfHexGridInHex(new P2(-1, -1), 8, 1, false);
	Hybmesh.Grid2D c34 = hm.addUnfHexGridInRect(new P2(0, 0), new P2(2, 3), 1, false);
	Hybmesh.Grid2D c35 = hm.addTriangleGrid(new P2(0, 0), new P2(1, 0), new P2(0, 1), 3, new int[]{1, 2, 3});
	checkdims(c33.dims(), new int[]{216, 306, 91});
	checkdims(c34.dims(), new int[]{28, 35, 8}); 
	checkdims(c35.dims(), new int[]{10, 15, 6}); 

	Hybmesh.Contour2D c36 = hm.createContour(new P2[]{new P2(0, 0), new P2(0, 1)}, null);
	Hybmesh.Contour2D c37 = hm.partitionContourConst(c36, 0.1, 30, false, -1, null, null, null, null);
	Hybmesh.Contour2D c38 = hm.createContour(new P2[]{new P2(0, 0), new P2(2, -0.1)}, null);
	Hybmesh.Contour2D c39 = hm.partitionContourConst(c38, 0.2, 30, false, -1, null, null, null, null);
	Hybmesh.Grid2D c40 = hm.addCustomRectGrid("linear", c37, c39, null, null, false);
	Hybmesh.Grid2D c41 = hm.addCustomRectGridHtfi(c37, c39,
			null, null,
			new double[]{1, 1, 1, 0.8}, false);
	hm.assignCallback(new BadCallback());
	try{
		Hybmesh.Grid2D c42_ = hm.addCustomRectGrid("orthogonal", c37, c39, null, null, false);
	} catch (Hybmesh.EUserInterrupt e){
		System.out.println("User interrupt catched");
	}
	hm.assignCallback(new GoodCallback());
	Hybmesh.Grid2D c42 = hm.addCustomRectGrid("orthogonal", c37, c39, null, null, false);
	hm.resetCallback();
	checkdims(c40.dims(), new int[]{121, 220, 100});
	checkdims(c40.dims(), c41.dims());
	checkdims(c40.dims(), c42.dims());
	try{
		Hybmesh.Grid2D c43_ = hm.addCircRectGrid(new P2(0, 0), 1, 0.05, 1.0, 1.0, "AAAlinear");
	} catch (Hybmesh.ERuntimeError e){
		System.out.println("Runtime error catched");
		System.out.println(e.getMessage());
	}
	Hybmesh.Grid2D c43 = hm.addCircRectGrid(new P2(0, 0), 1, 0.05, 1.0, 1.0, "linear");
	Hybmesh.Grid2D c44 = hm.addCircRectGrid(new P2(0, 0), 1, 0.05, 1.0, 1.0, "orthogonal_rect");
	checkdims(c43.dims(), c44.dims());
	Hybmesh.Grid2D c45 = hm.stripe(c39, new double[]{0, 0.01, 0.02, 0.05}, "no", null);
	Hybmesh.Grid2D c46 = hm.triangulateDomain(c45, null, null, null, "3");
	Hybmesh.Contour2D c47 = hm.gridBndToContour(c45, true);
	Hybmesh.Contour2D c48 = hm.gridBndToContour(c46, true);
	checkCond(Math.abs(c47.domainArea() - c48.domainArea())<1e-8);
	Hybmesh.Contour2D c49 = hm.addRectContour(new P2(0, 0), new P2(1, 1), null);
	Hybmesh.Contour2D c50 = hm.partitionContourConst(c49, 0.1, 30, false, -1, null, null, null, null);
	Hybmesh.Grid2D c51 = hm.pebiFill(c50, null, new double[]{0.5}, new P2[]{new P2(0.5, 0.5)});
	Hybmesh.Grid2D c52 = hm.buildBoundaryGrid1(c51, new double[]{0, 0.01, 0.02}, "left", null, null, null);
	Hybmesh.Grid2D c53 = hm.buildBoundaryGrid1(c51, new double[]{0, 0.01, 0.02}, "right", null, null, null);
	checkdims(c52.dims(), c53.dims());
	Hybmesh.Grid2D c54 = hm.excludeContours(c53, new Hybmesh.Object2D[]{c51}, "inner");
	Hybmesh.Grid2D c55 = hm.excludeContours(c53, new Hybmesh.Object2D[]{c51}, "outer");
	checkdims(c54.dims(), c53.dims());
	checkdims(c55.dims(), new int[]{0, 0, 0});
	double[] c56x=new double[11], c57x=new double[11];
	for (int i=0; i<11; ++i){
		c56x[i] = (i/10.0);
		c57x[i] = (0.3 + (double)i/30.0);
	}
	Hybmesh.Grid2D c56 = hm.addUnfRectGrid1(c56x, c56x, null);
	Hybmesh.Grid2D c57 = hm.addUnfRectGrid1(c57x, c57x, null);
	Hybmesh.Grid2D c58 = hm.uniteGrids1(c56, c57, 0.1, false, false, 0, "3");
	hm.stdoutVerbosity(3);
	Hybmesh.Grid2D c59 = hm.mapGrid(c58, c57, new P2[]{new P2(0, 0)}, new P2[]{new P2(0.3, 0.3)}, null, null, null, false, false);
	hm.stdoutVerbosity(0);
	hm.healGrid(c59, 30, 30);
	Hybmesh.Grid2D c60 = hm.excludeContours(c58, new Hybmesh.Object2D[]{c57}, "inner");
	Hybmesh.Grid3D c61 = hm.extrudeGrid(c60, new double[]{0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3}, 0, 0);
	Hybmesh.Surface3D c62 = hm.grid3BndToSurface(c61);
	hm.assignCallback(new GoodCallback());
	Hybmesh.Grid3D c63 = hm.tetrahedralFill(new Hybmesh.Object3D[]{c62});
	hm.resetCallback();
	Hybmesh.Grid3D c64 = hm.revolveGrid(c60, new P2(0, 0), new P2(0, 1),
			new double[]{0, 45, 90, 180}, 0, 0, false);
	Hybmesh.Grid3D c65 = hm.revolveGrid(c60, new P2(0, 0), new P2(0, 1),
			new double[]{180, 270, 360}, 0, 0, false);
	hm.export3DGridVtk(c64, "c64.vtk", "");
	hm.export3DGridVtk(c65, "c65.vtk", "");
	Hybmesh.Grid3D c66 = hm.mergeGrids3(c64, c65);
	hm.export3DGridVtk(c66, "c66.vtk", null);
	Hybmesh.Grid3D c67 = c66.deepcopy();
	c67.scale(20, 20, 20, new P3(0, 0, 0));
	Hybmesh.Surface3D c68 = hm.grid3BndToSurface(c66);
	Hybmesh.Surface3D c69 = hm.grid3BndToSurface(c67);
	checkCond(Math.abs(c68.domainVolume()-125*c69.domainVolume())<1e-8);
	c60.rotate(20, new P2(0, 0));
	c69.free();

	hm.addBoundaryType(1, "bnd1");
	hm.addBoundaryType(2, "bnd2");
	Hybmesh.Grid2D c70 = hm.addUnfRectGrid(new P2(0, 0), new P2(1, 1), 10, 10, new int[]{0, 1, 2, 3});
	hm.exportGridMsh(c70, "out.msh", new int[]{0}, new int[]{2}, new boolean[]{false});

	Map<Integer, Double> c71 = c70.skewness(-1);
	Map<Integer, Double> c72 = c70.skewness(0.1);
	checkCond(c71.size() == 100 && c72.size() == 0);

	c70.setBtypesAll(22);
	hm.exportContourVtk(c70, "out.vtk");
	double[] c73 = c70.rawVertices();
	checkCond(c73.length == 242);
	int[] c74 = c70.rawTab("bnd_bt");
	checkCond(c74.length == 80);
	int[] c75 = c70.rawTab("bnd");
	checkCond(c75.length == 40);
	int[] c76 = c70.rawTab("cell_edge");
	checkCond(c76.length == 400);
	int[] c77 = c70.rawTab("edge_cell");
	checkCond(c77.length == 440);
	int[] c78 = c70.rawTab("bt");
	boolean allof78=true;
	for(int i=0; i<c78.length; ++i){
		if (c78[i]!=22 && c78[i]!=0) allof78 = false;
	}
	checkCond(c78.length == 220 && allof78);

	int[] c79 = c70.rawTab("cell_dim");
	boolean allof79=true;
	for(int i=0; i<c79.length; ++i){
		if (c79[i]!=4) allof79 = false;
	}
	checkCond(c79.length == 100 && allof79);

	int[] c80 = c70.rawTab("cell_vert");
	checkCond(c80.length == 400);
	c70.setBtypes(1, new int[]{c75[0], c75[1], c75[2]});
	hm.exportGridTecplot(c70, "out.dat");

	Hybmesh.Grid3D c81 = hm.extrudeGrid(c70, new double[]{0, 1, 2}, 0, 0);
	Hybmesh.Surface3D c82 = hm.grid3BndToSurface(c81);
	double[] c83 = c81.rawVertices();
	double[] c84 = c82.rawVertices();
	checkCond(c83.length == 121*3*3);
	checkCond(c84.length == (121*2 + 40)*3);
	int[] c85 = c81.rawTab("cell_vert");
	checkCond(c85.length == 200*8);
	int[] c86 = c82.rawTab("face_vert");
	checkCond(c86.length == (200+80)*4);

	Hybmesh.Contour2D c87 = hm.gridBndToContour(c70, true);
	double[] c88 = c87.rawVertices();
	checkCond(c88.length == 12);
	int[] c89 = c87.rawTab("bt");
	checkCond(c89.length == 6);

	Hybmesh.Contour2D c90 = hm.addRectContour(new P2(1, 1.2), new P2(2, 2), new int[]{0, 1, 2, 3});
	Hybmesh.Contour2D c91 = hm.addRectContour(new P2(0, 0), new P2(1.5, 1.5), new int[]{0, 1, 2, 3});
	Hybmesh.Contour2D c92 = hm.uniteContours(new Hybmesh.Contour2D[]{c90, c91});
	Hybmesh.Contour2D[] c93 = hm.decomposeContour(c92);
	checkCond(c93.length == 3);
	Hybmesh.Contour2D c94 = hm.pickContour(new P2(-1, -1), c93);
	c94.getPoint(new P2(-1, -1), null, null);
	c94.getPoint(null, new P2(0.5, -1), null);
	c94.getPoint(null, null, new P2(1.0, -1));

	hm.exportAllHmd("out.hmd", "ascii");
	hm.removeAll();
	hm.importGridHmg("out.hmd", "Grid2D_1", false);
	hm.import3DGridHmg("out.hmd", "Grid3D_1", false);
	hm.removeAll();
}

public static void main(String[] args) throws Exception{
	System.out.println("Hello");

	Hybmesh.hybmesh_exec_path = "../../../src/py";

	try(Hybmesh hm = new Hybmesh()){
		pings(hm);
	}
}


}
