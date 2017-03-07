using System;
using P2 = Hybmesh.Point2;
using P3 = Hybmesh.Point3;

//==================== application
#pragma warning disable 0219

public class Program{

static void CheckCond(bool cond){
	if (!cond) throw new Exception("check failed");
}

static void CheckDims(int[] v1, int[] v2){
	if (v1.Length != v2.Length) throw new Exception("check failed");
	for (int i=0; i<v2.Length; ++i){
		if (v1[i] != v2[i]) throw new Exception("check failed");
	}
}

static int good_callback(string s1, string s2, double p1, double p2){
	Console.WriteLine(s1+" "+p1.ToString()+" --- "+s2+" "+p2.ToString());
	return 0;
}

static int bad_callback(string s1, string s2, double p1, double p2){
	Console.WriteLine(s1+" "+p1.ToString()+" --- "+s2+" "+p2.ToString());
	return (p1 > 0.7) ? 1 : 0;
}

public static void pings(Hybmesh hm){
	var g1 = hm.AddUnfRectGrid1(new double[]{1, 2, 3}, new double[]{2, 3, 4});
	CheckDims(g1.Dims(), new int[]{9, 12, 4});
	var c1 = hm.GridBndToContour(g1, true);
	var c2 = hm.GridBndToContour(g1, false);
	CheckDims(c1.Dims(), new int[]{4, 4});
	CheckDims(c2.Dims(), new int[]{8, 8});
	P2[] pv = new P2[4];
	pv[0] = new P2(0, 0);
	pv[1] = new P2(1, 2);
	pv[2] = new P2(3, 1);
	pv[3] = new P2(2, 0);
	var c3 = hm.CreateContour(pv, new int[]{0, 0, 1});
	var c4 = hm.CreateContour(pv, new int[]{3, 2, 1});
	CheckDims(c3.Dims(), new int[]{4, 3});
	CheckDims(c4.Dims(), new int[]{4, 3});
	var c5 = hm.CreateSplineContour(pv, new int[]{0, 0, 1});
	var c6 = hm.CreateSplineContour(pv, new int[]{3, 2, 1});
	CheckDims(c5.Dims(), new int[]{101, 100});
	CheckDims(c6.Dims(), new int[]{101, 100});
	var c7 = hm.ExtractSubcontours(c6, pv);
	CheckCond(c7.Length == 3);
	var c8 = hm.AddRectContour(new P2(1, 1), new P2(10, 10), new int[]{0, 1, 2, 3});
	CheckDims(c8.Dims(), new int[]{4, 4});
	var c9 = hm.AddCircContour(new P2(5.0, 5.0), 10, 32);
	CheckDims(c9.Dims(), new int[]{32, 32});
	double ar9 = c9.DomainArea();
	CheckCond(Math.Abs(ar9 - 3.1415926*100) < 5);
	var c10 = hm.SimplifyContour(c2);
	CheckDims(c10.Dims(), new int[]{4, 4});
	var c11 = hm.SimplifyContour(c9, 91);
	CheckDims(c11.Dims(), new int[]{4, 4});
	var c12 = hm.UniteContours(new Hybmesh.Contour2D[]{c8, c9});
	CheckDims(c12.Dims(), new int[]{36, 36});
	var c13 = hm.DecomposeContour(c12);
	CheckCond(c13.Length == 2);
	var c14 = hm.ClipDomain(c13[0], c13[1], "difference");
	var c15 = hm.AddRectContour(new P2(0, 0), new P2(1, 1));
	var c16 = hm.AddRectContour(new P2(10, 10), new P2(11, 11));
	var c17 = hm.ClipDomain(c15, c16, "intersection");
	var c18 = hm.ClipDomain(c13[1], c13[1], "difference");
	CheckDims(c17.Dims(), new int[]{0, 0});
	CheckDims(c18.Dims(), new int[]{0, 0});
	var c19 = hm.PartitionContourConst(c15, 0.01);
	CheckDims(c19.Dims(), new int[]{400, 400});
	var c20 = hm.PartitionContourConst(c15, 0.01);
	CheckDims(c19.Dims(), new int[]{400, 400});
	var c21 = hm.PartitionContourRefPoints(c15, new double[]{0.01, 0.5},
		new P2[]{new P2(0, 0), new P2(1, 1)},
		30, true, -1, null, null, new P2(0, 0), new P2(1, 1));
	CheckDims(c21.Dims(), new int[]{18, 18});
	var c22 = hm.PartitionContourRefLengths(c15, new double[]{0.01, 0.5},
		new double[]{0, 2},
		30, true, -1, null, null, new P2(0, 0), new P2(1, 1));
	CheckDims(c21.Dims(), c22.Dims());
	var c23 = hm.MatchedPartition(c15, 0.1, 0.5, null,
			new double[]{0.01}, new P2[]{new P2(0.1, 0.1)});
	CheckDims(c23.Dims(), new int[]{74, 74});
	var c24 = hm.PartitionSegment(0.5, 1.5, 0.1, 0.5, new double[]{0.01, 1.0});
	CheckCond(c24.Length == 5 && c24[0] == 0.5);

	var c25 = hm.CreateContour(new P2[]{new P2(0, 0), new P2(1, 0), new P2(2, 0)});
	var c26 = hm.CreateContour(new P2[]{new P2(2, 0.1), new P2(3, 0.1), new P2(4, 1)});
	var c27 = hm.ConnectSubcontours(new Hybmesh.Contour2D[]{c25, c26}, new int[]{});
	var c28 = hm.ConnectSubcontours(new Hybmesh.Contour2D[]{c25, c26}, new int[]{}, "yes");
	var c29 = hm.ConnectSubcontours(new Hybmesh.Contour2D[]{c25, c26}, new int[]{}, "force");
	CheckDims(c27.Dims(), new int[]{5, 4});
	CheckDims(c28.Dims(), new int[]{4, 4});
	CheckDims(c29.Dims(), new int[]{5, 5});

	var c30 = hm.AddUnfRectGrid1(new double[]{1, 2, 3}, new double[]{4, 5, 6, 7});
	CheckDims(c30.Dims(), new int[]{12, 17, 6});
	var c31 = hm.AddUnfCircGrid(new P2(0, 0), 10, 32, 5, 1.2, false);
	CheckDims(c31.Dims(), new int[]{160, 288, 129});
	var c32 = hm.AddUnfRingGrid(new P2(1, 1), 3, 7, 16, 2);
	CheckDims(c32.Dims(), new int[]{48, 80, 32});
	var c33 = hm.AddUnfHexGridInHex(new P2(-1, -1), 8, 1);
	var c34 = hm.AddUnfHexGridInRect(new P2(0, 0), new P2(2, 3), 1);
	var c35 = hm.AddTriangleGrid(new P2(0, 0), new P2(1, 0), new P2(0, 1), 3, new int[]{1, 2, 3});
	CheckDims(c33.Dims(), new int[]{216, 306, 91});
	CheckDims(c34.Dims(), new int[]{28, 35, 8}); 
	CheckDims(c35.Dims(), new int[]{10, 15, 6}); 

	var c36 = hm.CreateContour(new P2[]{new P2(0, 0), new P2(0, 1)});
	var c37 = hm.PartitionContourConst(c36, 0.1);
	var c38 = hm.CreateContour(new P2[]{new P2(0, 0), new P2(2, -0.1)});
	var c39 = hm.PartitionContourConst(c38, 0.2);
	var c40 = hm.AddCustomRectGrid("linear", c37, c39);
	var c41 = hm.AddCustomRectGridHtfi(c37, c39,
			null, null,
			new double[]{1, 1, 1, 0.8});
	hm.AssignCallback(bad_callback);
	try{
		var c42_ = hm.AddCustomRectGrid("orthogonal", c37, c39);
	} catch (Hybmesh.EUserInterrupt e){
		Console.WriteLine("Interrupt catched");
	}
	hm.AssignCallback(good_callback);
	var c42 = hm.AddCustomRectGrid("orthogonal", c37, c39);
	hm.ResetCallback();
	CheckDims(c40.Dims(), new int[]{121, 220, 100});
	CheckDims(c40.Dims(), c41.Dims());
	CheckDims(c40.Dims(), c42.Dims());
	try{
		var c43_ = hm.AddCircRectGrid(new P2(0, 0), -1, 0.05);
	} catch (Hybmesh.ERuntimeError e){
		Console.WriteLine("Runtime error catched");
		Console.WriteLine(e.Message);
	}

	var c43 = hm.AddCircRectGrid(new P2(0, 0), 1, 0.05);
	var c44 = hm.AddCircRectGrid(new P2(0, 0), 1, 0.05, 1.0, 1.0, "orthogonal_rect");
	CheckDims(c43.Dims(), c44.Dims());
	var c45 = hm.Stripe(c39, new double[]{0, 0.01, 0.02, 0.05});
	var c46 = hm.TriangulateDomain(c45, null, null, null, "3");
	var c47 = hm.GridBndToContour(c45);
	var c48 = hm.GridBndToContour(c46);
	CheckCond(Math.Abs(c47.DomainArea() - c48.DomainArea())<1e-8);
	var c49 = hm.AddRectContour(new P2(0, 0), new P2(1, 1));
	var c50 = hm.PartitionContourConst(c49, 0.1);
	var c51 = hm.PebiFill(c50, null, new double[]{0.5}, new P2[]{new P2(0.5, 0.5)});
	var c52 = hm.SimpleBoundaryGrid(c51, new double[]{0, 0.01, 0.02});
	var c53 = hm.SimpleBoundaryGrid(c51, new double[]{0, 0.01, 0.02}, "right");
	CheckDims(c52.Dims(), c53.Dims());
	var c54 = hm.ExcludeContours(c53, new Hybmesh.Object2D[]{c51}, "inner");
	var c55 = hm.ExcludeContours(c53, new Hybmesh.Object2D[]{c51}, "outer");
	CheckDims(c54.Dims(), c53.Dims());
	CheckDims(c55.Dims(), new int[]{0, 0, 0});
	double[] c56x=new double[11], c57x=new double[11];
	for (int i=0; i<11; ++i){
		c56x[i] = (i/10.0);
		c57x[i] = (0.3 + (double)i/30.0);
	}
	var c56 = hm.AddUnfRectGrid1(c56x, c56x);
	var c57 = hm.AddUnfRectGrid1(c57x, c57x);
	var c58 = hm.UniteGrids1(c56, c57, 0.1);
	hm.StdoutVerbosity(3);
	var c59 = hm.MapGrid(c58, c57, new P2[]{new P2(0, 0)}, new P2[]{new P2(0.3, 0.3)});
	hm.StdoutVerbosity(0);
	hm.HealGrid(c59, 30, 30);
	var c60 = hm.ExcludeContours(c58, new Hybmesh.Object2D[]{c57}, "inner");
	var c61 = hm.ExtrudeGrid(c60, new double[]{0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3});
	var c62 = hm.Grid3BndToSurface(c61);
	hm.AssignCallback(good_callback);
	var c63 = hm.TetrahedralFill(new Hybmesh.Object3D[]{c62});
	hm.ResetCallback();
	var c64 = hm.RevolveGrid(c60, new P2(0, 0), new P2(0, 1),
			new double[]{0, 45, 90, 180});
	var c65 = hm.RevolveGrid(c60, new P2(0, 0), new P2(0, 1),
			new double[]{180, 270, 360});
	hm.Export3DGridVtk(c64, "c64.vtk");
	hm.Export3DGridVtk(c65, "c65.vtk");
	var c66 = hm.MergeGrids3(c64, c65);
	hm.Export3DGridVtk(c66, "c66.vtk");
	var c67 = c66.Deepcopy();
	c67.Scale(20, 20, 20, new P3(0, 0, 0));
	var c68 = hm.Grid3BndToSurface(c66);
	var c69 = hm.Grid3BndToSurface(c67);
	CheckCond(Math.Abs(c68.DomainVolume()-125*c69.DomainVolume())<1e-8);
	c60.Rotate(20, new P2(0, 0));
	c69.Free();

	hm.AddBoundaryType(1, "bnd1");
	hm.AddBoundaryType(2, "bnd2");
	var c70 = hm.AddUnfRectGrid(new P2(0, 0), new P2(1, 1), 10, 10, new int[]{0, 1, 2, 3});
	hm.ExportGridMsh(c70, "out.msh", new int[]{0}, new int[]{2}, new bool[]{false});

	var c71 = c70.Skewness(-1);
	var c72 = c70.Skewness(0.1);
	CheckCond(c71.Length == 100 && c72.Length == 0);

	c70.SetBtypesAll(22);
	hm.ExportContourVtk(c70, "out.vtk");
	var c73 = c70.RawVertices();
	CheckCond(c73.Length == 242);
	var c74 = c70.RawTab("bnd_bt");
	CheckCond(c74.Length == 80);
	var c75 = c70.RawTab("bnd");
	CheckCond(c75.Length == 40);
	var c76 = c70.RawTab("cell_edge");
	CheckCond(c76.Length == 400);
	var c77 = c70.RawTab("edge_cell");
	CheckCond(c77.Length == 440);
	var c78 = c70.RawTab("bt");
	bool allof78=true;
	foreach (int i in c78){
		if (i!=22 && i!=0) allof78 = false;
	}
	CheckCond(c78.Length == 220 && allof78);

	var c79 = c70.RawTab("cell_dim");
	bool allof79=true;
	foreach (int i in c79){
		if (i!=4) allof79 = false;
	}
	CheckCond(c79.Length == 100 && allof79);

	var c80 = c70.RawTab("cell_vert");
	CheckCond(c80.Length == 400);
	c70.SetBtypes(1, new int[]{c75[0], c75[1], c75[2]});
	hm.ExportGridTecplot(c70, "out.dat");

	var c81 = hm.ExtrudeGrid(c70, new double[]{0, 1, 2});
	var c82 = hm.Grid3BndToSurface(c81);
	var c83 = c81.RawVertices();
	var c84 = c82.RawVertices();
	CheckCond(c83.Length == 121*3*3);
	CheckCond(c84.Length == (121*2 + 40)*3);
	var c85 = c81.RawTab("cell_vert");
	CheckCond(c85.Length == 200*8);
	var c86 = c82.RawTab("face_vert");
	CheckCond(c86.Length == (200+80)*4);

	var c87 = hm.GridBndToContour(c70);
	var c88 = c87.RawVertices();
	CheckCond(c88.Length == 12);
	var c89 = c87.RawTab("bt");
	CheckCond(c89.Length == 6);

	var c90 = hm.AddRectContour(new P2(1, 1.2), new P2(2, 2), new int[]{0, 1, 2, 3});
	var c91 = hm.AddRectContour(new P2(0, 0), new P2(1.5, 1.5), new int[]{0, 1, 2, 3});
	var c92 = hm.UniteContours(new Hybmesh.Contour2D[]{c90, c91});
	var c93 = hm.DecomposeContour(c92);
	CheckCond(c93.Length == 3);
	var c94 = hm.PickContour(new P2(-1, -1), c93);
	c94.GetPoint(new P2(-1, -1));
	c94.GetPoint(null, new P2(0.5, -1));
	c94.GetPoint(null, null, new P2(1.0, -1));

	hm.ExportAllHmd("out.hmd", "ascii");
	hm.RemoveAll();
	hm.ImportGridHmg("out.hmd", "Grid2D_1");
	hm.Import3DGridHmg("out.hmd", "Grid3D_1");
	hm.RemoveAll();
}
public static void Main(){
	Hybmesh.hybmesh_exec_path = "../../src/py/hybmesh.py";
	Hybmesh.hybmesh_lib_path = "../../build/bin";
	using (var hm = new Hybmesh()){
		pings(hm);
	}
	Console.WriteLine("Done");
}

}
