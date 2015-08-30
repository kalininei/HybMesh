#include "hybmesh_contours2d.h"
#include "hybmesh_contours2d.hpp"

int hybmesh_contours2d_ping(int a){
	return 2*a;
}

using namespace HMCont2D;


std::tuple<const Point*, const Point*, int> Contour::Edge(int i) const{
	const Point* p = Pnt(i);
	int b = btype.at(p);
	return std::tuple<const Point*, const Point*, int>(
		p, Pnt(i+1), b);
}

void Contour::AddPointToEnd(double x, double y, int b){
	shared_ptr<Point> newp(new Point(x, y));
	pts.push_back(newp);
	btype[newp.get()] = b;
}

void Contour::AddPointToEnd(Point xy, int b){
	AddPointToEnd(xy.x, xy.y, b);
}

bool Contour::HasCrosses(const Contour& c1, const Contour& c2){
	std::cout<<"DUMMY HasCrosses"<<std::endl;
	return false;
}
void Contour::PtsReallocate(){
	std::cout<<"DUMMY PtsReallocate"<<std::endl;
}

BoundingBox Contour::BuildBoundingBox() const{
	if (NumPoints() == 0) { return BoundingBox(0, 0, 0, 0); }
	auto ret = BoundingBox(pts[0]->x, pts[0]->y, pts[0]->x, pts[0]->y);
	for (auto p: pts) ret.WidenWithPoint(*p);
	return ret;
}

ClosedContour ClosedContour::DeepCopy() const{
	ClosedContour c2(*this);
	c2.PtsReallocate();
	return c2;
}

bool ClosedContour::IsValid() const{
	std::cout<<"DUMMY IsValid()"<<std::endl;
	return true;
}

vector<ClosedContour> ClosedContour::FromConnectedPoints(const vector<Point>& xy,
		const vector<int>& i0, const vector<int>& i1,
		const vector<int>& bnd){
	vector<ClosedContour> ret;
	typedef std::tuple<int, int, int> TEdge;
	auto g0 = [](const TEdge& e)->int{ return std::get<0>(e); };
	auto g1 = [](const TEdge& e)->int{ return std::get<1>(e); };
	auto gb = [](const TEdge& e)->int{ return std::get<2>(e); };
	//make a copy of edges from input
	std::set<TEdge> edges;
	for (int i=0; i<std::min(i0.size(), i1.size()); ++i){
		int b = i < bnd.size() ? bnd[i] : 0;
		edges.insert(std::make_tuple(i0[i], i1[i], b));
	}
	while (edges.size() > 0){
		//filter off closed contour: move'em from edges to econt
		vector<TEdge> econt;
		econt.push_back(*edges.begin()); edges.erase(edges.begin());
		int ilast = g1(econt[0]);
		while (1) {
			auto fnd = std::find_if(edges.begin(), edges.end(),
				[&](const TEdge& e)->bool{
					return (g0(e) == ilast || g1(e) == ilast);
				});
			if (fnd == edges.end()) break;
			ilast = (g0(*fnd) == ilast) ? g1(*fnd) : g0(*fnd);
			econt.push_back(*fnd); edges.erase(fnd); 
		}

		//check if contour is closed
		if (econt.size() < 2) continue;
		if (g0(econt[0]) != g0(econt.back()) && 
				g0(econt[0]) != g1(econt.back()) &&
				g1(econt[0]) != g0(econt.back()) &&
				g1(econt[0]) != g1(econt.back()))
			continue;

		//construct contour from sequence of edges
		vector<Point> vp;
		vector<int> vb;
		if (g1(econt[0]) == g0(econt.back()) ||
				g1(econt[0]) == g1(econt.back())){
			ilast = g1(econt[0]);
		} else {
			ilast = g0(econt[0]);
		}
		for (int i=0; i<econt.size(); ++i){
			vb.push_back(gb(econt[i]));
			if (g0(econt[i]) == ilast){
				vp.push_back(xy[g0(econt[i])]);
				ilast = g1(econt[i]);
			} else {
				vp.push_back(xy[g1(econt[i])]);
				ilast = g0(econt[i]);
			}
		}
		//assemble from sequence of nodes
		ret.push_back(ClosedContour());
		ClosedContour& c = ret.back();
		for (int i=0; i<vp.size(); ++i){
			c.AddPointToEnd(vp[i], vb[i]);
		}
	}
	return ret;
}

void ContourTree::AddContour(const ClosedContour& cont){
	//check for contour validity
	if (!cont.IsValid()) throw std::runtime_error("Invalid contour was added to ContourTree");
	//check if contour intersects other contours
	for (auto& c: conts){
		if (Contour::HasCrosses(cont, *c))
			throw std::runtime_error("Impossible to add contour to ContourTree due to intersections");
	}
	//add to contours list
	aa::add_shared(conts, cont.DeepCopy());
	//rebuild structure
	RebuildStructure();
}

void ContourTree::RebuildStructure(){
	std::cout<<"DUMMY RebuildStructure"<<std::endl;
}

int ContourTree::NumPoints() const{
	int ret = 0;
	for (auto& c: conts) ret += c->NumPoints();
	return ret;
}
int ContourTree::NumEdges() const{
	int ret = 0;
	for (auto& c: conts) ret += c->NumEdges();
	return ret;
}

BoundingBox ContourTree::BuildBoundingBox() const{
	vector<BoundingBox> bd;
	for (auto c: conts) bd.push_back(c->BuildBoundingBox());
	return BoundingBox(bd);
}


// ============================= C procedures
void* create_contour_tree(int Npnt, double* points, int Nedges, int* edges){
	try{
		vector<Point> p;
		for (int i=0; i<Npnt; ++i){
			p.push_back(Point(points[2*i], points[2*i+1]));
		}
		vector<int> i0, i1, b;
		for (int i=0; i<Nedges; ++i){
			i0.push_back(edges[3*i]);
			i1.push_back(edges[3*i+1]);
			b.push_back(edges[3*i+2]);
		}
		auto conts = ClosedContour::FromConnectedPoints(p, i0, i1, b);
		std::unique_ptr<ContourTree> ret(new ContourTree());
		for (auto& c: conts) ret->AddContour(c);
		return ret.release();
	} catch (std::exception& e) {
		std::cout<<e.what()<<std::endl;
		return 0;
	}
}

void free_contour_tree(void* tree){
	if (tree!=NULL) delete static_cast<ContourTree*>(tree);
}






