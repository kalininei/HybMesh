#include <string>
#include <iostream>
#include <thread>
#include <Hybmesh.hpp>

//parameters of basic grid: structured grid in a square
double x1=0.0;  //bottom left point
double y1=0.0;

int nx=10;      //segmentation along axis
int ny=10;

double Lx=1.0;  //sizes of the square
double Ly=1.0;

//parameters of secondary grid: ring-type grid
double x2=0.5;       //ring center
double y2=0.5;

double innerrad=0.1; //ring radii
double outerrad=0.3;

double arcstep=0.02; //step size along the arc and radius of a ring
double radstep=0.02;

//unite grids options
double bufsize=0.1;      //Size of the buffer.

bool empty_holes=true;   //Leave ring internals blank.

bool fix_bnd=false;      //Allow shift of secondary grid boundary vertices:
                         //   only if secondary and basic grid boundaries intersect.
                         //Allow remeshing of basic grid boundary:
                         //   only if buffer zone includes basic grid boundary

double angle0=5;         //least significant angle: only if grid boundaries are remeshed.

std::string filler="3";  //Fill buffer with triangles

//path to hybmesh executable
std::string hmexec="../../../src/py/hybmesh.py";

//builds secondary grid and saves result into outfn
void grid_secondary(std::string outfn){
	//initialize builder
	Hybmesh builder2(hmexec);
	
	//calculate number of segmens in arc and radius directions.
	int na = (2.0 * 3.1415926 * innerrad)/arcstep;
	int nr = (outerrad-innerrad)/radstep;

	//build uniform ring grid prototype assigning
	//5 as its boundary feature
	Hybmesh::Grid2D g2=builder2.add_unf_ring_grid(
		Hybmesh::Point2(x2, y2),
		innerrad, outerrad,
		na, nr, 1.0, {5, 5});

	//export to native ascii file.
	//Initializer list {g2} converts g2 to single valued vector.
	builder2.export_grid_hmg({g2}, outfn);
}


int main(){
	//initialize main builder
	Hybmesh builder1(hmexec);

	//launch secondary grid builder in a separate thread
	std::thread builder2(grid_secondary, "g2.hmg");

	//build basic grid in the main thread assigning boundary types
	//1, 2, 3, 4 for left, bottom, right and top sides respectively.
	Hybmesh::Grid2D g1=builder1.add_unf_rect_grid(
		Hybmesh::Point2(x1, y1),
		Hybmesh::Point2(x1+Lx, y1+Ly),
		nx, ny,
		{1, 2, 3, 4});

	//wait until secondary grid is built
	//and read it from the file
	builder2.join();
	std::vector<Hybmesh::Grid2D> g2v=builder1.import_grid_hmg("g2.hmg");
	Hybmesh::Grid2D g2=g2v[0];

	//turn on console output for unite operation
	builder1.stdout_verbosity(3);

	try{
		//unite
		Hybmesh::Grid2D result=builder1.unite_grids1(
			g1, g2, bufsize,
			empty_holes, fix_bnd, angle0, filler);

		//save the result
		builder1.export_grid_vtk(result, "result.vtk");

		std::cout<<"Done"<<std::endl;
	} catch (Hybmesh::ERuntimeError& e){
		//controlled build exception has occured
		std::cout<<"Failed"<<std::endl;
		std::cout<<e.what()<<std::endl;
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
