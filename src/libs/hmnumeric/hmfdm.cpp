#include "hmfdm.hpp"
using namespace HMFdm;

void LaplasSolver::set_predef_value(int i, int j, double val){
	int gi = glob_index(i, j);
	auto fnd = predefined_values.find(gi);
	if (fnd == predefined_values.end()){
		//cannot change boundary stencil after solution
		//was initialized
		assert(was_init() == false);
		predefined_values.emplace(gi, val);
	} else {
		fnd->second = val;
	}
}

void LaplasSolver::SetBndValues(Bnd b, const std::function<double(int, int)>& f){
	switch (b){
		case Bnd::All:
			SetBndValues(Bnd::Bottom, f);
			SetBndValues(Bnd::Right, f);
			SetBndValues(Bnd::Top, f);
			SetBndValues(Bnd::Left, f);
			break;
		case Bnd::Bottom:
			for (int i=0; i<x.size(); ++i)
				set_predef_value(i, 0, f(i, 0));
			return;
		case Bnd::Right:
			for (int i=0; i<y.size(); ++i)
				set_predef_value(Nx()-1, i, f(Nx()-1, i));
			return;
		case Bnd::Top:
			for (int i=0; i<x.size(); ++i)
				set_predef_value(i, Ny()-1, f(i, Ny()-1));
			return;
		case Bnd::Left:
			for (int i=0; i<y.size(); ++i)
				set_predef_value(0, i, f(0, i));
			return;
	};
}

void LaplasSolver::Solve(vector<double>& ans){
	if (was_init() == false) initialize();
	assemble_rhs();
	solver->Solve(rhs, ans);
}

void LaplasSolver::assemble_rhs(){
	rhs.resize(N());
	std::fill(rhs.begin(), rhs.end(), 0.0);
	for (auto& v: predefined_values) rhs[v.first] = v.second;
}

void LaplasSolver::initialize(){
	HMMath::Mat m;
	m.data.resize(N());
	//internal, left, right
	for (int j=1; j<Ny()-1; ++j){
		double hy1 = y[j] - y[j-1];
		double hy2 = y[j+1] - y[j];
		double Ay = 1.0/hy1;
		double By = 1.0/hy2;
		double Cy = 2.0/(hy1 + hy2);
		//internal
		for (int i=1; i<Nx()-1; ++i){
			double gi = glob_index(i, j);

			double hx1 = x[i] - x[i-1];
			double hx2 = x[i+1] - x[i];
			double Ax = 1.0/hx1;
			double Bx = 1.0/hx2;
			double Cx = 2.0/(hx1 + hx2);

			m.set(gi, gi, Cx*(Ax+Bx) + Cy*(Ay+By));
			m.set(gi, gi-1, -Cx*Ax);
			m.set(gi, gi+1, -Cx*Bx);
			m.set(gi, gi-Nx(), -Cy*Ay);
			m.set(gi, gi+Nx(), -Cy*By);
		}
		//left
		{
			int i = 0;
			double gi = glob_index(i, j);
			double hx2 = x[i+1] - x[i];
			double Bx = 1.0/hx2;
			double Cx = 2.0/(hx2);

			m.set(gi, gi, Cx*Bx + Cy*(Ay+By));
			m.set(gi, gi+1, -Cx*Bx);
			m.set(gi, gi-Nx(), -Cy*Ay);
			m.set(gi, gi+Nx(), -Cy*By);
		}
		//right
		{
			int i = Nx()-1;
			double gi = glob_index(i, j);

			double hx1 = x[i] - x[i-1];
			double Ax = 1.0/hx1;
			double Cx = 2.0/(hx1);

			m.set(gi, gi, Cx*Ax + Cy*(Ay+By));
			m.set(gi, gi-1, -Cx*Ax);
			m.set(gi, gi-Nx(), -Cy*Ay);
			m.set(gi, gi+Nx(), -Cy*By);
		}
	}
	//top/bottom
	for (int i=1; i<Nx()-1; ++i){
		double hx1 = x[i] - x[i-1];
		double hx2 = x[i+1] - x[i];
		double Ax = 1.0/hx1;
		double Bx = 1.0/hx2;
		double Cx = 2.0/(hx1 + hx2);
		//top
		{
			int j = Ny() - 1;
			double gi = glob_index(i, j);

			double hy1 = y[j] - y[j-1];
			double Ay = 1.0/hy1;
			double Cy = 2.0/(hy1);

			m.set(gi, gi, Cx*(Ax+Bx) + Cy*Ay);
			m.set(gi, gi-1, -Cx*Ax);
			m.set(gi, gi+1, -Cx*Bx);
			m.set(gi, gi-Nx(), -Cy*Ay);
		}
	
		//bottom
		{
			int j = 0;
			double gi = glob_index(i, j);

			double hy2 = y[j+1] - y[j];
			double By = 1.0/hy2;
			double Cy = 2.0/(hy2);

			m.set(gi, gi, Cx*(Ax+Bx) + Cy*By);
			m.set(gi, gi-1, -Cx*Ax);
			m.set(gi, gi+1, -Cx*Bx);
			m.set(gi, gi+Nx(), -Cy*By);
		}
	}
	//bottom-left
	{
		int i=0, j=0, gi=0;
		double hx2 = x[i+1] - x[i];
		double hy2 = y[j+1] - y[j];
		double Bx = 1.0/hx2;
		double By = 1.0/hy2;
		double Cx = 2.0/hx2;
		double Cy = 2.0/hy2;

		m.set(gi, gi, Cx*Bx + Cy*By);
		m.set(gi, gi+1, -Cx*Bx);
		m.set(gi, gi+Nx(), -Cy*By);
	}
	//bottom right
	{
		int i=Nx()-1, j=0, gi=Nx()-1;
		double hx1 = x[i] - x[i-1];
		double hy2 = y[j+1] - y[j];
		double Ax = 1.0/hx1;
		double By = 1.0/hy2;
		double Cx = 2.0/hx1;
		double Cy = 2.0/hy2;

		m.set(gi, gi, Cx*Ax + Cy*By);
		m.set(gi, gi-1, -Cx*Ax);
		m.set(gi, gi+Nx(), -Cy*By);
	}
	//top left
	{
		int i=0, j=Ny()-1;
		int gi=glob_index(i,j);
		double hx2 = x[i+1] - x[i];
		double hy1 = y[j] - y[j-1];
		double Bx = 1.0/hx2;
		double Ay = 1.0/hy1;
		double Cx = 2.0/hx2;
		double Cy = 2.0/hy1;

		m.set(gi, gi, Cx*Bx + Cy*Ay);
		m.set(gi, gi+1, -Cx*Bx);
		m.set(gi, gi-Nx(), -Cy*Ay);
	}
	//top right
	{
		int i=Nx()-1, j=Ny()-1, gi = N()-1;
		double hx1 = x[i] - x[i-1];
		double hy1 = y[j] - y[j-1];
		double Ax = 1.0/hx1;
		double Ay = 1.0/hy1;
		double Cx = 2.0/hx1;
		double Cy = 2.0/hy1;

		m.set(gi, gi, Cx*Ax + Cy*Ay);
		m.set(gi, gi-1, -Cx*Ax);
		m.set(gi, gi-Nx(), -Cy*Ay);
	}

	//prediefined values
	for (auto& v: predefined_values){
		auto& row = m.data[v.first];
		row.clear();
		row[v.first] = 1;
	}

	solver = HMMath::MatSolve::Factory(m);
}
