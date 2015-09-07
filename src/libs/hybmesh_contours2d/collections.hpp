#ifndef HMCONT2D_TPP_COLLECTIONS_H
#define HMCONT2D_TPP_COLLECTIONS_H

#include "hmproject.h"
#include "bgeom2d.h"
#include "edges.hpp"

namespace HMCont2D{ namespace Tpp{

//Basis implementation of shared pointer generator
template<class C>
struct ShpGenerator{
	typedef shared_ptr<C> Tsh;
	//create a deep copy of value x
	virtual shared_ptr<C> deepcopy(const C* x) const {
		return std::make_shared<C>(*x);
	}

	//allocate shared pointer
	virtual shared_ptr<C> allocate() const{
		return shared_ptr<C>(new C());
	}
};

//Generator for cases when storage type is basic for original info type.
template<class RealType, class StorageType,
	class = typename std::enable_if<std::is_base_of<StorageType, RealType>::value>::type
>
struct CastGenerator: public ShpGenerator<StorageType>{
	typedef shared_ptr<StorageType> Ssto;

	virtual Ssto deepcopy(const StorageType* x) const override{
		const RealType* a = static_cast<RealType*>(x);
		RealType* b = new RealType(*a);
		return Ssto(b);
	}

	virtual Ssto allocate() const override{
		return Ssto(new RealType());
	}
};

//Shortcut for using SFINAE
template<class T1, class T2>
using IsBase = typename std::enable_if<std::is_base_of<T1, T2>::value>::type;
template<class T1, class T2>
using IsNotBase = typename std::enable_if<!std::is_base_of<T1, T2>::value>::type;

//Basic collection type: vector of shared pointers
template<class C>
struct Collection{
	//structs
	typedef C Tvalue;
	typedef shared_ptr<C> Tentry;
	struct TDeepCopyResult{
		std::map<Tvalue*, Tentry> oldnew; //old object ptr -> new object ptr
		std::map<Tvalue*, Tentry> newold; //new object ptr -> old object ptr
	};
	typedef ShpGenerator<Tvalue> ShpGen;

	//data
	vector<Tentry> data;

	//This object is used for creation of deep copies of values stored in data.
	//It should be replaced if type of values stored in data differs from C.
	//The latter requires pointer casting using object of CastGenerator<Derived from C, C>
	ShpGen _generator;

	//get methods
	int size() const{ return data.size(); }
	typename vector<Tentry>::iterator begin() {return data.begin();}
	typename vector<Tentry>::iterator end() {return data.end();}
	typename vector<Tentry>::const_iterator begin() const {return data.cbegin();}
	typename vector<Tentry>::const_iterator end() const {return data.cend();}

	Tvalue& value(int i){ return *data[i]; }
	Tvalue value(int i) const { return data[i]; }
	Tvalue* pvalue(int i) const {return data[i].get(); }
	vector<Tvalue*> pvalues() const {
		vector<Tvalue*> ret;
		std::transform(begin(), end(), std::back_inserter(ret), [](shared_ptr<Tvalue> a){ return a.get(); });
		return ret;
	}
	bool contains(const Tvalue* a){
		return std::any_of(data.begin(), data.end(),
			[&a](const Tentry& e){ return e.get() == a;
		});
	}

	//set methods: templates are used because
	//data internal structure can differ from TValue (derived from it) declared here.
	template<class Inp>
	using Valid = IsBase<Tvalue, Inp>;

	template<class CInp, class=Valid<CInp>>
	CInp* add_value(const CInp& v){
		data.push_back(_generator.deepcopy(&v));
		return static_cast<CInp*>(pvalue(size()-1));
	}

	template<class CInp, class=Valid<CInp>>
	CInp* add_value(shared_ptr<CInp> v){
		data.push_back(v);
		return static_cast<CInp*>(pvalue(size()-1));
	}

	template<class CInp=C, class=Valid<CInp>>
	shared_ptr<CInp> pop_value(){
		shared_ptr<CInp> ret(data.back());
		data.resize(size() - 1);
		return ret;
	}

	//Methods
	template<class CInp, class=Valid<typename CInp::Tvalue>>
	void Unite(const CInp& c){
		std::copy(c.begin(), c.end(), std::back_inserter(data));
	}

	//static Methods
	//DeepCopy using predefined generator of shared object.
	template<class TTarget, class = IsBase<Collection, TTarget>>
	static TDeepCopyResult DeepCopy(const Collection& from, TTarget& to, ShpGen& gen) {
		TDeepCopyResult res;
		std::transform(from.data.begin(), from.data.end(), std::back_inserter(to.data),
			[&gen](const Tentry& e){ return gen.deepcopy(e.get()); }
		);
		std::transform(from.data.rbegin(), from.data.rend(), to.data.rbegin(),
				std::inserter(res.oldnew, res.oldnew.begin()),
			[](const Tentry& e1, const Tentry& e2){ return std::make_pair(e1.get(), e2); }
		);
		std::transform(from.data.rbegin(), from.data.rend(), to.data.rbegin(),
				std::inserter(res.newold, res.newold.begin()),
			[](const Tentry& e1, const Tentry& e2){ return std::make_pair(e2.get(), e1); }
		);
		return res;
	};

	//DeepCopy using generator in to._generator
	template<class TTarget, class = IsBase<Collection, TTarget>>
	static TDeepCopyResult DeepCopy(const Collection& from, TTarget& to){
		return DeepCopy(from, to, to._generator);
	}
};



} //Tpp


// ============= Basic points collection
struct PCollection: public Tpp::Collection<Point>{
	//get
	Point* point(int i) const { return pvalue(i); }

	//Methods
	static void SaveVtk(const PCollection& dt, const char* fn);

	//Scaling
	static ScaleBase Scale01(PCollection&);
	static void Scale(PCollection&, const ScaleBase& sc);
	static void Unscale(PCollection&, const ScaleBase& sc);
};




// ============= Basic edges collection
struct ECollection: public Tpp::Collection<Edge>{
	//get
	double length() const { return std::accumulate(data.begin(), data.end(), 0.0,
			[](double s, const Tentry& e){ return s + e->length(); });}
	Edge* edge(int i) const { return pvalue(i); }

	bool contains_point(const Point* a){
		return std::any_of(data.begin(), data.end(),
			[&a](const Tentry& e){ return e->contains(a); });
	}

	vector<Point*> all_points() const;

	//Methods
	static void SaveVtk(const ECollection& dt, const char* fn);
	static Point* FindClosestNode(const ECollection& dt, const Point& p);
};



};

#endif
