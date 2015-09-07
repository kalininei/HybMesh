#ifndef HYBMESH_CONTAINERS_H
#define HYBMESH_CONTAINERS_H
#include "hmproject.h"
#include "collections.hpp"

namespace HMCont2D{

//Presents EdgeCollection which owns all its points
template<class ECol>
struct Container: public ECol{
	typedef ECol TParent;
	typedef PCollection PCol;
	typedef ECollection::ShpGen EShpGen;
	typedef PCollection::ShpGen PShpGen;
	struct TDeepCopyResult{
		std::map<typename PCol::Tvalue*, typename PCol::Tentry> points_oldnew;
		std::map<typename PCol::Tvalue*, typename PCol::Tentry> points_newold;
		std::map<typename ECol::Tvalue*, typename ECol::Tentry> edges_oldnew;
		std::map<typename ECol::Tvalue*, typename ECol::Tentry> edges_newold;
	};
	struct TDeepEGenerator;

	PCol pdata;

	//get/set
	Point* point(int i){return pdata.pvalue(i);}

	//Methods
	template<class TFrom>
	void Unite(const TFrom& c){
		ECol::Unite(c);
		pdata.Unite(c.pdata);
	}
	static double Area(const Container& obj){ return ECol::Area(obj); }

	//Scaling
	static ScaleBase Scale01(Container& c){ return PCol::Scale01(c.pdata); }
	static void Scale(Container& c, const ScaleBase& sc) { PCol::Scale(c, sc); }
	static void Unscale(Container& c, const ScaleBase& sc) { PCol::Unscale(c, sc); }


	//DeepCopy: From another container/from another ECollection.
	template<class TFrom, class = Tpp::IsBase<TParent, TFrom>>
	static typename TFrom::TDeepCopyResult
	DeepCopy(const TFrom& from, Container& to, EShpGen& egen, PShpGen& pgen);

	//version with default generators (to._generator, to.pdata._generator)
	template<class TFrom, class = Tpp::IsBase<TParent, TFrom>> 
	static typename TFrom::TDeepCopyResult 
	DeepCopy(const TFrom& from, Container& to);
};


// ================================ DeepCopyImplementation
namespace ContainersDeepCopy{

//Edge shared_ptr Generator which changes points pointers 
template<class Cont>
struct TDeepEGenerator: public Cont::EShpGen{
	typedef typename Cont::EShpGen EShpGen;
	typedef typename Cont::PCol PCol;
	typedef typename Cont::TParent ECol;
	typedef std::map<typename PCol::Tvalue*, typename PCol::Tentry> TMap;
	EShpGen* pogen;
	TMap* pmp;
	TDeepEGenerator(EShpGen& ogen, TMap& mp): EShpGen(), pogen(&ogen), pmp(&mp){}

	shared_ptr<typename ECol::Tvalue> deepcopy(const typename ECol::Tvalue* x) const override {
		auto r = pogen->deepcopy(x);
		r->pstart = (*pmp)[r->pstart].get();
		r->pend = (*pmp)[r->pend].get();
		return r;
	}
};

//TFrom - ECollection type
template<class TFrom, class TTo, class = void>
struct Core{
	typedef typename TTo::EShpGen EShpGen;
	typedef typename TTo::PShpGen PShpGen;
	typedef typename TTo::PCol PCol;
	typedef TTo ECol;
	static
	typename TFrom::TDeepCopyResult
	exe(const TFrom& from, TTo& to, EShpGen& egen, PShpGen& pgen){
		//copy points
		vector<Point*> allpnt = from.all_points();
		std::map<Point*, shared_ptr<Point>> oldnew;
		std::transform(allpnt.begin(), allpnt.end(), std::back_inserter(to.pdata.data),
			[&](Point* p){
				auto newp = pgen.deepcopy(p);
				oldnew[p] = newp;
				return newp;
			}
		);
		//set generator so it can deal with new points pointers
		TDeepEGenerator<TTo> deepgen(egen, oldnew);
		//make edges copy
		return TFrom::DeepCopy(from, to, deepgen);
	};

};

//TForm - Another container
template<class TFrom, class TTo>
struct Core<
	TFrom, TTo,
	Tpp::IsBase<typename TTo::TParent, typename TFrom::TParent>
>{ 
	typedef typename TTo::EShpGen EShpGen;
	typedef typename TTo::PShpGen PShpGen;
	typedef typename TTo::PCol PCol;
	typedef typename TTo::TParent ECol;
	static 
	typename TFrom::TDeepCopyResult
	exe(const TFrom& from, TTo& to, EShpGen& egen, PShpGen& pgen){
		//make points copy
		auto res1 = PCol::DeepCopy(from.pdata, to.pdata, pgen);
		//set generator so it can deal with new points pointers
		TDeepEGenerator<TTo> deepgen(egen, res1.oldnew);
		//make edges copy
		auto res2 = ECol::DeepCopy(from, to, deepgen);
		//assemble result
		return typename TFrom::TDeepCopyResult {res1.oldnew, res1.newold, res2.oldnew, res2.newold};
	}
};


}//Deepcopy Core


//DeepCopy: From another container/from another ECollection.
template<class ECol>
template<class TFrom, class>
typename TFrom::TDeepCopyResult
Container<ECol>::DeepCopy(const TFrom& from, Container& to, EShpGen& egen, PShpGen& pgen){
	return ContainersDeepCopy::Core<TFrom, Container>::exe(from, to, egen, pgen);
}

//version with default generators (to._generator, to.pdata._generator)
template<class ECol>
template<class TFrom, class> 
typename TFrom::TDeepCopyResult 
Container<ECol>::DeepCopy(const TFrom& from, Container& to){
	return ContainersDeepCopy::Core<TFrom, Container>::exe(from, to, to._generator, to.pdata._generator);
}


} //namespace

#endif
