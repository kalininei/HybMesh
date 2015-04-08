#ifndef _AWTA_ADDALGO_H
#define _AWTA_ADDALGO_H
#include <type_traits>
#include <algorithm>
#include <numeric>
#include <vector>
#include <set>
#include <map>

namespace aa{

//================foreach
//foreach with multiple iterators increasing simultaneously
template<class Fun, typename Iter1, typename Iter2>
void foreach2(Iter1 i, Iter1 iend, Iter2 j, Fun&& method){
	while (i!=iend) method(*i++,*j++);
}

template<class Fun, typename Iter1, typename Iter2, typename Iter3>
void foreach3(Iter1 i, Iter1 iend, Iter2 j, Iter3 k, Fun&& method){
	while (i!=iend) method(*i++,*j++,*k++);
}
//foreach for all container elements
template<class Fun, class C1>
void Cforeach(C1&& c1, Fun&& method){
	auto i=c1.begin();
	while (i!=c1.end()) method(*i++);
}
template<class Fun, class C1, class C2>
void Cforeach2(C1&& c1, C2&& c2, Fun&& method){
	auto i=c1.begin();
	auto j=c2.begin();
	while (i!=c1.end()) method(*i++,*j++);
}
template<class Fun, class C1, class C2, class C3>
void Cforeach3(C1&& c1, C2&& c2, C3&& c3, Fun&& method){
	auto i=c1.begin();
	auto j=c2.begin();
	auto k=c3.begin();
	while (i!=c1.end()) method(*i++,*j++,*k++);
}
//variadic foreach
template<class Fun, typename Iter1, typename ... Iters>
void foreachN(Fun&& method, Iter1 i, Iter1 iend, Iters ... args){
	while (i!=iend) method(*i++, *args++ ...);
}
//foreach with index as a functor argument
template<class Fun, typename Iter1>
void foreach_ind(Iter1 i, Iter1 iend, Fun&& method){
	int ind=0;
	while (i!=iend) M(ind++,*i++);
}
template<class Fun, typename Iter1, typename Iter2>
void foreach2_ind(Iter1 i, Iter1 iend, Iter2 j, Fun&& method){
	int ind=0;
	while (i!=iend) method(ind++,*i++,*j++);
}
template<class Fun, typename Iter1, typename Iter2, typename Iter3>
void foreach3_ind(Iter1 i, Iter1 iend, Iter2 j, Iter3 k, Fun&& method){
	int ind=0;
	while (i!=iend) method(ind++,*i++,*j++,*k++);
}
//variadic foreach with index as a functor argument
template<class Fun, typename Iter1, typename ... IterN>
void foreachN_ind(Fun&& method, Iter1 i, Iter1 iend, IterN ... args){
	int ind=0;
	while (i!=iend) method(ind++,*i++, *args++ ... );
}

//==================== transform
//container and container
template<class Fun, class C1, class C2,
	class CHK1=typename std::remove_reference<C1>::type::allocator_type,
	class CHK2=typename std::remove_reference<C2>::type::allocator_type
>
void Ctransform(C1&& c1, C2&& c2, Fun&& method){
	std::transform(c1.begin(),c1.end(),c2.begin(),method); 
};
//container and input iterator
template<class Fun, class C1, typename Iter2,
	class CHK1=typename std::remove_reference<C1>::type::allocator_type,
	class CHK2=typename Iter2::iterator_category
>
void Ctransform(C1&& c1, Iter2 j, Fun&& method){
	std::transform(c1.begin(),c1.end(),j,method); 
};
//===================== accumulate
//using unary method
template<class Fun, class C1>
typename C1::value_type Caccumulate(C1&& c1, Fun&& method){
	return std::accumulate(c1.begin(),c1.end(),C1::value_type(0),
		[&method](double x, double y){
			return x+method(y);
		});
};

//====================== remove vector entries
template<class A>
void remove_entries(std::vector<A>& v, const std::set<int>& supi){
	std::vector<A> v1; v1.reserve(v.size()-supi.size());
	auto sit=supi.begin();
	for (size_t i=0; i<v.size(); ++i){
		if (sit==supi.end() || (int)i!=*sit) v1.push_back(v[i]);
		else (++sit);
	}
	v=v1;
}

//====================== Container<Data>  ->  Container<Data*>
template<template<class...> class Container, class Data,
	class PData=typename Container<Data>::iterator::pointer>  //(const Data*) for sets
Container<PData> pointers_container(Container<Data>& d){
	auto ret = Container<PData>();
	std::insert_iterator<Container<PData>> ins(ret, ret.begin());
	for (auto& v: d) ins = &v;
	return ret;
};

//====================== add a value to a Container<shared_ptr<Data>>
template<class Container, class Data=typename Container::value_type::element_type>
Data* add_shared(Container& c, const Data& d){
	auto p = new Data(d);
	std::inserter(c, c.end()) = typename Container::value_type(p);
	return p;
}


//====================== Function traits
template<class F>
struct function_traits: public function_traits<decltype(&F::operator())>{};

// function pointer
template<class R, class... Args>
struct function_traits<R(*)(Args...)> : public function_traits<R(Args...)> {};
template<class R, class... Args>
struct function_traits<R(&)(Args...)> : public function_traits<R(Args...)> {};
// functors
template<class R, class C, class... Args>
struct function_traits<R(C::*)(Args...)> : public function_traits<R(Args...)> {};
template<class R, class C, class... Args>
struct function_traits<R(C::*)(Args...) const> : public function_traits<R(Args...)> {};
 
//base
template<class R, class... Args>
struct function_traits<R(Args...)>{
	using result_type = R;

	static constexpr std::size_t arity = sizeof...(Args);

	template <std::size_t N>
	struct argument{
		using type = typename std::tuple_element<N,std::tuple<Args...>>::type;
	};
};

//======================= fill vector with Method(container values)
template<class Container, class Method, class ResTp = typename function_traits<Method>::result_type>
std::vector<ResTp> Cfill_vector(const Container& c, Method&& M){
	std::vector<ResTp> ret; ret.reserve(c.size());
	std::transform(c.begin(), c.end(), std::back_inserter(ret), 
			[&M](const typename Container::value_type& x){ return M(x); });
	return ret;
}

template<class Container, class OutContainer, class Method>
void Cfill_container(const Container& c, OutContainer& outc, Method&& M){
	std::transform(c.begin(), c.end(), std::inserter(outc, outc.end()), 
			[&M](const typename Container::value_type& x){ return M(x); });
}

// ====================== map with double as the key using epsilon compare
struct _MapComp{
	_MapComp(double _e): e(_e){}
	const double e;
	bool operator()(double a, double b){ return a+e<b; }
};

template <typename Arg>
class DoubleMap: public std::map<double, Arg, _MapComp>{
public:
	DoubleMap(double e): std::map<double, Arg, _MapComp>(_MapComp(e)){}
};

}//namespace

#endif
