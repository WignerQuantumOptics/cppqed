// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#pragma once

#include "FrameworkTraits.h"
#include "StateVector.h"


namespace quantumdynamics {

using namespace quantumdata;

/// TODO: this could be fused with ode::system, but let’s not go into that for the moment :)
template <typename H, typename StateIn, typename StateOut>
concept ode_time_independent_derivative = requires(const H& h, StateIn psi, StateOut dpsidt) { h(psi,dpsidt); };

template <typename H, typename StateIn, typename StateOut>
concept ode_time_dependent_derivative = requires(const H& h, double t, StateIn psi, StateOut dpsidt) { h(t,psi,dpsidt); };


namespace hamiltonian_ns {

inline const struct NoOp {/*LogTree label{"noOp"};*/} noOp;

template <typename H, size_t RANK>
concept time_independent_functional = ode_time_independent_derivative<H,StateVectorConstView<RANK>,StateVectorView<RANK>>;

template <typename H, size_t RANK>
concept one_time_dependent_functional = ode_time_dependent_derivative<H,StateVectorConstView<RANK>,StateVectorView<RANK>>;

template <typename H, size_t RANK>
concept two_time_dependent_functional = requires (const H& h, double t, StateVectorConstView<RANK> psi, StateVectorView<RANK> dpsidt, double t0) { h(t,psi,dpsidt,t0); };

template <typename H, size_t RANK>
concept functional = time_independent_functional<H,RANK> || one_time_dependent_functional<H,RANK> || two_time_dependent_functional<H,RANK> || std::same_as<std::decay_t<H>,NoOp>;

template <typename H, size_t RANK> constexpr bool recursiveTrait = [] {
  if constexpr (hana_sequence<H>)
    return !!hana::all_of(
      decltype(hana::transform(std::declval<H>(), hana::typeid_)){},
      [] <class T> (T) { return recursiveTrait<typename T::type,RANK>; } );
  else return functional<H,RANK> ;
} ();

} // hamiltonian_ns


template <typename H, size_t RANK>
concept hamiltonian = hamiltonian_ns::recursiveTrait<H,RANK>;

/// applying a Hamiltonian term is by default interpreted as |dpsidt>+=H|psi>/(i*hbar)
/** However, if indexing is done carefully, `psi` and `dpsidt` can refer to the same underlying data */
template <size_t RANK, hamiltonian<RANK> T>
void applyHamiltonian(const T& hamiltonian, double t, StateVectorConstView<RANK> psi, StateVectorView<RANK> dpsidt, double t0)
{
  if      constexpr (hana_sequence<T>) hana::for_each( hamiltonian, [&] (const auto& h) {applyHamiltonian(h,t,psi,dpsidt,t0);} );
  else if constexpr (hamiltonian_ns::  time_independent_functional<T,RANK>) hamiltonian(psi,dpsidt);
  else if constexpr (hamiltonian_ns::one_time_dependent_functional<T,RANK>) hamiltonian(t-t0,psi,dpsidt);
  else if constexpr (hamiltonian_ns::two_time_dependent_functional<T,RANK>) hamiltonian(t,psi,dpsidt,t0);
  else static_assert(always_false<T>::value, "Unsupported type in applyHamiltonian");
}


/// composition of two hamiltonian functionals, meant for interaction-type elements
/**
 * For simplicity, it always returns a hamiltonian_ns::two_time_dependent_functional
 * \note Here, we cannot really use the operator syntax as in MultiDiagonal as this function will not be found by ADL
 */
template <size_t RANK1, size_t RANK2>
auto compose(const hamiltonian_ns::functional<RANK1> auto& h1, const hamiltonian_ns::functional<RANK2> auto& h2);



/*
template <typename H, size_t RANK>
concept hamiltonian = hana_sequence<H> && !!hana::all_of(
  decltype(hana::transform(std::declval<H>(), hana::typeid_)){},
  []<class T>(T) { return hamiltonian_element<typename T::type,RANK>; });

*/



template <typename StateIn, typename StateOut>
using ODE_derivativeTimeIndependentFunctional = std::function<void(StateIn psi, StateOut dpsidt)> ;

template <typename StateIn, typename StateOut>
using ODE_derivativeTimeDependentFunctional = std::function<void(double t, StateIn psi, StateOut dpsidt)> ;


template <size_t RANK> using TimeIndependentTerm = ODE_derivativeTimeIndependentFunctional<StateVectorConstView<RANK>,StateVectorView<RANK>>;
template <size_t RANK> using TimeDependentTerm = ODE_derivativeTimeDependentFunctional<StateVectorConstView<RANK>,StateVectorView<RANK>>;


} // structure



#include "SparseMatrix.h"

#include "progressbar.hpp"

namespace quantumdynamics::hamiltonian_ns {


template <size_t RANK, hamiltonian<RANK> H>
quantumoperator::SparseMatrix vectorize(H&& h, Dimensions<RANK> d)
{
  quantumoperator::SparseMatrix::Elements elements;
  StateVector<RANK> psi{d,zeroInit}, dpsidt{d,zeroInit};
  auto psiView{psi.mutableView().dataView}, dpsidtView{dpsidt.mutableView().dataView};
  size_t dim=multiarray::calculateExtent(d);
  progressbar bar(dim);
  for ( size_t i=0; i<dim; (bar.update(), ++i) ) {
    psiView[i]=1.;
    applyHamiltonian(h,0.,psi,dpsidt.mutableView(),0.);
    for (size_t j=0; j<dim; ++j) if (dcomp v=dpsidtView[j]; abs(v)) elements.push_back({j,i,v});
    psiView[i]=0.; for (dcomp& v : dpsidtView) v=0.;
  }
  return {.elements{elements}};
}


} // quantumdynamics::hamiltonian_ns
