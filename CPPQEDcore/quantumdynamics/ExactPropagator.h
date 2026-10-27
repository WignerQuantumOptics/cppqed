// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#pragma once

#include "StateVector.h"

#include <valarray>


namespace quantumdynamics {

using namespace quantumdata;

namespace exact_propagator_ns {

inline const struct NoOp {/* LogTree label{"noOp"};*/} noOp;

template <typename H, size_t RANK>
concept one_time_dependent_functional = requires (const H& h, double t, StateVectorView<RANK> psi) { h(t,psi); };

template <typename H, size_t RANK>
concept two_time_dependent_functional = requires (const H& h, double t, StateVectorView<RANK> psi, double t0) { h(t,psi,t0); };

template <typename H, size_t RANK>
concept functional = one_time_dependent_functional<H,RANK> || two_time_dependent_functional<H,RANK> || std::same_as<std::decay_t<H>,NoOp> ;

} // exact_propagator_ns

/// applying the propagators (interpreted as replacing |psi> with U|psi>)
template <size_t RANK, exact_propagator_ns::functional<RANK> T>
void applyPropagator(const T& h, double t, StateVectorView<RANK> psi, double t0)
{
  if      constexpr (exact_propagator_ns::one_time_dependent_functional<T,RANK>) h(t-t0,psi);
  else if constexpr (exact_propagator_ns::two_time_dependent_functional<T,RANK>) h(t,psi,t0);
}


template <typename H, size_t RANK>
concept exact_propagator = /* labelled<H> && */ exact_propagator_ns::functional<H,RANK>;


/// A unary propagator that assumes that the operator that transforms between the pictures is diagonal
template<bool IS_TWO_TIME=false>
struct UnaryDiagonalPropagator
{
  using Diagonal = std::valarray<dcomp>;

  using UpdateFunctional = std::conditional_t<IS_TWO_TIME,std::function<void(double,Diagonal&,double)>,std::function<void(double,Diagonal&)>>;

  UnaryDiagonalPropagator(std::string l, size_t dim, UpdateFunctional updateDiagonal) : label(l), diagonal(dim), updateDiagonal_{updateDiagonal} {}

  void operator()(double t, StateVectorView<1> psi, double t0) const requires (IS_TWO_TIME);// {if (t!=t_ || t0!=t0_) {updateDiagonal_(t_=t,diagonal,t0_=t0);} psi*=diagonal;}

  void operator()(double t, StateVectorView<1> psi) const requires (!IS_TWO_TIME);// {if (t!=t_) {updateDiagonal_(t_=t,diagonal);} psi*=diagonal;}

  const std::string label;

  mutable Diagonal diagonal;

private:
  mutable double t_=0, t0_=0;

  const UpdateFunctional updateDiagonal_;

};



namespace exact_propagator_ns {


template < auto retainedAxes, size_t RANK, functional<std::size(retainedAxes)> T >
void broadcast(const T& h, double t, StateVectorView<RANK> psi, double t0, const std::vector<size_t>& offsets)
{
  for (auto&& psiElem : sliceRange<retainedAxes>(psi,offsets)) applyPropagator(h,t,psiElem,t0);
}


} // exact_propagator_ns


} // structure
