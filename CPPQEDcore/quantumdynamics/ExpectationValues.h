// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#pragma once

#include "LazyDensityOperator.h"

#include "TemporalDataPoint.h"


namespace quantumdynamics {

using namespace quantumdata;

namespace expectation_values_ns {

inline const struct NoOp
{
  LogTree label;
  
  auto operator() (auto) const {return hana::make_tuple();}

} noOp;


template <typename L, size_t RANK>
concept time_dependent_functional = requires (const L& l, LDO<StateVector,RANK> rho, double t)
{
  { l(t,rho) } -> temporal_data_point ;
};


template <typename L, size_t RANK>
concept time_independent_functional = requires (const L& l, LDO<StateVector,RANK> rho)
{ 
  { l(rho) } -> temporal_data_point ;
};


// Alternatively:
/* requires (const L& l) { requires (
requires (StateVectorConstView<RANK> rho) { { l(rho) } -> temporal_data_point ; } &&
requires (DensityOperatorConstView<RANK> rho) { { l(rho) } -> temporal_data_point ; } ) ; } ;*/

template <typename L, size_t RANK>
concept functional = time_dependent_functional<L,RANK> || time_independent_functional<L,RANK> ;

template <typename H, size_t RANK> constexpr bool recursiveTrait = [] {
  if constexpr (hana_sequence<H>)
    return !!hana::all_of(
      decltype(hana::transform(std::declval<H>(), hana::typeid_)){},
      [] <class T> (T) { return recursiveTrait<typename T::type,RANK>; } );
  else return functional<H,RANK> ;
} ();

} // expectation_values_ns


template <typename L, size_t RANK>
concept expectation_values = expectation_values_ns::recursiveTrait<L,RANK> ;


template <size_t RANK, expectation_values<RANK> EV>
auto calculateExpectationValues(const EV& ev, double t, lazy_density_operator<RANK> auto matrix)
{
  if      constexpr (hana_sequence<EV>) return hana::transform(ev, [&] (const auto& f) {return calculateExpectationValues<RANK>(f,t,matrix);}) ;
  else if constexpr (expectation_values_ns::time_dependent_functional<EV,RANK>) return ev(t,matrix);
  else if constexpr (expectation_values_ns::time_independent_functional<EV,RANK>) return ev(matrix);
  else static_assert(always_false<EV>::value, "Unsupported type in calculateExpectationValues");
}


} // structure
