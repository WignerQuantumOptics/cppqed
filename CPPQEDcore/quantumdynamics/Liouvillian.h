// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#pragma once

#include "ExactPropagator.h"
#include "FrameworkTraits.h"
#include "Hamiltonian.h"

#include "LazyDensityOperator.h"

#include <variant>

namespace quantumdynamics {


template <size_t RANK>
using TimeDependentJump = std::function<void(double t, StateVectorView<RANK> psi)>;

template <size_t RANK>
using TimeIndependentJump = std::function<void(StateVectorView<RANK> psi)>;

template <size_t RANK>
using Jump = std::variant<TimeDependentJump<RANK>,TimeIndependentJump<RANK>>;


template <size_t RANK>
void applyJump(Jump<RANK> jump, double t, StateVectorView<RANK> psi)
{
  std::visit(overload{
    [&] (const TimeDependentJump<RANK>& j) {j(t,psi);},
    [&] (const TimeIndependentJump<RANK>& j) {j(psi);}
  },jump);
}


template <size_t RANK>
using TimeDependentRate = std::function<double(double t, StateVectorConstView<RANK> psi)>;

template <size_t RANK>
using TimeIndependentRate = std::function<double(StateVectorConstView<RANK> psi)>;

template <size_t RANK>
using Rate = std::variant<TimeDependentRate<RANK>,TimeIndependentRate<RANK>>;


template <size_t RANK>
double calculateRate(Rate<RANK> rate, double t, StateVectorConstView<RANK> psi)
{
  return std::visit(overload{
    [&] (const TimeDependentRate<RANK>& r) {return r(t,psi);},
    [&] (const TimeIndependentRate<RANK>& r) {return r(psi);}
  },rate);
}


template <size_t RANK>
using TimeDependentSuperoperator = ODE_derivativeTimeDependentFunctional<DensityOperatorConstView<RANK>,DensityOperatorView<RANK>>;

template <size_t RANK>
using TimeIndependentSuperoperator = ODE_derivativeTimeIndependentFunctional<DensityOperatorConstView<RANK>,DensityOperatorView<RANK>>;

template <size_t RANK>
using Superoperator = std::variant<TimeDependentSuperoperator<RANK>,TimeIndependentSuperoperator<RANK>>;


template <size_t RANK>
void applySuperoperator(Superoperator<RANK> sup, double t, DensityOperatorConstView<RANK> rho, DensityOperatorView<RANK>& drhodt)
{
  return std::visit(overload{
    [&] (const TimeDependentSuperoperator<RANK>& s) {s(t,rho,drhodt);},
    [&] (const TimeIndependentSuperoperator<RANK>& s) {s(rho,drhodt);}
  },sup);
}


template <size_t RANK>
struct Lindblad
{
  static constexpr size_t N_RANK=RANK;

  std::string label;
  Jump<RANK> jump;
  Rate<RANK> rate;
  Superoperator<RANK> superoperator;
};


/*
 * liouvillian can be a concept like this:
 * template <typename L, size_t RANK> concept liouvillian = std::ranges::forward_range<L> && std::is_same_v<std::ranges::range_value_t<L>,Lindblad<RANK>>;
 * but `std::join_view`s are very limited, so probably we anyway need to join broadcasted subsystem Lindblads into a single monolithic container in composites
*/
template <size_t RANK> using Liouvillian = std::vector<Lindblad<RANK>>;


template <size_t RANK>
json::array toJSON_array(const Liouvillian<RANK>& li)
{
  json::array res(size(li));
  for (auto&& [r,l] : std::views::zip(res,li)) r=l.label;
  return res;
}


template <size_t RANK>
auto rateFromJump(double t, StateVectorConstView<RANK> psi, Jump<RANK> jump)
{
  StateVector<RANK> psiTemp{ copy(psi) };
  applyJump(jump,t,psiTemp.mutableView());
  return sqr(norm(psiTemp)); // we cannot return psiTemp because StateVector is not copyable
}


template <size_t RANK>
void superoperatorFromJump(double t, DensityOperatorConstView<RANK> rho, DensityOperatorView<RANK> drhodt, Jump<RANK> jump, const std::vector<size_t>& rowIterationOffsets)
{
  DensityOperator<RANK> rhoTemp{ copy(rho) }; // deep copy

  auto sr=sliceRange<compileTimeOrdinals<RANK>>(rhoTemp.mutableView(),rowIterationOffsets);

  auto unaryIteration=[&] { for (auto& psiTemp : sr) applyJump(jump,t,psiTemp); };

  unaryIteration(); hermitianConjugateSelf(rhoTemp.mutableView()); unaryIteration();

  for (auto& [t,o] : std::views::zip(drhodt.dataView,rhoTemp.dataView) ) t+=o;
}


auto concatenate(const auto&... lious)
{
  Liouvillian<
    [] <size_t RANK> ( const Liouvillian<RANK>&, const auto&... ) {return RANK;} (lious...)
    // TODO: this seems to be a very nice solution that should be generally applicable instead of metafunctions
  > result( ( size(lious) + ... ) );
  size_t index{};

  ((std::copy_n(lious.begin(), size(lious), result.begin() + index), index += size(lious)), ...);

  return result;
}



// /**
//  * The below concept-based structure mirrors ExpectationValues & Hamiltonian, but probably won’t have a lot of practical use.
//  */
//
// template <typename L, size_t RANK>
// concept time_dependent_jump = exact_propagator_ns::one_time_dependent_functional<L,RANK>;
//
// template <typename L, size_t RANK>
// concept time_independent_jump = requires(L&& l, StateVectorView<RANK> psi) { l(psi); };
//
// template <typename L, size_t RANK>
// concept jump = time_dependent_jump<L,RANK> || time_independent_jump<L,RANK>;
//
//
// /// Here, we could just use an alias to expectationvalues::functional, however, lindblad represents expressly just a single jump operator …
// /** … so that rate is just a single double. Accordingly, the label is only a single string as well. */
// template <typename L, size_t RANK>
// concept time_dependent_rate = requires (const L& l, double t, StateVectorConstView<RANK> psi) { { l(t,psi) } -> ::std::convertible_to<double>; } ;
//
// template <typename L, size_t RANK>
// concept time_independent_rate = requires (const L& l, StateVectorConstView<RANK> psi) { { l(psi) } -> ::std::convertible_to<double>; } ;
//
// template <typename L, size_t RANK> concept rate = time_dependent_rate<L,RANK> || time_independent_rate<L,RANK> ;
//
//
// template <typename L, size_t RANK>
// concept time_dependent_superoperator = ode_time_dependent_derivative<L,DensityOperatorConstView<RANK>,DensityOperatorView<RANK>>;
//
// template <typename L, size_t RANK>
// concept time_independent_superoperator = ode_time_independent_derivative<L,DensityOperatorConstView<RANK>,DensityOperatorView<RANK>>;
//
// template <typename L, size_t RANK>
// concept superoperator = time_dependent_superoperator<L,RANK> || time_independent_superoperator<L,RANK>;
//
//
// template <typename L, size_t RANK>
// concept lindblad_with_jump = labelled<L,std::string> && jump<L,RANK> ;
//
//
// template <typename L, size_t RANK>
// concept lindblad_with_rate = lindblad_with_jump<L,RANK> && rate<L,RANK> ;
//
//
// template <typename L, size_t RANK>
// concept lindblad_with_superoperator = lindblad_with_jump<L,RANK> && superoperator<L,RANK> ;
//
//
// template <typename L, size_t RANK>
// concept lindblad = lindblad_with_superoperator<L,RANK> || lindblad_with_rate<L,RANK> || lindblad_with_jump<L,RANK>;
//


} // structure
