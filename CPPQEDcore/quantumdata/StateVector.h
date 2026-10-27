// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFileDefault
#pragma once

#include <multidiagonal/MultiArrayComplex.h>


namespace quantumdata {


using namespace ::cppqedutils;

template <size_t RANK>
using Dimensions=Extents<RANK>;


using StorageType = multiarray::StorageType<dcomp>;


template<typename>
constexpr auto multiArrayRank_v=std::nullopt;

inline constexpr auto noInit = multiarray::noInit<dcomp>;

inline constexpr auto zeroInit = multiarray::zeroInit<dcomp>;


template <typename Derived>
struct VectorSpaceOperatorsGenerator
{

  friend auto operator+(const Derived& a, const Derived& b)
  {
    checkExtents(a,b,"VectorSpaceOperatorsGenerator plus");
    return Derived{getDimensions(a), [&] (size_t e) {
      auto r{noInit(e)}; for (auto&& [re,ae,be] : std::views::zip(r,a.dataView,b.dataView)) re=ae+be; return r;
    }};
  }


  friend auto operator-(const Derived& a, const Derived& b)
  {
    checkExtents(a,b,"VectorSpaceOperatorsGenerator minus");
    return Derived{getDimensions(a), [&] (size_t e) {
      auto r{noInit(e)}; for (auto&& [re,ae,be] : std::views::zip(r,a.dataView,b.dataView)) re=ae-be; return r;
    }};
  }


  friend auto operator*(dcomp v, const Derived& a)
  {
    return Derived{getDimensions(a), [&] (size_t e) {
      auto r{noInit(e)}; for (auto&& [re,ae] : std::views::zip(r,a.dataView)) re=v*ae; return r;
    }};
  }

};


template <size_t RANK> struct DensityOperator; // forward declaration


template <size_t RANK>
using StateVectorView=MultiArrayView<dcomp,RANK>;


template <size_t RANK>
using StateVectorConstView=MultiArrayConstView<dcomp,RANK>;



/// State vector of arbitrary arity
/**
 * Owns its data. Should be present on the highest level of the simulations in a single copy
 * (representing the initial condition, to be evolved by the quantum dynamics).
 */
template<size_t RANK>
struct StateVector : MultiArray<dcomp,RANK>, private VectorSpaceOperatorsGenerator<StateVector<RANK>>
{
  StateVector(const StateVector&) = delete; StateVector& operator=(const StateVector&) = delete;

  StateVector(StateVector&&) = default; StateVector& operator=(StateVector&&) = default;

  StateVector(MultiArray<dcomp,RANK> && ma) : MultiArray<dcomp,RANK>{std::move(ma)} {}

  StateVector(Dimensions<RANK> dimensions, auto initializer) : MultiArray<dcomp,RANK>{dimensions,initializer} {}

  explicit StateVector(Dimensions<RANK> dimensions)
    : StateVector{dimensions, [] (size_t e) {auto res{zeroInit(e)}; res[0]=1.; return res;}} {}

  /// Adds a dyad of the present object to `densityOperator`
  /**
   * This is done without actually forming the dyad in memory (so that this is not implemented in terms of StateVector::dyad),
   * which is important in situations when an average density operator is needed from an ensemble of state vectors, an example being quantumtrajectory::EnsembleMCWF. 
   */
  friend void addTo(const StateVector<RANK>& psi, DensityOperator<RANK>& rho, double weight=1.)
  {
    for (size_t size=psi.dataView.size(), i=0; i<size; ++i) for (size_t j=0; j<size; ++j)
      rho.dataStorage()[i+size*j]+=psi.dataView[i]*conj(psi.dataView[j]);
  }

};


template <size_t RANK>
double normSqr(const StateVector<RANK>& sv) {return frobeniusNormSqr(sv);}


template <size_t RANK>
double norm(const StateVector<RANK>& sv) {return frobeniusNorm(sv);}


template <size_t RANK>
double renorm(StateVector<RANK>& sv)
{
  double res=norm(sv);
  for (dcomp& v : sv.dataStorage()) v/=res;
  return res;
}


template <size_t RANK>
auto getDimensions(const StateVector<RANK>& psi) {return psi.extents;}


/// Forms a dyad with the argument, a rather expensive operation, implemented as a directProduct
template <size_t RANK>
auto dyad(const StateVector<RANK>& sv1, const StateVector<RANK>& sv2)
{
#ifndef NDEBUG
  if (sv1.extents != sv2.extents) throw std::runtime_error("Mismatch in StateVector::dyad dimensions");
#endif // NDEBUG
  return directProduct(sv1, sv2, [] (dcomp v1, dcomp v2) {return v1*conj(v2);} );
}


/// Creates the direct product, relying on the direct-product constructor
template<size_t RANK1, size_t RANK2>
inline auto operator*(const StateVector<RANK1>& psi1, const StateVector<RANK2>& psi2)
{
  return StateVector<RANK1+RANK2>(directProduct<RANK1,RANK2>(psi1,psi2));
}


/// Calculates the inner product, relying on StateVector::vectorView
template<size_t RANK>
dcomp braket(const StateVector<RANK>& psi1, const StateVector<RANK>& psi2);/*
{
  using blitz::tensor::i;
  linalg::CVector temp(psi1.getTotalDimension());
  temp=conj(psi1.vectorView()(i))*psi2.vectorView()(i);
  return sum(temp);
}*/


template <size_t RANK>
constexpr auto multiArrayRank_v<StateVector<RANK>> = RANK;


} // quantumdata


template <size_t RANK>
constexpr auto cppqedutils::passByValue_v<quantumdata::StateVector<RANK>> = false;


