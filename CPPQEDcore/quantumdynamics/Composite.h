// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#pragma once

#include "ExpectationValues.h"
#include "Liouvillian.h"

#include "SliceIterator.h"


using namespace quantumdynamics;

template <size_t RANK, size_t ... ra> requires ( sizeof...(ra) < RANK )
struct Broadcaster
{
  Broadcaster(Extents<RANK> extents) : offsets{cso<ra...>(extents)} {}

  const std::vector<size_t> offsets;

  Lindblad<RANK> operator()(const Lindblad<sizeof...(ra)>& l)
  {
    static constexpr std::array retainedAxes{ra...};
    static constexpr size_t RRANK = std::size(retainedAxes);

    return {

      .label{l.label} ,

      .jump{ [o=offsets,&l] (double t, StateVectorView<RANK> psi) {
        for (auto&& psiElem : sliceRange<retainedAxes>(psi,o)) applyJump(l.jump,t,psiElem);
      } } ,

      .rate{ [o=offsets,&l] (double t, StateVectorConstView<RANK> psi) {
        return partialTrace<retainedAxes,RANK>(LDO<StateVector,RANK>{psi},
                                               o,
                                               [&] (StateVectorConstView<RANK-std::size(retainedAxes)> psiElem) {return calculateRate(l.rate,t,psiElem); },
                                               std::plus{} );
      } } ,

      .superoperator{
        [ o=offsets, matrixOffsets=std::vector<size_t>{}, &l ] (double t, DensityOperatorConstView<RANK> rho, DensityOperatorView<RANK> drhodt) mutable {
          // matrixOffsets is populated when the lambda is first called, row offsets fastest as in storage (#90)
          if (!matrixOffsets.size()) {
            matrixOffsets.resize(sqr(o.size()));
            size_t extent=std::lround(std::sqrt(rho.dataView.size()));
            auto i=matrixOffsets.begin();
            for (auto v=o.begin(); v!=o.end(); ++v )
              for (auto u=o.begin(); u!=o.end(); (*i++) = (*u++) + extent * (*v) ) ;
          }

          for ( auto&& [rho,drhodt] : std::views::zip( sliceRange<extendedAxes<retainedAxes,RANK>>(rho,matrixOffsets),
                                                       sliceRange<extendedAxes<retainedAxes,RANK>>(drhodt,matrixOffsets) ) )
            applySuperoperator<RRANK>(l.superoperator,t,rho,drhodt) ;
        }
      }

    };
  }


  Liouvillian<RANK> operator()(const Liouvillian<sizeof...(ra)>& l)
  {
    Liouvillian<RANK> result( size(l) );
    std::ranges::transform( l, result.begin(), [this] (const Lindblad<sizeof...(ra)>& v) {return (*this)(v);} );
    return result;
  }


  // TODO: taking S here as const reference is inconsistent with the definitions of the hamiltonian and expectation_values concepts
  template<typename S>
  auto operator() (const S& s)
  {
    static constexpr std::array retainedAxes{ra...};
    static constexpr size_t RRANK = std::size(retainedAxes);

    if      constexpr (hamiltonian<S,RRANK>)
      return [o=offsets,&s] ( double t, StateVectorConstView<RANK> psi, StateVectorView<RANK> dpsidt, double t0 ) {
        for ( auto&& [psi,dpsidt] : std::views::zip( sliceRange<retainedAxes>(psi,o), sliceRange<retainedAxes>(dpsidt,o) ) )
          applyHamiltonian(s,t,psi,dpsidt,t0);
      };
    else if constexpr (expectation_values<S,RRANK>)
      return [o=offsets,&s] ( double t, lazy_density_operator<RANK> auto matrix) {
        return partialTrace<retainedAxes,RANK>( matrix, o,
                                                [&] (auto psiElem) {return calculateExpectationValues<RRANK>(s,t,psiElem); },
                                                plusTDP{} );
      };
    else static_assert(always_false<S>::value, "Unsupported type in Broadcaster");
  }

};

