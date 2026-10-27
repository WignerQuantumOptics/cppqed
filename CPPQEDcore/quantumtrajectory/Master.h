// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#pragma once

#include "QuantumSystemDynamics.h"
#include "QuantumTrajectory.h"



namespace quantumtrajectory {


/// Auxiliary tools to Master
namespace master { typedef ode::Pars<> Pars; } // master


/// Master equation evolution from a \link DensityOperator density-operator\endlink initial condition
/**
 * \see \ref masterequation.
 * 
 * \note The ODE driver underlying this class needs to store several (typically 6–7, this depends on the chosen driver) density-operator instants.
 *
 */
template <size_t RANK,
          quantum_system_dynamics<RANK> QSD,
          ode::engine<StorageType> OE>
struct Master
{
  Master(auto&& qsd, auto&& rho, auto&& oe)
    : qsd{FWD(qsd)}, rho{FWD(rho)}, oe{FWD(oe)},
      rowIterationOffsets_{calculateSlicesOffsets<compileTimeOrdinals<RANK>>(this->rho.extents)} {}

  double time=0., time0=0.;
  QSD qsd;
  DensityOperator<RANK> rho;
  OE oe;

  friend double getDtDid(const Master& m) {return getDtDid(m.oe);}

  friend double getTime(const Master& m) {return m.time;}

  /// TODO: put here the system-specific things
  friend LogTree logIntro(const Master& m)
  {
    return {{"Master",{{"odeEngine",logIntro(m.oe)},{"System","TAKE FROM SYSTEM"}}}};
    // streamCharacteristics(m.qsd,streamParameters(m.qsd,streamIntro(m.oe,os)<<"Solving Master equation."<<endl<<endl ) )<<endl;
  }

  friend LogTree logOutro(const Master& m) {return logOutro(m.oe);}

  friend LogTree step(Master& m, double deltaT)
  {
    // TODO: think over what happens here when there is no Hamiltonian to apply – it’s probably OK, since at least one superoperator is always there
    // Moreover, I’m not sure anymore that these algorithms should be prepared for such special cases.
    auto res = step ( m.oe, deltaT, [&] (const StorageType& rhoRaw, StorageType& drhodtRaw, double t)
    {
      DensityOperatorConstView<RANK> rho{m.rho.extents,m.rho.strides,0,rhoRaw};
      DensityOperatorView<RANK> drhodt{m.rho.extents,m.rho.strides,0,drhodtRaw.begin(),std::ranges::fill(drhodtRaw,0)};

      hamiltonian_ns::broadcast<compileTimeOrdinals<RANK>>(getHa(m.qsd),t,rho,drhodt,m.time0,m.rowIterationOffsets_);

      twoTimesRealPartOfSelf(drhodt);

      for (auto&& lindblad : getLi(m.qsd))
        applySuperoperator<RANK>(lindblad.superoperator, t, rho, drhodt);

    },m.time,m.rho.dataStorage());

    // exact propagation
    exact_propagator_ns::broadcast<compileTimeOrdinals<RANK>>(getEx(m.qsd),m.time,m.rho.mutableView(),m.time0,m.rowIterationOffsets_);
    hermitianConjugateSelf(m.rho.mutableView());
    exact_propagator_ns::broadcast<compileTimeOrdinals<RANK>>(getEx(m.qsd),m.time,m.rho.mutableView(),m.time0,m.rowIterationOffsets_);
    conj(m.rho);
    m.time0=m.time;

    // The following "smoothening" of rho_ has proven necessary for the algorithm to remain stable:
    // We make the approximately Hermitian and normalized rho_ exactly so.

    twoTimesRealPartOfSelf(m.rho.mutableView());

    //std::cerr<<trace(m.rho)<<" ";

    renorm(m.rho);

    //std::cerr<<trace(m.rho)<<std::endl;

    return res;
  }

  friend auto temporalDataPoint(const Master& m)
  {
    return calculateExpectationValues<RANK>( getEV(m.qsd), m.time, LDO<DensityOperator,RANK>(m.rho) );
  }

  friend iarchive& readFromArrayOnlyArchive(Master& m, iarchive& iar) {return iar & m.rho;} // MultiArray can be (de)serialized

  /** structure of Master archives:
  * metaData – array – time – ( odeStepper – odeLogger – dtDid – dtTry )
  * state should precede time in order to be compatible with array-only archives
  */
  template <typename Archive>
  friend auto& stateIO(Master& m, Archive& ar)
  {
    stateIO(m.oe, ar & m.rho & m.time);
    // The following is a no-op in the case of state output, since time0 is equated with time @ the end of each step, however, it is an essential step for state input!
    m.time0=m.time;
    return ar;
  }


private:
  const std::vector<size_t> rowIterationOffsets_;

  // const DensityOperatorStreamer<RANK,V> dos_;

};


template<typename QSD, typename DO, typename OE>
Master(QSD&& , DO&& , OE&& ) -> Master< multiArrayRank_v<DO>/2, QSD, OE >;



namespace master {


template <template<typename> class OE, typename QSD, typename DO>
auto make(QSD&& qsd, DO&& state, const Pars& p)
{
  constexpr size_t RANK=multiArrayRank_v<std::decay_t<DO>>/2;
  using ODE=OE<StorageType>;

  double iDt=initialTimeStep(getFreqs(qsd)); // precalculate, since qsd gets forwarded (i.e., potentially moved)

  return Master<RANK,QSD,ODE>{FWD(qsd),FWD(state),ODE{iDt,p.epsRel,p.epsAbs}};
}


} // master


} // quantumtrajectory


template <size_t RANK, typename QSD, typename OE>
struct cppqedutils::trajectory::MakeSerializationMetadata<quantumtrajectory::Master<RANK,QSD,OE>>
{
  static auto _() {return SerializationMetadata{"CArray","Master",2*RANK};}
};


/*


template<int RANK, typename ODE_Engine, typename V>
std::ostream& quantumtrajectory::Master<RANK,ODE_Engine,V>::streamParameters(std::ostream& os) const
{
  using namespace std;

  streamCharacteristics(sys_, sys_->streamParameters(ode_.streamParameters(os)<<"Solving Master equation."<<endl<<endl ) )<<endl;

  if (const auto li=castLi(sys_)) {
    os<<"Decay channels:\n";
    {
      size_t i=0;
      li->streamKey(os,i);
    }
    os<<"Explicit superoperator calculations: ";
    DensityOperator rhotemp(rho_.getDimensions());
    {
      int n=0;
      for (int i=0; i<li->nAvr(); ++i)
        try {
          li->actWithSuperoperator(0,rho_.getArray(),rhotemp.getArray(),i);
          os<<i<<' '; ++n;
        }
        catch (const SuperoperatorNotImplementedException&) {}
      if (!n) os<<"none";
    }
    os<<endl;

  }
  
  return os;
}

*/
