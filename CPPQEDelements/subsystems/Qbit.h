// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{Defines the qbit-bundle (tackling the dynamics of a single qbit)}
#pragma once

#include "ExpectationValues.h"
#include "Liouvillian.h"

#include <multidiagonal/MultiDiagonal.h>

#include "Pars.h"


namespace qbit {

using namespace ::quantumdynamics;

inline auto diagonalH(dcomp z) { return [=] (StateVectorConstView<1> psi, StateVectorView<1> dpsidt) {dpsidt(1)-=z*psi(1);}; }

inline auto offDiagonalH(dcomp eta) { return [=] (StateVectorConstView<1> psi, StateVectorView<1> dpsidt) {dpsidt(0)+=-conj(eta)*psi(1); dpsidt(1)+=eta*psi(0);}; }

inline auto fullH(dcomp z, dcomp eta) { return [diag=diagonalH(z),offDiag=offDiagonalH(eta)] (StateVectorConstView<1> psi, StateVectorView<1> dpsidt) {
  diag(psi,dpsidt); offDiag(psi,dpsidt);
};}

// ::cppqedutils::LogTree label(decltype(diagonalH(dcomp{}))) { return {"diagonalH"}; }
// ::cppqedutils::LogTree label(decltype(offDiagonalH(dcomp{}))) { return {"offDiagonalH"}; }
// ::cppqedutils::LogTree label(decltype(fullH(dcomp{},dcomp{}))) { return {"fullH"}; }


UnaryDiagonalPropagator<> propagator(dcomp z);

TimeIndependentJump<1> sigmaJump(double gamma_m);
TimeIndependentJump<1> sigmaPlusJump(double gamma_p);
TimeIndependentJump<1> sigma_zJump(double gamma_phi);

TimeIndependentSuperoperator<1> sigmaSuperoperator(double gamma_m);
TimeIndependentSuperoperator<1> sigmaPlusSuperoperator(double gamma_p);
TimeIndependentSuperoperator<1> sigma_zSuperoperator(double gamma_p);

TimeIndependentRate<1> sigmaRate(double gamma_m);
TimeIndependentRate<1> sigmaPlusRate(double gamma_p);
/*TimeIndependentRate<1> sigma_zRate(double gamma_p);*/

inline Lindblad<1> loss(double gamma_m) {return {"loss", sigmaJump(gamma_m), sigmaRate(gamma_m), sigmaSuperoperator(gamma_m)};}
inline Lindblad<1> gain(double gamma_p) {return {"gain", sigmaPlusJump(gamma_p), sigmaPlusRate(gamma_p), sigmaPlusSuperoperator(gamma_p)};}
// Lindblad<1> dephasing(double gamma_phi) {return {"dephasing", sigma_zJump(gamma_phi), sigma_zRate(gamma_phi), sigma_zSuperoperator(gamma_phi)};}


static constexpr auto expectationValues = [] (lazy_density_operator<1> auto rho)
{
  return hana::make_tuple( slv<"population0">(_(rho,0)), slv<"polarization">(_(rho,0,1)) );
};


inline auto make(dcomp zSch/*, dcomp zI*/, dcomp eta, double gamma_m, double gamma_p/*, double gamma_phi*/, const json::object& descr)
{
  Liouvillian<1> liouvillian;

  if (gamma_m) liouvillian.push_back(loss(gamma_m));
  if (gamma_p) liouvillian.push_back(gain(gamma_m));
  // if (gamma_phi) liouvillian.emplace(dephasing(gamma_phi));

  return std::make_tuple(2, fullH(zSch,eta), liouvillian, expectationValues, descr); //makeHamiltonianCollection<1>(diagonalH(zSch)/*,{"diagI",propagator(zI)}*/,offDiagonalH(eta)),
}


inline auto make(double delta, dcomp eta, double gamma_m, double gamma_p, const json::object& descr)
{
  return make(dcomp{gamma_m-gamma_p,-delta},eta,gamma_m,gamma_p,descr);
}



struct Pars : ::parameters::JSONizable
{
  double delta, gamma_m, gamma_p;
  dcomp eta, init;

  Pars(popl::OptionParser& op, std::string mod="")
  {
    using namespace ::parameters;
    add(mod,op,tjc,"Qubit",
      _("delta","detuning",-10.,delta),
      _("eta","drive",0.,eta),
      _("gamma_m","decay rate",10.,gamma_m),
      _("gamma_p","gain rate",0.,gamma_p),
      _("init","initial condition",dcomp(0.),init,noJSON));
  }

};


StateVector<1> state0();
StateVector<1> state1();
StateVector<1> init(dcomp psi1);

inline StateVector<1> init(const Pars& p) {return init(p.init);}


inline auto make(const Pars& p)
{
  return make(p.delta,p.eta,p.gamma_m,p.gamma_p,p.jsonize());
}


// operators are defined here as MultiDiagonals, which is an overkill, of course
/**
 * TODO: define qbit operators as Sigmas, and define composition of MultiDiagonal with Sigma
 */
/// @brief Spin-½ operators, in the convention of `spin` for 2s = 1: |0⟩ is down, |1⟩ is up (excited), so
/// σ⁺ = |1⟩⟨0| is the lower diagonal and σz = diag(−½, +½); `sz()` and `spin::sz(1)` agree.
namespace multidiagonal {

using namespace ::quantumdynamics;
using MultiDiagonal = ::quantumoperator::MultiDiagonal<1> ;


MultiDiagonal splus();

MultiDiagonal sminus();

MultiDiagonal sx();
MultiDiagonal sy();


MultiDiagonal sz();

}

} // qbit

