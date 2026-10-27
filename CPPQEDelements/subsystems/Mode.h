// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#pragma once

#include "ExpectationValues.h"
#include "Liouvillian.h"

#include <multidiagonal/MultiDiagonal.h>

#include "Pars.h"



namespace mode {


using namespace ::quantumdynamics;
using MultiDiagonal = ::quantumoperator::MultiDiagonal<1> ;


MultiDiagonal aOp(size_t cutoff);

MultiDiagonal aDagOp(size_t cutoff);

MultiDiagonal nOp(size_t cutoff);


inline double photonnumber(StateVectorConstView<1> psi)
{
  double res=0;
  for (size_t n=0; n<psi.extents[0]; ++n) res+=n*sqrAbs(psi(n));
  return res;
}


inline auto hamiltonian(size_t cutoff, dcomp z, double omegaKerr, dcomp eta)
{
  MultiDiagonal ham, a{aOp(cutoff)}, aDag{aDagOp(cutoff)};

  if (abs(z)) ham-=z*nOp(cutoff);
  if (omegaKerr) ham+=(omegaKerr/1.i) * (aDag|aDag|a|a);
  if (abs(eta)) ham+=twoITimesImagPartOf(eta*aDag); // drive H = i(ηa† − η*a), so H/i = ηa† − η*a

  return ham;
}

inline auto hamiltonian(size_t cutoff, double delta, double omegaKerr, dcomp eta, double kappa, double nTh)
{
  return hamiltonian(cutoff,
                     {kappa*(2*nTh+1),-delta},
                     omegaKerr,eta);
}


TimeIndependentJump<1> aJump(double fact);
TimeIndependentJump<1> aDagJump(double fact);

TimeIndependentSuperoperator<1> aSuperoperator(double fact);
TimeIndependentSuperoperator<1> aDagSuperoperator(double fact);


Lindblad<1> photonLoss(double kappa, double nTh);
Lindblad<1> photonGain(double kappa, double nTh);


static constexpr auto expectationValues = [] (lazy_density_operator<1> auto rho)
{
  double pn=0, pnn=0; dcomp a=0;
  for (size_t n=1; n<rho.extents[0]; ++n) {
    pn +=n*  _(rho,n);
    pnn+=n*n*_(rho,n);
    a+=sqrt(n)*_(rho,n,n-1);
  }
  return hana::make_tuple(slv<"photon number">(pn),slv<"photon number square">(pnn),slv<"ladder operator">(a));
};


inline auto make(size_t cutoff, double delta, double omegaKerr, dcomp eta, double kappa, double nTh, json::object descr)
{
  Liouvillian<1> liouvillian;

  dcomp z{kappa*(2*nTh+1),-delta};

  if (kappa) {
    liouvillian.push_back(photonLoss(kappa,nTh));
    if (nTh) liouvillian.push_back(photonGain(kappa,nTh));
  }

  return std::make_tuple( cutoff, hamiltonian(cutoff,z,omegaKerr,eta), liouvillian, expectationValues, descr );

}



struct Pars : ::parameters::JSONizable
{
  size_t cutoff, initFock;
  double delta, omegaKerr, kappa, nTh;
  dcomp eta, init;

  Pars(popl::OptionParser& op, std::string mod="")
  {
    using namespace ::parameters;
    add(mod,op,tjc,"Mode",
        _("cutoff","Fock space cutoff",10,cutoff),
        _("delta","detuning",-10.,delta),
        _("omegaKerr","Kerr constant",0.,omegaKerr),
        _("eta","drive amplitude",dcomp(0),eta),
        _("kappa","decay rate",10.,kappa),
        _("nTh","thermal photon number",0.,nTh),
        _("initFock","Fock state initial condition",0,initFock,noJSON),
        _("init","Coherent state initial condition",dcomp(0),init,noJSON)
        );
  }

    // minitFock(p.add<size_t>("minitFock",mod,"Mode initial Fock state",0)),
    // minit(p.add<dcomp>("minit",mod,"Mode initial field",0)),

};


inline auto make(const Pars& p)
{
  return make(p.cutoff,p.delta,p.omegaKerr,p.eta,p.kappa,p.nTh,p.jsonize());
}


/// Coherent state
/**
 * The implementation relies on mathutils::coherentElement, which works also for high Fock-state elements
 *
 * \note The user has to take care that `alpha` is not too large for the given `cutoff` (rule of thumb: `cutoff>|alpha|^2`)
 */
StateVector<1> coherent(dcomp alpha, ///< amplitude
                        size_t cutoff ///< cutoff
                       );

StateVector<1> fock(size_t n, size_t dim, double phase=0);

/// Dispatcher for initial condition
StateVector<1> init(const Pars&);



/*


UnaryDiagonalPropagator<> propagator(dcomp z);


TimeIndependentSuperoperator<1> sigmaSuperoperator(double gamma_m);
TimeIndependentSuperoperator<1> sigmaPlusSuperoperator(double gamma_p);
TimeIndependentSuperoperator<1> sigma_zSuperoperator(double gamma_p);

TimeIndependentRate<1> sigmaRate(double gamma_m);
TimeIndependentRate<1> sigmaPlusRate(double gamma_p);

Lindblad<1> loss(double gamma_m) {return {"loss", sigmaJump(gamma_m), sigmaRate(gamma_m), sigmaSuperoperator(gamma_m)};}
Lindblad<1> gain(double gamma_p) {return {"gain", sigmaPlusJump(gamma_p), sigmaPlusRate(gamma_p), sigmaPlusSuperoperator(gamma_p)};}
// Lindblad<1> dephasing(double gamma_phi) {return {"dephasing", sigma_zJump(gamma_phi), sigma_zRate(gamma_phi), sigma_zSuperoperator(gamma_phi)};}


static constexpr auto expectationValues = [] (lazy_density_operator<1> auto rho) { return hana::make_tuple(real(_(rho,0,0)),_(rho,0,1)); };

::cppqedutils::LogTree label(decltype(expectationValues)) { return {"population","polarization"}; }


StateVector<1> state0();// {return mode::fock(0,2);}
StateVector<1> state1();// {return mode::fock(1,2);}
StateVector<1> init(dcomp psi1);
*/

} // mode
