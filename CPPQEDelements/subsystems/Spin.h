// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#pragma once

#include "Liouvillian.h"

#include <multidiagonal/MultiDiagonal.h>

#include "Pars.h"



namespace spin {


using namespace ::quantumdynamics;
using MultiDiagonal = ::quantumoperator::MultiDiagonal<1> ;


/// @brief Spin-s operators on the basis index i = m + s, i.e. the number of quanta above m = −s.
/// This is the Holstein–Primakoff map onto a harmonic-oscillator mode with i = n:
/// S_z = n − s, S⁺ = a†√(2s − a†a), S⁻ = √(2s − a†a) a, so S⁻ ≈ √(2s) a for large s, with `mode::aOp`.
/// S⁺ raises m, so it is the lower diagonal (offset −1). Spin ½ coincides with `qbit::multidiagonal`.
MultiDiagonal splus(size_t twoS, size_t dim=0);

MultiDiagonal sminus(size_t twoS, size_t dim=0);

MultiDiagonal sx(size_t twoS, size_t dim=0);
MultiDiagonal sy(size_t twoS, size_t dim=0);


MultiDiagonal sz(size_t twoS, size_t dim=0);


/*
inline auto sn(size_t dim, size_t twoS, dcomp z, bool isExact, double theta, double phi)
{
  return sin(theta)*(cos(phi)*sx(dim,twoS,z,isExact)+sin(phi)*sy(dim,twoS,z,isExact))+cos(theta)*sz(dim,twoS);
}

Tridiagonal sn(Ptr);


struct Pars
{
  size_t &twoS, &dim;
  double &theta, &phi;
  double &omega, &gamma;

  Pars(parameters::Table&, const std::string& ="");

};


class Liouvillian
  : public quantumdynamics::ElementLiouvillian<1,1,false>
{
protected:
  Liouvillian(size_t twoS, double gamma) : quantumdynamics::ElementLiouvillian<1,1,false>("Spin","S-"), twoS_(twoS), gamma_(gamma) {}
  
private:
  void   doActWithJ(quantumdynamics::NoTime,       quantumdynamics::freesystem::StateVectorLow     &) const ;
  double rate      (quantumdynamics::NoTime, const quantumdynamics::freesystem::LazyDensityOperator& matrix) const ; // {return -1.;}

  const size_t twoS_;
  const double gamma_;
  
};

*/
} // spin


/** \todo Implement some spherical-coordinates class and use it here */ /*
class SpinBase
  : public quantumdynamics::Free, 
    public quantumdynamics::ElementAveraged<1>
{
public:
  typedef quantumdynamics::freesystem::LazyDensityOperator LazyDensityOperator;

  SpinBase(size_t twoS, double theta, double phi, double omega, double gamma, size_t dim=0);

  size_t getTwoS() const {return twoS_;}

  double getTheta() const {return theta_;}
  double getPhi  () const {return   phi_;}

  double getOmega() const {return omega_;}
  double getGamma() const {return gamma_;}

  const dcomp get_z() const {return dcomp(gamma_,omega_);} ///< This plays analogous role to \f$z\f$ in Mode

private:
  void process_v(spin::Averages&) const;

  const spin::Averages average_v(quantumdynamics::NoTime, const LazyDensityOperator&) const;

  const size_t twoS_;

  const double theta_, phi_, omega_, gamma_, s_;

};

// NOTE if p.gamma is nonzero, Spin & SpinSch will contain the “Hamiltonian” part of the loss!

class Spin 
  : public SpinBase,
    public quantumdynamics::FreeExact<false>
{
public:

  Spin(const spin::Pars& p) 
    : SpinBase(p.twoS,0.,0.,p.omega,p.gamma,p.dim), FreeExact(getTotalDimension())
      // here, we need to assume that the Hamiltonian is diagonal in the given basis, that is why we set theta=phi=0.
  {}

private:
  void updateU(quantumdynamics::OneTime) const;

  bool applicableInMaster_v() const {return !getGamma();}

};


class SpinSch
  : public SpinBase,
    public quantumoperator::TridiagonalHamiltonian<1,false>
{
public:

  SpinSch(const spin::Pars& p) 
    : SpinBase(p.twoS,p.theta,p.phi,p.omega,p.gamma,p.dim),
      quantumoperator::TridiagonalHamiltonian<1,false>(-get_z()*spin::sn(getDimension(),p.twoS,get_z(),false,p.theta,p.phi))
  {
    getParsStream()<<"Schrodinger picture."<<std::endl;
  }

};


class DissipativeSpin
  : public Spin,
    public spin::Liouvillian
{
public:
  
  DissipativeSpin(const spin::Pars& p) : Spin(p), spin::Liouvillian(p.twoS,p.gamma) {}
  
};


class DissipativeSpinSch
  : public SpinSch,
    public spin::Liouvillian
{
public:
  
  DissipativeSpinSch(const spin::Pars& p) : SpinSch(p), spin::Liouvillian(p.twoS,p.gamma) {}
  
};
*/
