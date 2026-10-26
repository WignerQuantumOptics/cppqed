// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "Qbit.h"


using namespace cppqedutils; using namespace quantumdynamics;


UnaryDiagonalPropagator<> qbit::propagator(dcomp z)
{
  return {"diagonalPropagator", 2 , [=] (double t, auto& d) {d[0]=1; d[1]=exp(-z*t);} };
}


TimeIndependentJump<1> qbit::sigmaJump(double gamma_m) { return [=] (StateVectorView<1> psi) {psi(0)=sqrt(2.*gamma_m)*psi(1); psi(1)=0;}; }

TimeIndependentJump<1> qbit::sigmaPlusJump(double gamma_p) { return [=] (StateVectorView<1> psi) {psi(1)=sqrt(2.*gamma_p)*psi(0); psi(0)=0;}; }

TimeIndependentJump<1> qbit::sigma_zJump(double gamma_phi) { return [=] (StateVectorView<1> psi) {double fact=sqrt(2.*gamma_phi); psi(0)*=fact; psi(1)*=-fact;}; }


TimeIndependentSuperoperator<1> qbit::sigmaSuperoperator(double gamma_m)
{
  return [=] (DensityOperatorConstView<1> rho, DensityOperatorView<1> drhodt) {drhodt(0,0)+=2.*gamma_m*rho(1,1);};
}

TimeIndependentSuperoperator<1> qbit::sigmaPlusSuperoperator(double gamma_p)
{
  return [=] (DensityOperatorConstView<1> rho, DensityOperatorView<1> drhodt) {drhodt(1,1)+=2.*gamma_p*rho(0,0);};
}


TimeIndependentRate<1> qbit::sigmaRate(double gamma_m) {return [=] (StateVectorConstView<1> psi) {return 2.*gamma_m*sqrAbs(psi(1));}; }


TimeIndependentRate<1> qbit::sigmaPlusRate(double gamma_p) {return [=] (StateVectorConstView<1> psi) {return 2.*gamma_p*sqrAbs(psi(0));}; }



qbit::multidiagonal::MultiDiagonal qbit::multidiagonal::splus()
{
  MultiDiagonal res;
  res.diagonals.emplace( MultiDiagonal::Offsets{-1}, multiarray::fromStorage( [=] {auto res{noInit(1)}; res[0]=1; return res;} () ) ); // σ⁺ = |1⟩⟨0|: lower
  return res;
}

qbit::multidiagonal::MultiDiagonal qbit::multidiagonal::sminus() {return hermitianConjugateOf(splus());}

qbit::multidiagonal::MultiDiagonal qbit::multidiagonal::sx() {return (splus()+sminus())/2;}
qbit::multidiagonal::MultiDiagonal qbit::multidiagonal::sy() {return (splus()-sminus())/2i;}

qbit::multidiagonal::MultiDiagonal qbit::multidiagonal::sz()
{
  MultiDiagonal res;
  res.diagonals.emplace( MultiDiagonal::Offsets{0}, multiarray::fromStorage( [=] {auto res{noInit(2)}; res[0]=-.5; res[1]=.5; return res;} () ) ); // as spin ½: |1⟩ is up (#87)
  return res;
}


StateVector<1> qbit::init(dcomp psi1)
{
  StateVector<1> res({2},noInit);
  res(0)=sqrt(1-sqrAbs(psi1)); res(1)=psi1;
  return res;
}
