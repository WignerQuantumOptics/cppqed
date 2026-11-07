// Copyright András Vukics 2006–2026. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "Spin.h"


using namespace quantumdynamics;

namespace md = ::quantumoperator::multidiagonal;


void decideDimension(size_t twoS, size_t& dim) {dim = (dim && dim<twoS+1) ? dim : twoS+1;}


spin::MultiDiagonal spin::splus(size_t twoS, size_t dim)
{
  decideDimension(twoS,dim);
  std::vector<dcomp> values(dim-1);
  for (size_t i=0; i<values.size(); ++i) values[i]=sqrt((twoS-i)*(i+1.)); // i=m+s (m is the magnetic quantum number)
  return md::diagonal(-1,std::move(values)); // S⁺ raises m: lower
}


spin::MultiDiagonal spin::sminus(size_t twoS, size_t dim) {return hermitianConjugateOf(splus(twoS,dim));}

spin::MultiDiagonal spin::sx(size_t twoS, size_t dim) {return (splus(twoS,dim)+sminus(twoS,dim))/2;}

spin::MultiDiagonal spin::sy(size_t twoS, size_t dim) {return (splus(twoS,dim)-sminus(twoS,dim))/2i;}

spin::MultiDiagonal spin::sz(size_t twoS, size_t dim)
{
  decideDimension(twoS,dim);
  std::vector<dcomp> values(dim);
  for (size_t i=0; i<values.size(); ++i) values[i]=i-twoS/2.; // i=m+s (m is the magnetic quantum number)
  return md::diagonal(0,std::move(values));
}



/*
namespace spin {

namespace {

typedef Tridiagonal::Diagonal Diagonal;

}




Diagonal mainDiagonal(size_t dim, size_t twoS, dcomp z)
{
  Diagonal diagonal(dim);
  diagonal=blitz::tensor::;
  return Diagonal(z*diagonal);
}




Tridiagonal sn(Ptr spin)
{
  return sn(spin->getDimension(),spin->getTwoS(),spin->get_z(),bool(dynamic_pointer_cast<const quantumdynamics::FreeExact<false>>(spin)),spin->getTheta(),spin->getPhi());
} 


Tridiagonal sz(size_t dim, size_t twoS)
{
  Diagonal diagonal(dim);
  return Tridiagonal(diagonal=blitz::tensor::i-twoS/2.);
}

Tridiagonal sz(Ptr spin) {return sz(spin->getDimension(),spin->getTwoS());}


Tridiagonal splus(Ptr spin) {return splus(spin->getDimension(),spin->getTwoS(),spin->get_z(),bool(dynamic_pointer_cast<const quantumdynamics::FreeExact<false>>(spin)));}



Pars::Pars(parameters::Table& p, const std::string& mod)
  : twoS(p.addTitle("Spin",mod).add<size_t>("twoS",mod,"2*s, the size of the spin (dimension: twoS+1 or spinDim)",1)),
    dim (p.add<size_t>("spinDim",mod,"the dimension of the truncated spin Hilbert space",0)),
    theta(p.add("theta",mod,"Spin orientation inclination",0.)),
    phi  (p.add("phi"  ,mod,"Spin orientation azimuth    ",0.)),
    omega(p.add("omega",mod,"Spin precession frequency",1.)),
    gamma(p.add("gamma",mod,"Spin decay rate"          ,1.))
{}



} // spin


SpinBase::SpinBase(size_t twoS, double theta, double phi, double omega, double gamma, size_t dim) 
  : Free(spin::decideDimension(twoS,dim),{RF{"omega",omega,1},RF{"gamma",gamma,1}}),
    quantumdynamics::ElementAveraged<1>("Spin",{"<sz>","<sz^2>","real(<s^+>)","imag(\")","|Psi(dim-1)|^2"}), twoS_(twoS), theta_(theta), phi_(phi), omega_(omega), gamma_(gamma), s_(twoS/2.)
{
  getParsStream()<<"Spin "<<s_<<endl;
  getParsStream()<<"theta="<<theta_<<" phi="<<phi_<<endl;
}



const spin::Averages SpinBase::average_v(quantumdynamics::NoTime, const LazyDensityOperator& matrix) const
{
  auto averages(initializedAverages());

  for (int n=0; n<matrix.getDimension(); n++) {

    double diag=matrix(n);
    averages(0)+=(n-s_)*diag;
    averages(1)+=sqr(n-s_)*diag;

    if (n<matrix.getDimension()-1) {
      dcomp temp(sqrt((twoS_-n)*(n+1))*matrix(n)(n+1));
      averages(2)+=real(temp);
      averages(3)+=imag(temp);
    }

  }

  averages(4)=matrix(matrix.getDimension()-1);

  return averages;

}


void SpinBase::process_v(spin::Averages& averages) const
{
  averages(1)-=sqr(averages(0));
}


void Spin::updateU(quantumdynamics::OneTime dtDid) const
{
  Diagonal& factors(getDiagonal());
  for (int i=0; i<factors.size(); i++)
    factors(i)=exp(-dtDid*i*get_z());
}


void spin::Liouvillian::doActWithJ(quantumdynamics::NoTime, quantumdynamics::freesystem::StateVectorLow& psi) const
{
  double fact=sqrt(2.*gamma_);
  int ubound=psi.ubound(0);
  for (int n=0; n<ubound; ++n)
    psi(n)=fact*sqrt((n+1)*(twoS_-n))*psi(n+1);
  psi(ubound)=0;
}


double spin::Liouvillian::rate(quantumdynamics::NoTime, const quantumdynamics::freesystem::LazyDensityOperator& matrix) const
{
  double res=0;
  for (size_t n=1; n<matrix.getDimension(); ++n)
    res+=2.*gamma_*n*(twoS_-n+1)*matrix(n);
  return res;
}

*/
