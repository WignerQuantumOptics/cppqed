// Copyright András Vukics 2026. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
// The subsystem operators ported to signed offsets (#23), against dense matrices built independently.
#include "Mode.h"
#include "Qbit.h"
#include "Spin.h"

#include <catch2/catch_test_macros.hpp>

#include <cmath>

using namespace cppqedutils;
using quantumoperator::MultiDiagonal;

namespace {

constexpr double tol = 1e-12;

bool approxEqual(const CMatrix& A, const CMatrix& B) {return A.rows()==B.rows() && A.cols()==B.cols() && (A-B).norm() < tol;}

/// The oracle: probe with unit vectors (as in the cppqed-multidiagonal test suite).
CMatrix densify(const MultiDiagonal<1>& md)
{
  const size_t D = calculateAndCheckDimensions(md)[0];
  CMatrix res{CMatrix::Zero(D,D)};
  MultiArray<dcomp,1> psi{{D}, multiarray::zeroInit<dcomp>}, dpsidt{{D}, multiarray::zeroInit<dcomp>};
  for (size_t k=0; k<D; ++k) {
    psi(k)=1.;
    for (dcomp& v : dpsidt.dataStorage()) v=0.;
    md(0., psi, dpsidt);
    res.col(k) = vectorize(dpsidt);
    psi(k)=0.;
  }
  return res;
}

CMatrix commutator(const CMatrix& A, const CMatrix& B) {return A*B-B*A;}

} // anonymous


TEST_CASE("Mode: a, a† and n are the truncated ladder and number operators — #23", "[mode][operators]")
{
  const size_t cutoff=6;
  CMatrix a{CMatrix::Zero(cutoff,cutoff)}, n{CMatrix::Zero(cutoff,cutoff)};
  for (size_t k=0; k+1<cutoff; ++k) a(k,k+1)=std::sqrt(double(k+1));   // ⟨k|a|k+1⟩ = √(k+1)
  for (size_t k=0; k<cutoff; ++k) n(k,k)=double(k);

  CHECK(approxEqual(densify(mode::aOp(cutoff)), a));
  CHECK(approxEqual(densify(mode::aDagOp(cutoff)), a.adjoint()));
  CHECK(approxEqual(densify(mode::nOp(cutoff)), n));
}


TEST_CASE("Qbit: σ⁺ = |1⟩⟨0| raises, σ⁻ lowers, σx and σy as their combinations, σz = diag(−½,+½) — #23, #87", "[qbit][operators]")
{
  CMatrix sp{CMatrix::Zero(2,2)}, sz{CMatrix::Zero(2,2)}; sp(1,0)=1.; sz(0,0)=-.5; sz(1,1)=.5;
  const CMatrix sm=sp.adjoint();

  CHECK(approxEqual(densify(qbit::multidiagonal::splus()), sp));
  CHECK(approxEqual(densify(qbit::multidiagonal::sminus()), sm));
  CHECK(approxEqual(densify(qbit::multidiagonal::sx()), ((sp+sm)/2.).eval()));
  CHECK(approxEqual(densify(qbit::multidiagonal::sy()), ((sp-sm)/2i).eval()));
  CHECK(approxEqual(densify(qbit::multidiagonal::sz()), sz));
}


TEST_CASE("Qbit: [σz, σ⁺] = σ⁺, and the qubit is spin ½ — #87", "[qbit][spin][operators]")
{
  const CMatrix sz=densify(qbit::multidiagonal::sz()), sp=densify(qbit::multidiagonal::splus());
  CHECK(approxEqual(commutator(sz,sp), sp));
  CHECK(approxEqual(sz, densify(spin::sz(1))));
  CHECK(approxEqual(sp, densify(spin::splus(1))));
}


TEST_CASE("Spin: S⁺ and Sz match |s,m⟩ matrix elements; [Sz,S±] = ±S± — #23", "[spin][operators]")
{
  for (size_t twoS : {1uz, 2uz, 5uz}) {
    const size_t dim=twoS+1;
    const double s=twoS/2.;
    CMatrix sp{CMatrix::Zero(dim,dim)}, sz{CMatrix::Zero(dim,dim)};
    for (size_t i=0; i<dim; ++i) {
      const double m=i-s;                                                 // index i ↔ magnetic quantum number m
      sz(i,i)=m;
      if (i+1<dim) sp(i+1,i)=std::sqrt(s*(s+1)-m*(m+1));                  // S⁺|s,m⟩ = √(s(s+1)−m(m+1)) |s,m+1⟩
    }
    CAPTURE(twoS);
    CHECK(approxEqual(densify(spin::splus(twoS)), sp));
    CHECK(approxEqual(densify(spin::sminus(twoS)), sp.adjoint()));
    CHECK(approxEqual(densify(spin::sz(twoS)), sz));
    CHECK(approxEqual(commutator(sz,sp), sp));
  }
}


TEST_CASE("Spin: Holstein–Primakoff onto the mode, S_z = n − s and S⁻ = √(2s − a†a) a — #87", "[spin][mode][operators]")
{
  for (size_t twoS : {1uz, 2uz, 5uz}) {
    const size_t dim=twoS+1;                                              // the mode truncated at 2s+1 quanta
    const CMatrix a=densify(mode::aOp(dim)), n=densify(mode::nOp(dim)), one=CMatrix::Identity(dim,dim);
    const CMatrix root=(twoS*one-n).cwiseSqrt();                         // √(2s − a†a): diagonal, entries ≥ 0
    CAPTURE(twoS);
    CHECK(approxEqual(densify(spin::sz(twoS)), (n-(twoS/2.)*one).eval()));
    CHECK(approxEqual(densify(spin::sminus(twoS)), (root*a).eval()));
    CHECK(approxEqual(densify(spin::splus(twoS)), (a.adjoint()*root).eval()));
  }
}
