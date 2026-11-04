// Copyright András Vukics 2026. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
// Slicing along retained axes: the order of the slices (#90), and views with an offset of their own (#34).
#include "Composite.h"
#include "Mode.h"

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <numeric>

using namespace cppqedutils;

namespace {

/// Over a rank-4 array whose elements hold their own addresses: the slices start at increasing addresses, and together they cover every element once
template <auto ra>
void checkMemoryOrder(Extents<4> ext)
{
  std::vector<size_t> data(multiarray::calculateExtent(ext));
  std::iota(data.begin(),data.end(),0uz);
  const MultiArrayView<const size_t,4> view{ext,multiarray::calculateStrides(ext),0,std::span<const size_t>{data}};

  std::vector<size_t> starts, addresses;
  for (auto&& slice : sliceRange<ra>(view)) {
    starts.push_back(slice.dataView[slice.offset]);
    for (Extents<std::size(ra)> idx{}; idx[0]!=slice.extents[0]; incrementMultiIndex(idx,slice.extents)) addresses.push_back(slice(idx));
  }

  std::vector<size_t> ascending{starts}; std::ranges::sort(ascending);
  CHECK(starts==ascending);
  std::ranges::sort(addresses);
  CHECK(addresses==data);
}

/// Rank-3 states as the slices at index 1 of the last axis of rank-4 arrays, so that their offset is not 0
struct SliceOfLargerState
{
  static constexpr Extents<4> ext4{3,4,5,2};
  static constexpr Extents<3> ext3{3,4,5};
  const Extents<4> strides4{multiarray::calculateStrides(ext4)};
  const Extents<3> strides3{strides4[0],strides4[1],strides4[2]};
  const size_t offset{strides4[3]};

  MultiArray<dcomp,4> psi{ext4,[] (size_t e) {
    std::vector<dcomp> res(e); for (size_t k=0; k<e; ++k) res[k]=dcomp{std::cos(0.3*k),std::sin(0.7*k)}; return res;}};

  MultiArrayConstView<dcomp,3> psiSlice() const {return {ext3,strides3,offset,psi.dataView};}
};

} // anonymous


TEST_CASE("calculateSlicesOffsets: the slices come in memory order, first dummy axis fastest — #90", "[slicing]")
{
  // extents {2,3,4} have strides {1,2,6}; with axis 1 retained, the dummy axes 0 and 2 have strides 1 and 6
  CHECK(cso<1>(Extents<3>{2,3,4})==std::vector<size_t>{0,1,6,7,12,13,18,19});

  const Extents<4> ext{2,3,4,5};
  checkMemoryOrder<retainedAxes<0>>(ext);
  checkMemoryOrder<retainedAxes<2>>(ext);
  checkMemoryOrder<retainedAxes<3>>(ext);
  checkMemoryOrder<retainedAxes<1,3>>(ext);
  checkMemoryOrder<retainedAxes<3,0>>(ext);
  checkMemoryOrder<retainedAxes<0,1,2>>(ext);
}


TEST_CASE("sliceRange: the first slice of a view with an offset starts at that offset — #34", "[slicing]")
{
  // a 2×3 view at offset 1 of the storage 0…6: element (i,j) holds 1+i+2j
  const std::vector<double> data{0,1,2,3,4,5,6};
  const MultiArrayView<const double,2> view{{2,3},{1,2},1,std::span<const double>{data}};

  // a complete loop also used to read one past the end of the offsets (ASan: heap-buffer-overflow)
  std::vector<std::vector<double>> slices;
  for (auto&& slice : sliceRange<retainedAxes<1>>(view)) slices.push_back({slice(0),slice(1),slice(2)});

  CHECK(slices==std::vector<std::vector<double>>{{1,3,5},{2,4,6}});
}


TEST_CASE("Broadcaster: a term applied to a slice of a larger state acts on that slice only — #34, #90", "[slicing][composite]")
{
  const SliceOfLargerState s;
  const auto A{mode::aOp(4)+2.*mode::aDagOp(4)+mode::nOp(4)};

  MultiArray<dcomp,4> dpsidt{s.ext4,multiarray::zeroInit<dcomp>}, expected{s.ext4,multiarray::zeroInit<dcomp>};

  Broadcaster<3,1> b{s.ext3};
  b(A)(0.,s.psiSlice(),MultiArrayView<dcomp,3>{s.ext3,s.strides3,s.offset,std::span<dcomp>{dpsidt.dataStorage()}},0.);

  // the same operator applied to each slice, built by hand
  for (size_t i0=0; i0<s.ext3[0]; ++i0) for (size_t i2=0; i2<s.ext3[2]; ++i2) {
    const size_t o=s.offset+i0*s.strides4[0]+i2*s.strides4[2];
    A(0.,MultiArrayConstView<dcomp,1>{{s.ext3[1]},{s.strides4[1]},o,s.psi.dataView},
      MultiArrayView<dcomp,1>{{s.ext3[1]},{s.strides4[1]},o,std::span<dcomp>{expected.dataStorage()}});
  }

  CHECK(dpsidt.dataStorage()==expected.dataStorage());
}


TEST_CASE("partialTrace of a slice of a larger state — #34", "[slicing]")
{
  const SliceOfLargerState s;

  // population of level 1 on axis 1
  const double p1=quantumdata::partialTrace<retainedAxes<1>>(quantumdata::LDO<quantumdata::StateVector,3>{s.psiSlice()}, cso<1>(s.ext3),
                                                              [] (quantumdata::LDO<quantumdata::StateVector,1> l) {return quantumdata::_(l,1uz);},
                                                              std::plus{});

  double expected=0;
  for (size_t i0=0; i0<s.ext3[0]; ++i0) for (size_t i2=0; i2<s.ext3[2]; ++i2) expected+=sqrAbs(s.psi(i0,1,i2,1));

  CHECK(std::abs(p1-expected)<1e-12*expected);
}
