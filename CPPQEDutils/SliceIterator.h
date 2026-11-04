// Copyright András Vukics 2022–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#pragma once


#include "FrameworkTraits.h"

#include <multidiagonal/MultiArray.h>

#include <boost/operators.hpp>



namespace cppqedutils {
  

template<size_t... i>
constexpr std::array retainedAxes{i...};

template <auto retainedAxes, size_t RANK>
constexpr std::array extendedAxes = concatenate(retainedAxes, [rA=retainedAxes] mutable {for (size_t& i : rA) i+=RANK; return rA;} () );


namespace slice_iterator {


/// Checking the consistency of template arguments for use in slicing
/**
 * - size of `RetainedAxes` must not be larger than `RANK`
 * - `RetainedAxes` must not contain duplicated elements
 * - all elements of `RetainedAxes` must be smaller than `RANK`
 */
template <auto retainedAxes, size_t RANK>
constexpr bool consistent = (
  RANK > 1 && std::size(retainedAxes) <= RANK &&
  (*std::ranges::max_element(retainedAxes) < RANK) &&
  [us(retainedAxes)] mutable {
    std::ranges::sort(us);
    const auto [first, last] = std::ranges::unique(us);
    return (first==last);
} () );


/// Filters in the indices corresponding to a subsystem
/**
 * When the size of `retainedAxes` equals `RANK`, this is a transposition
 * 
 * TODO: accept several arrays at the same time (since this is how it is used by clients)
 */
template <auto retainedAxes, size_t RANK>
auto filterIn(Extents<RANK> idx) requires ( std::size(retainedAxes) <= RANK )
{
  Extents<std::size(retainedAxes)> res;
  std::ranges::fold_left(retainedAxes, res.begin(), [&] (auto iterator, auto ind) {*iterator=idx[ind]; return ++iterator;});
  return res;
}

  
/// Filters out the indices corresponding to a subsystem (specified by `retainedAxes`)
/** TODO: accept several arrays at the same time (since this is how it is used by clients) */
template <auto retainedAxes, size_t RANK>
auto filterOut(Extents<RANK> idx) requires ( std::size(retainedAxes) <= RANK )
{
  Extents<RANK-std::size(retainedAxes)> res;
  
  std::ranges::fold_left(std::views::iota(0ul,RANK), res.begin(), [&] (auto iterator, auto ind) {
    if ( ! std::ranges::contains(retainedAxes, ind) ) *iterator++ = idx[ind];
    return iterator;
  } );
  
  return res;

}


/// Calculates the slices offsets purely from the (dummy) extents and strides
/**
 * The slices come in the order of storage, first dummy axis fastest, so that consecutive slices lie close in memory (#90).
 *
 * Any indexing offset that the MultiArrayView might itself have is added at the point of indexing/slicing
 */
template <auto retainedAxes, size_t RANK>
auto calculateSlicesOffsets(Extents<RANK> extents, Extents<RANK> strides)
{
  Extents<RANK-std::size(retainedAxes)>
    dummyExtents{filterOut<retainedAxes>(extents)},
    dummyStrides{filterOut<retainedAxes>(strides)},
    idx; idx.fill(0);

  std::vector<size_t> res(multiarray::calculateExtent(dummyExtents));

  for (size_t& offset : res) {
    offset = std::ranges::fold_left_first( std::views::zip(idx,dummyStrides) | std::views::transform( [] (auto ids) {return get<0>(ids)*get<1>(ids);} ) , std::plus{}).value_or(0uz);
    // first index fastest, unlike incrementMultiIndex
    for (size_t k=0; k<idx.size() && ++idx[k]==dummyExtents[k]; ++k) idx[k]=0;
  }

  return res;
}


} // slice_iterator



/// To be used by Composites
template <auto retainedAxes, size_t RANK>
auto calculateSlicesOffsets(Extents<RANK> extents)
{
  return slice_iterator::calculateSlicesOffsets<retainedAxes>(extents,multiarray::calculateStrides(extents));
}


/// TODO: how to integrate this with the above overloads of calculateSlicesOffsets
template <size_t... i, size_t RANK> requires (sizeof...(i)<RANK)
auto cso(Extents<RANK> extents)
{
  return calculateSlicesOffsets<retainedAxes<i...>>(extents);
}


template <auto retainedAxes, typename T, size_t RANK>
auto sliceRangeSimple(MultiArrayView<T,RANK> mav, const std::vector<size_t>& offsets)
requires slice_iterator::consistent<retainedAxes,RANK>
{
  using SliceType=MultiArrayView<T,std::size(retainedAxes)>;
  
  std::vector<SliceType> res(offsets.size());
  
  for (auto&& [slice,sliceOffset] : std::views::zip(res,offsets) )
    slice=SliceType{slice_iterator::filterIn<retainedAxes>(mav.extents),slice_iterator::filterIn<retainedAxes>(mav.strides),mav.offset+sliceOffset,mav.dataView};
  
  return res;
}


template <auto retainedAxes, typename T, size_t RANK>
auto sliceRangeSimple(MultiArrayView<T,RANK> mav)
{
  return sliceRangeSimple<retainedAxes>(mav,slice_iterator::calculateSlicesOffsets<retainedAxes>(mav.extents,mav.strides));
}


namespace slice_iterator {

template <typename O>
concept slice_iteration_offsets = std::is_same_v<O,std::span<const size_t>> || std::is_same_v<O,std::vector<size_t>> ;

} // slice_iterator


/// TODO: think over semantics here, can MultiArrayView be the reference type?
template <slice_iterator::slice_iteration_offsets Offsets, typename T, size_t RRANK> // RRANK stands for retained rank
class SliceIterator : public boost::forward_iterator_helper<SliceIterator<Offsets,T,RRANK>,MultiArrayView<T,RRANK>>
/// TODO: Unfortunately, std::iterator has been deprecated as of C++17, and there doesn’t seem to be a replacement in STL
{
public:
  using member_iterator = std::conditional_t<
    std::is_same_v<Offsets,std::span<const size_t>>,
    std::span<const size_t>::iterator,
    std::vector<size_t>::const_iterator>;
  
  SliceIterator(const SliceIterator&) = default; SliceIterator(SliceIterator&&) = default; SliceIterator() = default;
  SliceIterator& operator=(const SliceIterator&) = default; SliceIterator& operator=(SliceIterator&&) = default;
  
  SliceIterator(Extents<RRANK> e, Extents<RRANK> s, size_t originalOffset, member_iterator i, std::span<T> dv) :
    iter{i}, originalOffset_{originalOffset}, mav_{e,s,originalOffset,dv} {}

  /// The slice’s offset is set here, since on construction and increment `iter` may be the end (#34)
  auto& operator*() const {mav_.offset=*iter+originalOffset_; return mav_;}

  SliceIterator& operator++() {++iter; return *this;}
  
  member_iterator iter;

  friend bool operator==(const SliceIterator& i0, const member_iterator& i1) {return i0.iter==i1;}

  friend bool operator==(const SliceIterator& i0, const SliceIterator& i1) {return i0.iter==i1.iter;}
  
private:
  size_t originalOffset_;
  
  mutable MultiArrayView<T,RRANK> mav_; // the reason why we store a full MultiArrayView is simply that operator* can return a reference
  
};



template <typename Offsets, typename T, size_t RRANK> // RRANK stands for retained rank
class SliceRangeBase
{
public:
  using iterator=SliceIterator<Offsets,T,RRANK>;
  
  SliceRangeBase(const SliceRangeBase&) = delete; SliceRangeBase& operator=(const SliceRangeBase&) = delete; 
  SliceRangeBase(SliceRangeBase&&) = default; SliceRangeBase& operator=(SliceRangeBase&&) = default; // SliceRangeBase() = default;
  
  
/*  SliceRangeBase(Extents<RRANK> e, Extents<RRANK> s, Offsets os, std::span<T> dv) requires std::is_same_v<Offsets,std::span<const size_t>>
    : begin{e,s,os.begin(),dv}, end{os.end()}, offsets_{os} {}*/

  SliceRangeBase(Extents<RRANK> e, Extents<RRANK> s, size_t originalOffset, const Offsets& os, std::span<T> dv) // requires std::is_same_v<Offsets,std::vector<size_t>>
    : offsets_{os}, b_{e,s,originalOffset,offsets_.begin(),dv}, e_{offsets_.end()} {}
    
  auto begin() const {return b_;}
  auto end() const {return e_;}
    
private:
  Offsets offsets_;

  iterator b_;
  
  typename iterator::member_iterator e_;
  
};


/**
 * For the case when Offsets are stored outside of the class.
 * E.g. Composite might calculate and store the offsets for all its subsystems
 */
template <typename T, size_t RRANK>
using SliceRangeReferencing = SliceRangeBase<std::span<const size_t>,T,RRANK>;


/** For the case when Offsets are stored within the class. */
template <typename T, size_t RRANK> // RRANK stands for retained rank
using SliceRangeOwning = SliceRangeBase<std::vector<size_t>,T,RRANK>;



template <auto retainedAxes, typename T, size_t RANK>
auto sliceRange(MultiArrayView<T,RANK> mav, const std::vector<size_t>& offsets) requires slice_iterator::consistent<retainedAxes,RANK>
{
  return SliceRangeReferencing<T,std::size(retainedAxes)>{
    slice_iterator::filterIn<retainedAxes>(mav.extents),
    slice_iterator::filterIn<retainedAxes>(mav.strides),
    mav.offset,
    std::span<const size_t>(offsets),
    mav.dataView};
}


template <auto retainedAxes, typename T, size_t RANK>
auto sliceRange(MultiArrayView<T,RANK> mav, std::vector<size_t>&& offsets) requires slice_iterator::consistent<retainedAxes,RANK>
{
  return SliceRangeOwning<T,std::size(retainedAxes)>{
    slice_iterator::filterIn<retainedAxes>(mav.extents),
    slice_iterator::filterIn<retainedAxes>(mav.strides),
    mav.offset,
    std::move(offsets),
    mav.dataView};
}


template <auto retainedAxes, typename T, size_t RANK>
auto sliceRange(MultiArrayView<T,RANK> mav) requires slice_iterator::consistent<retainedAxes,RANK>
{
  return sliceRange<retainedAxes>(mav,slice_iterator::calculateSlicesOffsets<retainedAxes>(mav.extents,mav.strides));
}


template <auto axes, typename T, size_t RANK>
MultiArrayView<T,RANK> transpose(MultiArrayView<T,RANK> mav) requires ( std::size(axes) == RANK )
{
  mav.extents=filterIn<axes>(mav.extents); mav.strides=filterIn<axes>(mav.strides);
  return mav;
}



template <auto retainedAxes, typename F, typename... ARGS>
void broadcastFor(F&& f, ARGS&&... args)
{
  auto arrayArgs=hana::filter(hana::make_tuple(args...),
                              multilambda{
                                [] <typename T, size_t RANK> (MultiArrayView<T,RANK>) {return std::true_type{};},
                                [] (const auto&) {return std::false_type{};}
                              });

  // check extents

  // create offset vector

  // create ranges
}



namespace multi_index_range_solution_by_chatgpt {


template <size_t RANK>
struct MultiIndexIterator {
  using ExtentsType = Extents<RANK>;
  using value_type = ExtentsType;
  using difference_type = std::ptrdiff_t;
  using pointer = ExtentsType*;
  using reference = ExtentsType&;
  using iterator_category = std::forward_iterator_tag;

  mutable ExtentsType idx;
  ExtentsType extents;

  MultiIndexIterator(ExtentsType idx, ExtentsType extents)
    : idx(std::move(idx)), extents(std::move(extents))
  {}

  reference operator*() const { return idx; }
  pointer operator->() const { return &idx; }

  MultiIndexIterator& operator++() {
    increment(hana::llong_c<RANK - 1>);
    return *this;
  }

  MultiIndexIterator operator++(int) {
    MultiIndexIterator prev = *this;
    ++(*this);
    return prev;
  }

  friend auto operator<=>(const MultiIndexIterator& lhs, const MultiIndexIterator& rhs) = default;

  friend bool operator!=(const MultiIndexIterator& lhs, size_t end) {return lhs.idx[0]!=end;}

private:
  void increment(auto n) {
    using namespace hana::literals;
    if constexpr (n != 0_c) {
      if (idx[n] == extents[n] - 1) {
        idx[n] = 0;
        increment(n - 1_c);
      } else {
        ++idx[n];
      }
    } else {
      ++idx[0];
    }
  }
};

template <size_t RANK>
struct MultiIndexRange {
  using ExtentsType = Extents<RANK>;

  ExtentsType extents;

  explicit MultiIndexRange(ExtentsType extents)
    : extents(std::move(extents))
  {}

  MultiIndexIterator<RANK> begin() const {
    ExtentsType idx{};
    return MultiIndexIterator<RANK>(idx, extents);
  }

  size_t end() const { return extents[0]; }
};


// It doesn’t seem to be a valid range somehow, so for example this doesn’t compile
// static_assert( std::ranges::range<MultiIndexRange<3>> );
// This is ok though:
static_assert( std::input_iterator<MultiIndexIterator<3>> );


} // multi_index_range_solution_by_chatgpt

} // cppqedutils

