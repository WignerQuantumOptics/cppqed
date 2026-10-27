// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#pragma once

#include "FrameworkTraits.h"

#include "popl.hpp"

#include <list>

inline std::string parsedCommandLine="";

popl::OptionParser optionParser(std::string pre="", std::string post="");

void parse(popl::OptionParser&, int argc, const char* const argv[]);


namespace parameters {

using namespace cppqedutils;

using to_json_converter=std::list<std::function<void(json::object&)>> ;

template <typename T>
void add(popl::OptionParser& op, std::string option, std::string description, T defaultValue, T& binding)
{
  op.add<popl::Value<T>>("",option,description,defaultValue,&binding);
}


template <typename T>
void add(popl::OptionParser& op, std::string option, std::string mod, std::string description, T defaultValue, T& binding)
{
  add(op,option+mod,description,defaultValue,binding);
}


template <typename T>
void add(popl::OptionParser& op, to_json_converter& tjc, std::string option, std::string mod, std::string description, T defaultValue, T& binding)
{
  add(op,option,mod,description,defaultValue,binding) ;
  tjc.push_back( [o=option,&b=binding] (json::object& v) {v.emplace(o,json::value_from(b));} );
}


// inline void add(popl::OptionParser& op, std::string option, std::string description, bool* binding)
// {
//   op.add<popl::Switch>("",option,description,binding);
// }
//
//
// inline void add(popl::OptionParser& op, std::string option, std::string mod, std::string description, bool* binding)
// {
//   add(op,option+mod,description,binding);
// }



/// TODO: insert this more directly into popl
inline void addTitle(popl::OptionParser& op, std::string title, std::string mod = "")
{
  op.add<popl::Value<std::string>>("","\n### "+title+mod,"###","###");
}


inline void addTitle(popl::OptionParser& op, to_json_converter& tjc, std::string title, std::string mod = "")
{
  addTitle(op,title,mod);
  tjc.push_back( [t=title] (json::object& v) {v.emplace("######",t);} );
}


template <typename T, std::convertible_to<T> U>
auto _(std::string option, std::string description, const U& defaultValue, T& binding)
{
  return std::tuple<std::string,std::string,T,T&>{option,description,defaultValue,binding};
}


inline const struct NoJSON {} noJSON;

template <typename T, std::convertible_to<T> U>
auto _(std::string option, std::string description, const U& defaultValue, T& binding, NoJSON)
{
  return std::tuple<std::string,std::string,T,T&,NoJSON>{option,description,defaultValue,binding,noJSON};
}


template <typename T, std::convertible_to<T> U>
auto _(to_json_converter& tjc, std::string option, std::string description, const U& defaultValue, T& binding)
{
  return std::tuple<to_json_converter&,std::string,std::string,T,T&>{tjc,option,description,defaultValue,binding};
}



/// `T` is either a string (title), or a tuple
template <typename... T> requires ( ... && ( decltype( ::cppqedutils::multilambda {
  [] <typename U, typename V> (const std::tuple<to_json_converter&,std::string,std::string,U,V&> &) { return std::is_convertible<U,V>{}; },
  [] <typename U, typename V> (const std::tuple<std::string,std::string,U,V&> &) { return std::is_convertible<U,V>{}; },
  [] <typename U> (const U &) { return std::is_convertible<U,std::string>{}; }
  } (std::declval<T>()))::value ) )
void add_dispatch(std::string mod, popl::OptionParser& op, const T&... t)
{
  ::cppqedutils::multilambda worker {
    [&op,mod] <typename U, typename V> (const std::tuple<to_json_converter&,std::string,std::string,U,V&> & t ) {
      add(op,std::get<0>(t),std::get<1>(t),mod,std::get<2>(t),std::get<3>(t),std::get<4>(t));
    },
    [&op,mod] <typename U, typename V> (const std::tuple<std::string,std::string,U,V&> & t ) {
      add(op,std::get<0>(t),mod,std::get<1>(t),std::get<2>(t),std::get<3>(t));
    },
    [&op,mod] (const std::string& title) {addTitle(op,title,mod);}
  };

  (worker(t), ...);
}



template <typename... T> requires ( ... && ( decltype( ::cppqedutils::multilambda {
  [] <typename U, typename V> (const std::tuple<std::string,std::string,U,V&> &) { return std::is_convertible<U,V>{}; },
  [] <typename U, typename V> (const std::tuple<std::string,std::string,U,V&,NoJSON> &) { return std::is_convertible<U,V>{}; },
  [] <typename U> (const U &) { return std::is_convertible<U,std::string>{}; }
  } (std::declval<T>()))::value ) )
void add_dispatch(std::string mod, popl::OptionParser& op, to_json_converter& tjc, const T&... t)
{
  ::cppqedutils::multilambda worker {
    [&op,&tjc,mod] <typename U, typename V> (const std::tuple<std::string,std::string,U,V&> & t ) {
      add(op,tjc,std::get<0>(t),mod,std::get<1>(t),std::get<2>(t),std::get<3>(t));
    },
    [&op,&tjc,mod] <typename U, typename V> (const std::tuple<std::string,std::string,U,V&,NoJSON> & t ) {
      add(op,std::get<0>(t),mod,std::get<1>(t),std::get<2>(t),std::get<3>(t));
    },
    [&op,&tjc,mod] (const std::string& title) {addTitle(op,tjc,title,mod);}
  };

  (worker(t), ...);
}


struct JSONizable
{
  to_json_converter tjc;

  /// JSONize
  json::object jsonize() const
  {
    json::object res{};
    for (const auto& f : tjc) f(res);
    return res;
  }

};


} // parameters


void add(std::string mod, popl::OptionParser& op, parameters::to_json_converter& tjc, const auto&... t) {parameters::add_dispatch(mod,op,tjc,t...);}

void add(popl::OptionParser& op, parameters::to_json_converter& tjc, const auto&... t) {parameters::add_dispatch("",op,tjc,t...);}


void add(std::string mod, popl::OptionParser& op, const auto&... t) {parameters::add_dispatch(mod,op,t...);}

void add(popl::OptionParser& op, const auto&... t) {parameters::add_dispatch("",op,t...);}

