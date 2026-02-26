// Copyright (C) 2018 Fabio Rosa
//
// This file is inspired to the valarray non-member begin/end overload
// The intent of this function is to allow range for loops to work with Eigen lib arrays, not to provide container semantics.
//
// EigenUtils is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// EigenUtils is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You can check a copy of the GNU General Public License at
// <http://www.gnu.org/licenses/>.

#include "Eigen/Dense"

#ifndef ___EigenIteratorHelper___
#define ___EigenIteratorHelper___

// MSVC tip: use /Zc:__cplusplus compiler language option
#if __cplusplus >= 201103L
  /**
   *  @brief  Return an iterator pointing to the first element of
   *          the valarray.
   *  @param  __va  valarray.
   */
  template<class _Tp>
    inline _Tp*
    begin(Eigen::Array<_Tp, -1, 1>& __va)
    { return std::addressof(__va(0)); }

  /**
   *  @brief  Return an iterator pointing to the first element of
   *          the const valarray.
   *  @param  __va  valarray.
   */
  template<class _Tp>
    inline const _Tp*
    begin(const Eigen::Array<_Tp, -1, 1>& __va)
    { return std::addressof(__va(0)); }

  /**
   *  @brief  Return an iterator pointing to one past the last element of
   *          the valarray.
   *  @param  __va  valarray.
   */
  template<class _Tp>
    inline _Tp*
    end(Eigen::Array<_Tp, -1, 1>& __va)
    { return std::addressof(__va(0)) + __va.size(); }

  /**
   *  @brief  Return an iterator pointing to one past the last element of
   *          the const valarray.
   *  @param  __va  valarray.
   */
  template<class _Tp>
    inline const _Tp*
    end(const Eigen::Array<_Tp, -1, 1>& __va)
    { return std::addressof(__va(0)) + __va.size(); }

#endif

#endif