// Copyright (c) 2018 by University Paris-Est Marne-la-Vallee
// GeometricExplicit.hpp
// This file is part of the Garamon for e2ga.
// Authors: Stephane Breuils and Vincent Nozick
// Conctat: vincent.nozick@u-pem.fr
//
// Licence MIT
// A a copy of the MIT License is given along with this program

/// \file GeometricExplicit.hpp
/// \author Stephane Breuils, Vincent Nozick
/// \brief Explicit precomputed per grades geometric products of e2ga.


#ifndef E2GA_GEOMETRIC_PRODUCT_EXPLICIT_HPP__
#define E2GA_GEOMETRIC_PRODUCT_EXPLICIT_HPP__
#pragma once

#include <Eigen/Core>

#include "e2ga/Mvec.hpp"
#include "e2ga/Constants.hpp"


/*!
 * @namespace e2ga
 */
namespace e2ga {
    template<typename T> class Mvec;

    
    template<typename T>
	std::array<std::array<std::array<std::function<void(const Eigen::Matrix<T, Eigen::Dynamic, 1> & , const Eigen::Matrix<T, Eigen::Dynamic, 1> & , Eigen::Matrix<T, Eigen::Dynamic, 1>&)>, 3>, 3>, 3> geometricFunctionsContainer =
	{{
		{{
			{{{},{},{}}},
			{{{},{},{}}},
			{{{},{},{}}}
		}},
		{{
			{{{},{},{}}},
			{{{},{},{}}},
			{{{},{},{}}}
		}},
		{{
			{{{},{},{}}},
			{{{},{},{}}},
			{{{},{},{}}}
		}}
	}};

}/// End of Namespace

#endif // E2GA_GEOMETRIC_PRODUCT_EXPLICIT_HPP__