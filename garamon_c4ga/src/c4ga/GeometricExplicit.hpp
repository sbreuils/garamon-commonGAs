// Copyright (c) 2018 by University Paris-Est Marne-la-Vallee
// GeometricExplicit.hpp
// This file is part of the Garamon for c4ga.
// Authors: Stephane Breuils and Vincent Nozick
// Conctat: vincent.nozick@u-pem.fr
//
// Licence MIT
// A a copy of the MIT License is given along with this program

/// \file GeometricExplicit.hpp
/// \author Stephane Breuils, Vincent Nozick
/// \brief Explicit precomputed per grades geometric products of c4ga.


#ifndef C4GA_GEOMETRIC_PRODUCT_EXPLICIT_HPP__
#define C4GA_GEOMETRIC_PRODUCT_EXPLICIT_HPP__
#pragma once

#include <Eigen/Core>

#include "c4ga/Mvec.hpp"
#include "c4ga/Constants.hpp"


/*!
 * @namespace c4ga
 */
namespace c4ga {
    template<typename T> class Mvec;

    /// \brief Compute the geometric product between two homogeneous multivectors mv1 (grade 2) and mv2 (grade 2). 
	/// \tparam the type of value that we manipulate, either float or double or something else.
	/// \param mv1 - the first homogeneous multivector of grade 2 represented as an Eigen::VectorXd
	/// \param mv2 - the second homogeneous multivector of grade 2 represented as a Eigen::VectorXd
	/// \param mv3 - the result of mv1 mv2 whose grade is 2
	template<typename T>
	void geometric_2_2_2(const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv1, const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv2, Eigen::Matrix<T, Eigen::Dynamic, 1>& mv3){
		mv3.coeffRef(0) +=  mv1.coeff(0)*mv2.coeff(4) - mv1.coeff(1)*mv2.coeff(5) - mv1.coeff(2)*mv2.coeff(6) - mv1.coeff(3)*mv2.coeff(7) - mv1.coeff(4)*mv2.coeff(0) + mv1.coeff(5)*mv2.coeff(1) + mv1.coeff(6)*mv2.coeff(2) + mv1.coeff(7)*mv2.coeff(3);
		mv3.coeffRef(1) +=  mv1.coeff(0)*mv2.coeff(5) + mv1.coeff(1)*mv2.coeff(4) - mv1.coeff(2)*mv2.coeff(9) - mv1.coeff(3)*mv2.coeff(10) - mv1.coeff(4)*mv2.coeff(1) - mv1.coeff(5)*mv2.coeff(0) + mv1.coeff(9)*mv2.coeff(2) + mv1.coeff(10)*mv2.coeff(3);
		mv3.coeffRef(2) +=  mv1.coeff(0)*mv2.coeff(6) + mv1.coeff(1)*mv2.coeff(9) + mv1.coeff(2)*mv2.coeff(4) - mv1.coeff(3)*mv2.coeff(12) - mv1.coeff(4)*mv2.coeff(2) - mv1.coeff(6)*mv2.coeff(0) - mv1.coeff(9)*mv2.coeff(1) + mv1.coeff(12)*mv2.coeff(3);
		mv3.coeffRef(3) +=  mv1.coeff(0)*mv2.coeff(7) + mv1.coeff(1)*mv2.coeff(10) + mv1.coeff(2)*mv2.coeff(12) + mv1.coeff(3)*mv2.coeff(4) - mv1.coeff(4)*mv2.coeff(3) - mv1.coeff(7)*mv2.coeff(0) - mv1.coeff(10)*mv2.coeff(1) - mv1.coeff(12)*mv2.coeff(2);
		mv3.coeffRef(4) +=  mv1.coeff(0)*mv2.coeff(8) + mv1.coeff(1)*mv2.coeff(11) + mv1.coeff(2)*mv2.coeff(13) + mv1.coeff(3)*mv2.coeff(14) - mv1.coeff(8)*mv2.coeff(0) - mv1.coeff(11)*mv2.coeff(1) - mv1.coeff(13)*mv2.coeff(2) - mv1.coeff(14)*mv2.coeff(3);
		mv3.coeffRef(5) += -mv1.coeff(0)*mv2.coeff(11) + mv1.coeff(1)*mv2.coeff(8) - mv1.coeff(6)*mv2.coeff(9) - mv1.coeff(7)*mv2.coeff(10) - mv1.coeff(8)*mv2.coeff(1) + mv1.coeff(9)*mv2.coeff(6) + mv1.coeff(10)*mv2.coeff(7) + mv1.coeff(11)*mv2.coeff(0);
		mv3.coeffRef(6) += -mv1.coeff(0)*mv2.coeff(13) + mv1.coeff(2)*mv2.coeff(8) + mv1.coeff(5)*mv2.coeff(9) - mv1.coeff(7)*mv2.coeff(12) - mv1.coeff(8)*mv2.coeff(2) - mv1.coeff(9)*mv2.coeff(5) + mv1.coeff(12)*mv2.coeff(7) + mv1.coeff(13)*mv2.coeff(0);
		mv3.coeffRef(7) += -mv1.coeff(0)*mv2.coeff(14) + mv1.coeff(3)*mv2.coeff(8) + mv1.coeff(5)*mv2.coeff(10) + mv1.coeff(6)*mv2.coeff(12) - mv1.coeff(8)*mv2.coeff(3) - mv1.coeff(10)*mv2.coeff(5) - mv1.coeff(12)*mv2.coeff(6) + mv1.coeff(14)*mv2.coeff(0);
		mv3.coeffRef(8) +=  mv1.coeff(4)*mv2.coeff(8) + mv1.coeff(5)*mv2.coeff(11) + mv1.coeff(6)*mv2.coeff(13) + mv1.coeff(7)*mv2.coeff(14) - mv1.coeff(8)*mv2.coeff(4) - mv1.coeff(11)*mv2.coeff(5) - mv1.coeff(13)*mv2.coeff(6) - mv1.coeff(14)*mv2.coeff(7);
		mv3.coeffRef(9) += -mv1.coeff(1)*mv2.coeff(13) + mv1.coeff(2)*mv2.coeff(11) - mv1.coeff(5)*mv2.coeff(6) + mv1.coeff(6)*mv2.coeff(5) - mv1.coeff(10)*mv2.coeff(12) - mv1.coeff(11)*mv2.coeff(2) + mv1.coeff(12)*mv2.coeff(10) + mv1.coeff(13)*mv2.coeff(1);
		mv3.coeffRef(10) += -mv1.coeff(1)*mv2.coeff(14) + mv1.coeff(3)*mv2.coeff(11) - mv1.coeff(5)*mv2.coeff(7) + mv1.coeff(7)*mv2.coeff(5) + mv1.coeff(9)*mv2.coeff(12) - mv1.coeff(11)*mv2.coeff(3) - mv1.coeff(12)*mv2.coeff(9) + mv1.coeff(14)*mv2.coeff(1);
		mv3.coeffRef(11) +=  mv1.coeff(4)*mv2.coeff(11) - mv1.coeff(5)*mv2.coeff(8) + mv1.coeff(8)*mv2.coeff(5) + mv1.coeff(9)*mv2.coeff(13) + mv1.coeff(10)*mv2.coeff(14) - mv1.coeff(11)*mv2.coeff(4) - mv1.coeff(13)*mv2.coeff(9) - mv1.coeff(14)*mv2.coeff(10);
		mv3.coeffRef(12) += -mv1.coeff(2)*mv2.coeff(14) + mv1.coeff(3)*mv2.coeff(13) - mv1.coeff(6)*mv2.coeff(7) + mv1.coeff(7)*mv2.coeff(6) - mv1.coeff(9)*mv2.coeff(10) + mv1.coeff(10)*mv2.coeff(9) - mv1.coeff(13)*mv2.coeff(3) + mv1.coeff(14)*mv2.coeff(2);
		mv3.coeffRef(13) +=  mv1.coeff(4)*mv2.coeff(13) - mv1.coeff(6)*mv2.coeff(8) + mv1.coeff(8)*mv2.coeff(6) - mv1.coeff(9)*mv2.coeff(11) + mv1.coeff(11)*mv2.coeff(9) + mv1.coeff(12)*mv2.coeff(14) - mv1.coeff(13)*mv2.coeff(4) - mv1.coeff(14)*mv2.coeff(12);
		mv3.coeffRef(14) +=  mv1.coeff(4)*mv2.coeff(14) - mv1.coeff(7)*mv2.coeff(8) + mv1.coeff(8)*mv2.coeff(7) - mv1.coeff(10)*mv2.coeff(11) + mv1.coeff(11)*mv2.coeff(10) - mv1.coeff(12)*mv2.coeff(13) + mv1.coeff(13)*mv2.coeff(12) - mv1.coeff(14)*mv2.coeff(4);
	}


	/// \brief Compute the geometric product between two homogeneous multivectors mv1 (grade 2) and mv2 (grade 3). 
	/// \tparam the type of value that we manipulate, either float or double or something else.
	/// \param mv1 - the first homogeneous multivector of grade 2 represented as an Eigen::VectorXd
	/// \param mv2 - the second homogeneous multivector of grade 3 represented as a Eigen::VectorXd
	/// \param mv3 - the result of mv1 mv2 whose grade is 3
	template<typename T>
	void geometric_2_3_3(const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv1, const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv2, Eigen::Matrix<T, Eigen::Dynamic, 1>& mv3){
		mv3.coeffRef(0) += -mv1.coeff(0)*mv2.coeff(6) + mv1.coeff(1)*mv2.coeff(3) + mv1.coeff(2)*mv2.coeff(10) + mv1.coeff(3)*mv2.coeff(11) - mv1.coeff(4)*mv2.coeff(0) - mv1.coeff(6)*mv2.coeff(4) - mv1.coeff(7)*mv2.coeff(5) + mv1.coeff(9)*mv2.coeff(1) + mv1.coeff(10)*mv2.coeff(2);
		mv3.coeffRef(1) += -mv1.coeff(0)*mv2.coeff(8) - mv1.coeff(1)*mv2.coeff(10) + mv1.coeff(2)*mv2.coeff(3) + mv1.coeff(3)*mv2.coeff(13) - mv1.coeff(4)*mv2.coeff(1) + mv1.coeff(5)*mv2.coeff(4) - mv1.coeff(7)*mv2.coeff(7) - mv1.coeff(9)*mv2.coeff(0) + mv1.coeff(12)*mv2.coeff(2);
		mv3.coeffRef(2) += -mv1.coeff(0)*mv2.coeff(9) - mv1.coeff(1)*mv2.coeff(11) - mv1.coeff(2)*mv2.coeff(13) + mv1.coeff(3)*mv2.coeff(3) - mv1.coeff(4)*mv2.coeff(2) + mv1.coeff(5)*mv2.coeff(5) + mv1.coeff(6)*mv2.coeff(7) - mv1.coeff(10)*mv2.coeff(0) - mv1.coeff(12)*mv2.coeff(1);
		mv3.coeffRef(3) += -mv1.coeff(1)*mv2.coeff(12) - mv1.coeff(2)*mv2.coeff(14) - mv1.coeff(3)*mv2.coeff(15) + mv1.coeff(5)*mv2.coeff(6) + mv1.coeff(6)*mv2.coeff(8) + mv1.coeff(7)*mv2.coeff(9) - mv1.coeff(11)*mv2.coeff(0) - mv1.coeff(13)*mv2.coeff(1) - mv1.coeff(14)*mv2.coeff(2);
		mv3.coeffRef(4) +=  mv1.coeff(0)*mv2.coeff(10) - mv1.coeff(1)*mv2.coeff(8) + mv1.coeff(2)*mv2.coeff(6) + mv1.coeff(3)*mv2.coeff(16) - mv1.coeff(4)*mv2.coeff(4) - mv1.coeff(5)*mv2.coeff(1) + mv1.coeff(6)*mv2.coeff(0) - mv1.coeff(10)*mv2.coeff(7) + mv1.coeff(12)*mv2.coeff(5);
		mv3.coeffRef(5) +=  mv1.coeff(0)*mv2.coeff(11) - mv1.coeff(1)*mv2.coeff(9) - mv1.coeff(2)*mv2.coeff(16) + mv1.coeff(3)*mv2.coeff(6) - mv1.coeff(4)*mv2.coeff(5) - mv1.coeff(5)*mv2.coeff(2) + mv1.coeff(7)*mv2.coeff(0) + mv1.coeff(9)*mv2.coeff(7) - mv1.coeff(12)*mv2.coeff(4);
		mv3.coeffRef(6) +=  mv1.coeff(0)*mv2.coeff(12) - mv1.coeff(2)*mv2.coeff(17) - mv1.coeff(3)*mv2.coeff(18) - mv1.coeff(5)*mv2.coeff(3) + mv1.coeff(8)*mv2.coeff(0) + mv1.coeff(9)*mv2.coeff(8) + mv1.coeff(10)*mv2.coeff(9) - mv1.coeff(13)*mv2.coeff(4) - mv1.coeff(14)*mv2.coeff(5);
		mv3.coeffRef(7) +=  mv1.coeff(0)*mv2.coeff(13) + mv1.coeff(1)*mv2.coeff(16) - mv1.coeff(2)*mv2.coeff(9) + mv1.coeff(3)*mv2.coeff(8) - mv1.coeff(4)*mv2.coeff(7) - mv1.coeff(6)*mv2.coeff(2) + mv1.coeff(7)*mv2.coeff(1) - mv1.coeff(9)*mv2.coeff(5) + mv1.coeff(10)*mv2.coeff(4);
		mv3.coeffRef(8) +=  mv1.coeff(0)*mv2.coeff(14) + mv1.coeff(1)*mv2.coeff(17) - mv1.coeff(3)*mv2.coeff(19) - mv1.coeff(6)*mv2.coeff(3) + mv1.coeff(8)*mv2.coeff(1) - mv1.coeff(9)*mv2.coeff(6) + mv1.coeff(11)*mv2.coeff(4) + mv1.coeff(12)*mv2.coeff(9) - mv1.coeff(14)*mv2.coeff(7);
		mv3.coeffRef(9) +=  mv1.coeff(0)*mv2.coeff(15) + mv1.coeff(1)*mv2.coeff(18) + mv1.coeff(2)*mv2.coeff(19) - mv1.coeff(7)*mv2.coeff(3) + mv1.coeff(8)*mv2.coeff(2) - mv1.coeff(10)*mv2.coeff(6) + mv1.coeff(11)*mv2.coeff(5) - mv1.coeff(12)*mv2.coeff(8) + mv1.coeff(13)*mv2.coeff(7);
		mv3.coeffRef(10) +=  mv1.coeff(0)*mv2.coeff(17) - mv1.coeff(1)*mv2.coeff(14) + mv1.coeff(2)*mv2.coeff(12) + mv1.coeff(7)*mv2.coeff(16) - mv1.coeff(8)*mv2.coeff(4) - mv1.coeff(10)*mv2.coeff(13) + mv1.coeff(11)*mv2.coeff(1) + mv1.coeff(12)*mv2.coeff(11) - mv1.coeff(13)*mv2.coeff(0);
		mv3.coeffRef(11) +=  mv1.coeff(0)*mv2.coeff(18) - mv1.coeff(1)*mv2.coeff(15) + mv1.coeff(3)*mv2.coeff(12) - mv1.coeff(6)*mv2.coeff(16) - mv1.coeff(8)*mv2.coeff(5) + mv1.coeff(9)*mv2.coeff(13) + mv1.coeff(11)*mv2.coeff(2) - mv1.coeff(12)*mv2.coeff(10) - mv1.coeff(14)*mv2.coeff(0);
		mv3.coeffRef(12) +=  mv1.coeff(4)*mv2.coeff(12) - mv1.coeff(6)*mv2.coeff(17) - mv1.coeff(7)*mv2.coeff(18) - mv1.coeff(8)*mv2.coeff(6) + mv1.coeff(9)*mv2.coeff(14) + mv1.coeff(10)*mv2.coeff(15) + mv1.coeff(11)*mv2.coeff(3) - mv1.coeff(13)*mv2.coeff(10) - mv1.coeff(14)*mv2.coeff(11);
		mv3.coeffRef(13) +=  mv1.coeff(0)*mv2.coeff(19) - mv1.coeff(2)*mv2.coeff(15) + mv1.coeff(3)*mv2.coeff(14) + mv1.coeff(5)*mv2.coeff(16) - mv1.coeff(8)*mv2.coeff(7) - mv1.coeff(9)*mv2.coeff(11) + mv1.coeff(10)*mv2.coeff(10) + mv1.coeff(13)*mv2.coeff(2) - mv1.coeff(14)*mv2.coeff(1);
		mv3.coeffRef(14) +=  mv1.coeff(4)*mv2.coeff(14) + mv1.coeff(5)*mv2.coeff(17) - mv1.coeff(7)*mv2.coeff(19) - mv1.coeff(8)*mv2.coeff(8) - mv1.coeff(9)*mv2.coeff(12) + mv1.coeff(11)*mv2.coeff(10) + mv1.coeff(12)*mv2.coeff(15) + mv1.coeff(13)*mv2.coeff(3) - mv1.coeff(14)*mv2.coeff(13);
		mv3.coeffRef(15) +=  mv1.coeff(4)*mv2.coeff(15) + mv1.coeff(5)*mv2.coeff(18) + mv1.coeff(6)*mv2.coeff(19) - mv1.coeff(8)*mv2.coeff(9) - mv1.coeff(10)*mv2.coeff(12) + mv1.coeff(11)*mv2.coeff(11) - mv1.coeff(12)*mv2.coeff(14) + mv1.coeff(13)*mv2.coeff(13) + mv1.coeff(14)*mv2.coeff(3);
		mv3.coeffRef(16) +=  mv1.coeff(1)*mv2.coeff(19) - mv1.coeff(2)*mv2.coeff(18) + mv1.coeff(3)*mv2.coeff(17) - mv1.coeff(5)*mv2.coeff(13) + mv1.coeff(6)*mv2.coeff(11) - mv1.coeff(7)*mv2.coeff(10) - mv1.coeff(11)*mv2.coeff(7) + mv1.coeff(13)*mv2.coeff(5) - mv1.coeff(14)*mv2.coeff(4);
		mv3.coeffRef(17) +=  mv1.coeff(4)*mv2.coeff(17) - mv1.coeff(5)*mv2.coeff(14) + mv1.coeff(6)*mv2.coeff(12) - mv1.coeff(8)*mv2.coeff(10) - mv1.coeff(10)*mv2.coeff(19) - mv1.coeff(11)*mv2.coeff(8) + mv1.coeff(12)*mv2.coeff(18) + mv1.coeff(13)*mv2.coeff(6) - mv1.coeff(14)*mv2.coeff(16);
		mv3.coeffRef(18) +=  mv1.coeff(4)*mv2.coeff(18) - mv1.coeff(5)*mv2.coeff(15) + mv1.coeff(7)*mv2.coeff(12) - mv1.coeff(8)*mv2.coeff(11) + mv1.coeff(9)*mv2.coeff(19) - mv1.coeff(11)*mv2.coeff(9) - mv1.coeff(12)*mv2.coeff(17) + mv1.coeff(13)*mv2.coeff(16) + mv1.coeff(14)*mv2.coeff(6);
		mv3.coeffRef(19) +=  mv1.coeff(4)*mv2.coeff(19) - mv1.coeff(6)*mv2.coeff(15) + mv1.coeff(7)*mv2.coeff(14) - mv1.coeff(8)*mv2.coeff(13) - mv1.coeff(9)*mv2.coeff(18) + mv1.coeff(10)*mv2.coeff(17) - mv1.coeff(11)*mv2.coeff(16) - mv1.coeff(13)*mv2.coeff(9) + mv1.coeff(14)*mv2.coeff(8);
	}


	/// \brief Compute the geometric product between two homogeneous multivectors mv1 (grade 2) and mv2 (grade 4). 
	/// \tparam the type of value that we manipulate, either float or double or something else.
	/// \param mv1 - the first homogeneous multivector of grade 2 represented as an Eigen::VectorXd
	/// \param mv2 - the second homogeneous multivector of grade 4 represented as a Eigen::VectorXd
	/// \param mv3 - the result of mv1 mv2 whose grade is 4
	template<typename T>
	void geometric_2_4_4(const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv1, const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv2, Eigen::Matrix<T, Eigen::Dynamic, 1>& mv3){
		mv3.coeffRef(0) +=  mv1.coeff(0)*mv2.coeff(7) - mv1.coeff(1)*mv2.coeff(4) + mv1.coeff(2)*mv2.coeff(2) - mv1.coeff(3)*mv2.coeff(10) - mv1.coeff(4)*mv2.coeff(0) + mv1.coeff(7)*mv2.coeff(6) - mv1.coeff(10)*mv2.coeff(3) + mv1.coeff(12)*mv2.coeff(1);
		mv3.coeffRef(1) +=  mv1.coeff(0)*mv2.coeff(8) - mv1.coeff(1)*mv2.coeff(5) + mv1.coeff(2)*mv2.coeff(10) + mv1.coeff(3)*mv2.coeff(2) - mv1.coeff(4)*mv2.coeff(1) - mv1.coeff(6)*mv2.coeff(6) + mv1.coeff(9)*mv2.coeff(3) - mv1.coeff(12)*mv2.coeff(0);
		mv3.coeffRef(2) +=  mv1.coeff(2)*mv2.coeff(11) + mv1.coeff(3)*mv2.coeff(12) - mv1.coeff(6)*mv2.coeff(7) - mv1.coeff(7)*mv2.coeff(8) + mv1.coeff(9)*mv2.coeff(4) + mv1.coeff(10)*mv2.coeff(5) - mv1.coeff(13)*mv2.coeff(0) - mv1.coeff(14)*mv2.coeff(1);
		mv3.coeffRef(3) +=  mv1.coeff(0)*mv2.coeff(9) - mv1.coeff(1)*mv2.coeff(10) - mv1.coeff(2)*mv2.coeff(5) + mv1.coeff(3)*mv2.coeff(4) - mv1.coeff(4)*mv2.coeff(3) + mv1.coeff(5)*mv2.coeff(6) - mv1.coeff(9)*mv2.coeff(1) + mv1.coeff(10)*mv2.coeff(0);
		mv3.coeffRef(4) += -mv1.coeff(1)*mv2.coeff(11) + mv1.coeff(3)*mv2.coeff(13) + mv1.coeff(5)*mv2.coeff(7) - mv1.coeff(7)*mv2.coeff(9) - mv1.coeff(9)*mv2.coeff(2) + mv1.coeff(11)*mv2.coeff(0) + mv1.coeff(12)*mv2.coeff(5) - mv1.coeff(14)*mv2.coeff(3);
		mv3.coeffRef(5) += -mv1.coeff(1)*mv2.coeff(12) - mv1.coeff(2)*mv2.coeff(13) + mv1.coeff(5)*mv2.coeff(8) + mv1.coeff(6)*mv2.coeff(9) - mv1.coeff(10)*mv2.coeff(2) + mv1.coeff(11)*mv2.coeff(1) - mv1.coeff(12)*mv2.coeff(4) + mv1.coeff(13)*mv2.coeff(3);
		mv3.coeffRef(6) +=  mv1.coeff(0)*mv2.coeff(10) + mv1.coeff(1)*mv2.coeff(9) - mv1.coeff(2)*mv2.coeff(8) + mv1.coeff(3)*mv2.coeff(7) - mv1.coeff(4)*mv2.coeff(6) - mv1.coeff(5)*mv2.coeff(3) + mv1.coeff(6)*mv2.coeff(1) - mv1.coeff(7)*mv2.coeff(0);
		mv3.coeffRef(7) +=  mv1.coeff(0)*mv2.coeff(11) + mv1.coeff(3)*mv2.coeff(14) - mv1.coeff(5)*mv2.coeff(4) + mv1.coeff(6)*mv2.coeff(2) - mv1.coeff(8)*mv2.coeff(0) - mv1.coeff(10)*mv2.coeff(9) + mv1.coeff(12)*mv2.coeff(8) - mv1.coeff(14)*mv2.coeff(6);
		mv3.coeffRef(8) +=  mv1.coeff(0)*mv2.coeff(12) - mv1.coeff(2)*mv2.coeff(14) - mv1.coeff(5)*mv2.coeff(5) + mv1.coeff(7)*mv2.coeff(2) - mv1.coeff(8)*mv2.coeff(1) + mv1.coeff(9)*mv2.coeff(9) - mv1.coeff(12)*mv2.coeff(7) + mv1.coeff(13)*mv2.coeff(6);
		mv3.coeffRef(9) +=  mv1.coeff(0)*mv2.coeff(13) + mv1.coeff(1)*mv2.coeff(14) - mv1.coeff(6)*mv2.coeff(5) + mv1.coeff(7)*mv2.coeff(4) - mv1.coeff(8)*mv2.coeff(3) - mv1.coeff(9)*mv2.coeff(8) + mv1.coeff(10)*mv2.coeff(7) - mv1.coeff(11)*mv2.coeff(6);
		mv3.coeffRef(10) += -mv1.coeff(0)*mv2.coeff(14) + mv1.coeff(1)*mv2.coeff(13) - mv1.coeff(2)*mv2.coeff(12) + mv1.coeff(3)*mv2.coeff(11) - mv1.coeff(8)*mv2.coeff(6) + mv1.coeff(11)*mv2.coeff(3) - mv1.coeff(13)*mv2.coeff(1) + mv1.coeff(14)*mv2.coeff(0);
		mv3.coeffRef(11) +=  mv1.coeff(4)*mv2.coeff(11) + mv1.coeff(7)*mv2.coeff(14) - mv1.coeff(8)*mv2.coeff(7) - mv1.coeff(10)*mv2.coeff(13) + mv1.coeff(11)*mv2.coeff(4) + mv1.coeff(12)*mv2.coeff(12) - mv1.coeff(13)*mv2.coeff(2) - mv1.coeff(14)*mv2.coeff(10);
		mv3.coeffRef(12) +=  mv1.coeff(4)*mv2.coeff(12) - mv1.coeff(6)*mv2.coeff(14) - mv1.coeff(8)*mv2.coeff(8) + mv1.coeff(9)*mv2.coeff(13) + mv1.coeff(11)*mv2.coeff(5) - mv1.coeff(12)*mv2.coeff(11) + mv1.coeff(13)*mv2.coeff(10) - mv1.coeff(14)*mv2.coeff(2);
		mv3.coeffRef(13) +=  mv1.coeff(4)*mv2.coeff(13) + mv1.coeff(5)*mv2.coeff(14) - mv1.coeff(8)*mv2.coeff(9) - mv1.coeff(9)*mv2.coeff(12) + mv1.coeff(10)*mv2.coeff(11) - mv1.coeff(11)*mv2.coeff(10) + mv1.coeff(13)*mv2.coeff(5) - mv1.coeff(14)*mv2.coeff(4);
		mv3.coeffRef(14) +=  mv1.coeff(4)*mv2.coeff(14) - mv1.coeff(5)*mv2.coeff(13) + mv1.coeff(6)*mv2.coeff(12) - mv1.coeff(7)*mv2.coeff(11) + mv1.coeff(8)*mv2.coeff(10) - mv1.coeff(11)*mv2.coeff(9) + mv1.coeff(13)*mv2.coeff(8) - mv1.coeff(14)*mv2.coeff(7);
	}


	/// \brief Compute the geometric product between two homogeneous multivectors mv1 (grade 2) and mv2 (grade 5). 
	/// \tparam the type of value that we manipulate, either float or double or something else.
	/// \param mv1 - the first homogeneous multivector of grade 2 represented as an Eigen::VectorXd
	/// \param mv2 - the second homogeneous multivector of grade 5 represented as a Eigen::VectorXd
	/// \param mv3 - the result of mv1 mv2 whose grade is 5
	template<typename T>
	void geometric_2_5_5(const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv1, const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv2, Eigen::Matrix<T, Eigen::Dynamic, 1>& mv3){
		mv3.coeffRef(0) += -mv1.coeff(0)*mv2.coeff(4) + mv1.coeff(1)*mv2.coeff(3) - mv1.coeff(2)*mv2.coeff(2) + mv1.coeff(3)*mv2.coeff(1) - mv1.coeff(4)*mv2.coeff(0);
		mv3.coeffRef(1) += -mv1.coeff(3)*mv2.coeff(5) + mv1.coeff(7)*mv2.coeff(4) - mv1.coeff(10)*mv2.coeff(3) + mv1.coeff(12)*mv2.coeff(2) - mv1.coeff(14)*mv2.coeff(0);
		mv3.coeffRef(2) +=  mv1.coeff(2)*mv2.coeff(5) - mv1.coeff(6)*mv2.coeff(4) + mv1.coeff(9)*mv2.coeff(3) - mv1.coeff(12)*mv2.coeff(1) + mv1.coeff(13)*mv2.coeff(0);
		mv3.coeffRef(3) += -mv1.coeff(1)*mv2.coeff(5) + mv1.coeff(5)*mv2.coeff(4) - mv1.coeff(9)*mv2.coeff(2) + mv1.coeff(10)*mv2.coeff(1) - mv1.coeff(11)*mv2.coeff(0);
		mv3.coeffRef(4) +=  mv1.coeff(0)*mv2.coeff(5) - mv1.coeff(5)*mv2.coeff(3) + mv1.coeff(6)*mv2.coeff(2) - mv1.coeff(7)*mv2.coeff(1) + mv1.coeff(8)*mv2.coeff(0);
		mv3.coeffRef(5) +=  mv1.coeff(4)*mv2.coeff(5) - mv1.coeff(8)*mv2.coeff(4) + mv1.coeff(11)*mv2.coeff(3) - mv1.coeff(13)*mv2.coeff(2) + mv1.coeff(14)*mv2.coeff(1);
	}


	/// \brief Compute the geometric product between two homogeneous multivectors mv1 (grade 3) and mv2 (grade 2). 
	/// \tparam the type of value that we manipulate, either float or double or something else.
	/// \param mv1 - the first homogeneous multivector of grade 3 represented as an Eigen::VectorXd
	/// \param mv2 - the second homogeneous multivector of grade 2 represented as a Eigen::VectorXd
	/// \param mv3 - the result of mv1 mv2 whose grade is 3
	template<typename T>
	void geometric_3_2_3(const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv1, const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv2, Eigen::Matrix<T, Eigen::Dynamic, 1>& mv3){
		mv3.coeffRef(0) +=  mv1.coeff(0)*mv2.coeff(4) - mv1.coeff(1)*mv2.coeff(9) - mv1.coeff(2)*mv2.coeff(10) - mv1.coeff(3)*mv2.coeff(1) + mv1.coeff(4)*mv2.coeff(6) + mv1.coeff(5)*mv2.coeff(7) + mv1.coeff(6)*mv2.coeff(0) - mv1.coeff(10)*mv2.coeff(2) - mv1.coeff(11)*mv2.coeff(3);
		mv3.coeffRef(1) +=  mv1.coeff(0)*mv2.coeff(9) + mv1.coeff(1)*mv2.coeff(4) - mv1.coeff(2)*mv2.coeff(12) - mv1.coeff(3)*mv2.coeff(2) - mv1.coeff(4)*mv2.coeff(5) + mv1.coeff(7)*mv2.coeff(7) + mv1.coeff(8)*mv2.coeff(0) + mv1.coeff(10)*mv2.coeff(1) - mv1.coeff(13)*mv2.coeff(3);
		mv3.coeffRef(2) +=  mv1.coeff(0)*mv2.coeff(10) + mv1.coeff(1)*mv2.coeff(12) + mv1.coeff(2)*mv2.coeff(4) - mv1.coeff(3)*mv2.coeff(3) - mv1.coeff(5)*mv2.coeff(5) - mv1.coeff(7)*mv2.coeff(6) + mv1.coeff(9)*mv2.coeff(0) + mv1.coeff(11)*mv2.coeff(1) + mv1.coeff(13)*mv2.coeff(2);
		mv3.coeffRef(3) +=  mv1.coeff(0)*mv2.coeff(11) + mv1.coeff(1)*mv2.coeff(13) + mv1.coeff(2)*mv2.coeff(14) - mv1.coeff(6)*mv2.coeff(5) - mv1.coeff(8)*mv2.coeff(6) - mv1.coeff(9)*mv2.coeff(7) + mv1.coeff(12)*mv2.coeff(1) + mv1.coeff(14)*mv2.coeff(2) + mv1.coeff(15)*mv2.coeff(3);
		mv3.coeffRef(4) += -mv1.coeff(0)*mv2.coeff(6) + mv1.coeff(1)*mv2.coeff(5) + mv1.coeff(4)*mv2.coeff(4) - mv1.coeff(5)*mv2.coeff(12) - mv1.coeff(6)*mv2.coeff(2) + mv1.coeff(7)*mv2.coeff(10) + mv1.coeff(8)*mv2.coeff(1) - mv1.coeff(10)*mv2.coeff(0) - mv1.coeff(16)*mv2.coeff(3);
		mv3.coeffRef(5) += -mv1.coeff(0)*mv2.coeff(7) + mv1.coeff(2)*mv2.coeff(5) + mv1.coeff(4)*mv2.coeff(12) + mv1.coeff(5)*mv2.coeff(4) - mv1.coeff(6)*mv2.coeff(3) - mv1.coeff(7)*mv2.coeff(9) + mv1.coeff(9)*mv2.coeff(1) - mv1.coeff(11)*mv2.coeff(0) + mv1.coeff(16)*mv2.coeff(2);
		mv3.coeffRef(6) += -mv1.coeff(0)*mv2.coeff(8) + mv1.coeff(3)*mv2.coeff(5) + mv1.coeff(4)*mv2.coeff(13) + mv1.coeff(5)*mv2.coeff(14) - mv1.coeff(8)*mv2.coeff(9) - mv1.coeff(9)*mv2.coeff(10) - mv1.coeff(12)*mv2.coeff(0) + mv1.coeff(17)*mv2.coeff(2) + mv1.coeff(18)*mv2.coeff(3);
		mv3.coeffRef(7) += -mv1.coeff(1)*mv2.coeff(7) + mv1.coeff(2)*mv2.coeff(6) - mv1.coeff(4)*mv2.coeff(10) + mv1.coeff(5)*mv2.coeff(9) + mv1.coeff(7)*mv2.coeff(4) - mv1.coeff(8)*mv2.coeff(3) + mv1.coeff(9)*mv2.coeff(2) - mv1.coeff(13)*mv2.coeff(0) - mv1.coeff(16)*mv2.coeff(1);
		mv3.coeffRef(8) += -mv1.coeff(1)*mv2.coeff(8) + mv1.coeff(3)*mv2.coeff(6) - mv1.coeff(4)*mv2.coeff(11) + mv1.coeff(6)*mv2.coeff(9) + mv1.coeff(7)*mv2.coeff(14) - mv1.coeff(9)*mv2.coeff(12) - mv1.coeff(14)*mv2.coeff(0) - mv1.coeff(17)*mv2.coeff(1) + mv1.coeff(19)*mv2.coeff(3);
		mv3.coeffRef(9) += -mv1.coeff(2)*mv2.coeff(8) + mv1.coeff(3)*mv2.coeff(7) - mv1.coeff(5)*mv2.coeff(11) + mv1.coeff(6)*mv2.coeff(10) - mv1.coeff(7)*mv2.coeff(13) + mv1.coeff(8)*mv2.coeff(12) - mv1.coeff(15)*mv2.coeff(0) - mv1.coeff(18)*mv2.coeff(1) - mv1.coeff(19)*mv2.coeff(2);
		mv3.coeffRef(10) +=  mv1.coeff(0)*mv2.coeff(13) - mv1.coeff(1)*mv2.coeff(11) + mv1.coeff(4)*mv2.coeff(8) - mv1.coeff(11)*mv2.coeff(12) - mv1.coeff(12)*mv2.coeff(2) + mv1.coeff(13)*mv2.coeff(10) + mv1.coeff(14)*mv2.coeff(1) - mv1.coeff(16)*mv2.coeff(7) - mv1.coeff(17)*mv2.coeff(0);
		mv3.coeffRef(11) +=  mv1.coeff(0)*mv2.coeff(14) - mv1.coeff(2)*mv2.coeff(11) + mv1.coeff(5)*mv2.coeff(8) + mv1.coeff(10)*mv2.coeff(12) - mv1.coeff(12)*mv2.coeff(3) - mv1.coeff(13)*mv2.coeff(9) + mv1.coeff(15)*mv2.coeff(1) + mv1.coeff(16)*mv2.coeff(6) - mv1.coeff(18)*mv2.coeff(0);
		mv3.coeffRef(12) += -mv1.coeff(3)*mv2.coeff(11) + mv1.coeff(6)*mv2.coeff(8) + mv1.coeff(10)*mv2.coeff(13) + mv1.coeff(11)*mv2.coeff(14) - mv1.coeff(12)*mv2.coeff(4) - mv1.coeff(14)*mv2.coeff(9) - mv1.coeff(15)*mv2.coeff(10) + mv1.coeff(17)*mv2.coeff(6) + mv1.coeff(18)*mv2.coeff(7);
		mv3.coeffRef(13) +=  mv1.coeff(1)*mv2.coeff(14) - mv1.coeff(2)*mv2.coeff(13) + mv1.coeff(7)*mv2.coeff(8) - mv1.coeff(10)*mv2.coeff(10) + mv1.coeff(11)*mv2.coeff(9) - mv1.coeff(14)*mv2.coeff(3) + mv1.coeff(15)*mv2.coeff(2) - mv1.coeff(16)*mv2.coeff(5) - mv1.coeff(19)*mv2.coeff(0);
		mv3.coeffRef(14) += -mv1.coeff(3)*mv2.coeff(13) + mv1.coeff(8)*mv2.coeff(8) - mv1.coeff(10)*mv2.coeff(11) + mv1.coeff(12)*mv2.coeff(9) + mv1.coeff(13)*mv2.coeff(14) - mv1.coeff(14)*mv2.coeff(4) - mv1.coeff(15)*mv2.coeff(12) - mv1.coeff(17)*mv2.coeff(5) + mv1.coeff(19)*mv2.coeff(7);
		mv3.coeffRef(15) += -mv1.coeff(3)*mv2.coeff(14) + mv1.coeff(9)*mv2.coeff(8) - mv1.coeff(11)*mv2.coeff(11) + mv1.coeff(12)*mv2.coeff(10) - mv1.coeff(13)*mv2.coeff(13) + mv1.coeff(14)*mv2.coeff(12) - mv1.coeff(15)*mv2.coeff(4) - mv1.coeff(18)*mv2.coeff(5) - mv1.coeff(19)*mv2.coeff(6);
		mv3.coeffRef(16) +=  mv1.coeff(4)*mv2.coeff(14) - mv1.coeff(5)*mv2.coeff(13) + mv1.coeff(7)*mv2.coeff(11) + mv1.coeff(10)*mv2.coeff(7) - mv1.coeff(11)*mv2.coeff(6) + mv1.coeff(13)*mv2.coeff(5) - mv1.coeff(17)*mv2.coeff(3) + mv1.coeff(18)*mv2.coeff(2) - mv1.coeff(19)*mv2.coeff(1);
		mv3.coeffRef(17) += -mv1.coeff(6)*mv2.coeff(13) + mv1.coeff(8)*mv2.coeff(11) + mv1.coeff(10)*mv2.coeff(8) - mv1.coeff(12)*mv2.coeff(6) + mv1.coeff(14)*mv2.coeff(5) + mv1.coeff(16)*mv2.coeff(14) - mv1.coeff(17)*mv2.coeff(4) - mv1.coeff(18)*mv2.coeff(12) + mv1.coeff(19)*mv2.coeff(10);
		mv3.coeffRef(18) += -mv1.coeff(6)*mv2.coeff(14) + mv1.coeff(9)*mv2.coeff(11) + mv1.coeff(11)*mv2.coeff(8) - mv1.coeff(12)*mv2.coeff(7) + mv1.coeff(15)*mv2.coeff(5) - mv1.coeff(16)*mv2.coeff(13) + mv1.coeff(17)*mv2.coeff(12) - mv1.coeff(18)*mv2.coeff(4) - mv1.coeff(19)*mv2.coeff(9);
		mv3.coeffRef(19) += -mv1.coeff(8)*mv2.coeff(14) + mv1.coeff(9)*mv2.coeff(13) + mv1.coeff(13)*mv2.coeff(8) - mv1.coeff(14)*mv2.coeff(7) + mv1.coeff(15)*mv2.coeff(6) + mv1.coeff(16)*mv2.coeff(11) - mv1.coeff(17)*mv2.coeff(10) + mv1.coeff(18)*mv2.coeff(9) - mv1.coeff(19)*mv2.coeff(4);
	}


	/// \brief Compute the geometric product between two homogeneous multivectors mv1 (grade 3) and mv2 (grade 3). 
	/// \tparam the type of value that we manipulate, either float or double or something else.
	/// \param mv1 - the first homogeneous multivector of grade 3 represented as an Eigen::VectorXd
	/// \param mv2 - the second homogeneous multivector of grade 3 represented as a Eigen::VectorXd
	/// \param mv3 - the result of mv1 mv2 whose grade is 2
	template<typename T>
	void geometric_3_3_2(const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv1, const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv2, Eigen::Matrix<T, Eigen::Dynamic, 1>& mv3){
		mv3.coeffRef(0) += -mv1.coeff(0)*mv2.coeff(6) - mv1.coeff(1)*mv2.coeff(8) - mv1.coeff(2)*mv2.coeff(9) - mv1.coeff(4)*mv2.coeff(10) - mv1.coeff(5)*mv2.coeff(11) + mv1.coeff(6)*mv2.coeff(0) - mv1.coeff(7)*mv2.coeff(13) + mv1.coeff(8)*mv2.coeff(1) + mv1.coeff(9)*mv2.coeff(2) + mv1.coeff(10)*mv2.coeff(4) + mv1.coeff(11)*mv2.coeff(5) + mv1.coeff(13)*mv2.coeff(7);
		mv3.coeffRef(1) +=  mv1.coeff(0)*mv2.coeff(3) + mv1.coeff(1)*mv2.coeff(10) + mv1.coeff(2)*mv2.coeff(11) - mv1.coeff(3)*mv2.coeff(0) - mv1.coeff(4)*mv2.coeff(8) - mv1.coeff(5)*mv2.coeff(9) - mv1.coeff(7)*mv2.coeff(16) + mv1.coeff(8)*mv2.coeff(4) + mv1.coeff(9)*mv2.coeff(5) - mv1.coeff(10)*mv2.coeff(1) - mv1.coeff(11)*mv2.coeff(2) + mv1.coeff(16)*mv2.coeff(7);
		mv3.coeffRef(2) += -mv1.coeff(0)*mv2.coeff(10) + mv1.coeff(1)*mv2.coeff(3) + mv1.coeff(2)*mv2.coeff(13) - mv1.coeff(3)*mv2.coeff(1) + mv1.coeff(4)*mv2.coeff(6) + mv1.coeff(5)*mv2.coeff(16) - mv1.coeff(6)*mv2.coeff(4) - mv1.coeff(7)*mv2.coeff(9) + mv1.coeff(9)*mv2.coeff(7) + mv1.coeff(10)*mv2.coeff(0) - mv1.coeff(13)*mv2.coeff(2) - mv1.coeff(16)*mv2.coeff(5);
		mv3.coeffRef(3) += -mv1.coeff(0)*mv2.coeff(11) - mv1.coeff(1)*mv2.coeff(13) + mv1.coeff(2)*mv2.coeff(3) - mv1.coeff(3)*mv2.coeff(2) - mv1.coeff(4)*mv2.coeff(16) + mv1.coeff(5)*mv2.coeff(6) - mv1.coeff(6)*mv2.coeff(5) + mv1.coeff(7)*mv2.coeff(8) - mv1.coeff(8)*mv2.coeff(7) + mv1.coeff(11)*mv2.coeff(0) + mv1.coeff(13)*mv2.coeff(1) + mv1.coeff(16)*mv2.coeff(4);
		mv3.coeffRef(4) += -mv1.coeff(0)*mv2.coeff(12) - mv1.coeff(1)*mv2.coeff(14) - mv1.coeff(2)*mv2.coeff(15) - mv1.coeff(4)*mv2.coeff(17) - mv1.coeff(5)*mv2.coeff(18) - mv1.coeff(7)*mv2.coeff(19) + mv1.coeff(12)*mv2.coeff(0) + mv1.coeff(14)*mv2.coeff(1) + mv1.coeff(15)*mv2.coeff(2) + mv1.coeff(17)*mv2.coeff(4) + mv1.coeff(18)*mv2.coeff(5) + mv1.coeff(19)*mv2.coeff(7);
		mv3.coeffRef(5) +=  mv1.coeff(1)*mv2.coeff(17) + mv1.coeff(2)*mv2.coeff(18) + mv1.coeff(3)*mv2.coeff(6) - mv1.coeff(4)*mv2.coeff(14) - mv1.coeff(5)*mv2.coeff(15) - mv1.coeff(6)*mv2.coeff(3) - mv1.coeff(13)*mv2.coeff(16) + mv1.coeff(14)*mv2.coeff(4) + mv1.coeff(15)*mv2.coeff(5) + mv1.coeff(16)*mv2.coeff(13) - mv1.coeff(17)*mv2.coeff(1) - mv1.coeff(18)*mv2.coeff(2);
		mv3.coeffRef(6) += -mv1.coeff(0)*mv2.coeff(17) + mv1.coeff(2)*mv2.coeff(19) + mv1.coeff(3)*mv2.coeff(8) + mv1.coeff(4)*mv2.coeff(12) - mv1.coeff(7)*mv2.coeff(15) - mv1.coeff(8)*mv2.coeff(3) + mv1.coeff(11)*mv2.coeff(16) - mv1.coeff(12)*mv2.coeff(4) + mv1.coeff(15)*mv2.coeff(7) - mv1.coeff(16)*mv2.coeff(11) + mv1.coeff(17)*mv2.coeff(0) - mv1.coeff(19)*mv2.coeff(2);
		mv3.coeffRef(7) += -mv1.coeff(0)*mv2.coeff(18) - mv1.coeff(1)*mv2.coeff(19) + mv1.coeff(3)*mv2.coeff(9) + mv1.coeff(5)*mv2.coeff(12) + mv1.coeff(7)*mv2.coeff(14) - mv1.coeff(9)*mv2.coeff(3) - mv1.coeff(10)*mv2.coeff(16) - mv1.coeff(12)*mv2.coeff(5) - mv1.coeff(14)*mv2.coeff(7) + mv1.coeff(16)*mv2.coeff(10) + mv1.coeff(18)*mv2.coeff(0) + mv1.coeff(19)*mv2.coeff(1);
		mv3.coeffRef(8) +=  mv1.coeff(6)*mv2.coeff(12) + mv1.coeff(8)*mv2.coeff(14) + mv1.coeff(9)*mv2.coeff(15) - mv1.coeff(10)*mv2.coeff(17) - mv1.coeff(11)*mv2.coeff(18) - mv1.coeff(12)*mv2.coeff(6) - mv1.coeff(13)*mv2.coeff(19) - mv1.coeff(14)*mv2.coeff(8) - mv1.coeff(15)*mv2.coeff(9) + mv1.coeff(17)*mv2.coeff(10) + mv1.coeff(18)*mv2.coeff(11) + mv1.coeff(19)*mv2.coeff(13);
		mv3.coeffRef(9) +=  mv1.coeff(0)*mv2.coeff(14) - mv1.coeff(1)*mv2.coeff(12) + mv1.coeff(5)*mv2.coeff(19) + mv1.coeff(6)*mv2.coeff(8) - mv1.coeff(7)*mv2.coeff(18) - mv1.coeff(8)*mv2.coeff(6) - mv1.coeff(11)*mv2.coeff(13) + mv1.coeff(12)*mv2.coeff(1) + mv1.coeff(13)*mv2.coeff(11) - mv1.coeff(14)*mv2.coeff(0) + mv1.coeff(18)*mv2.coeff(7) - mv1.coeff(19)*mv2.coeff(5);
		mv3.coeffRef(10) +=  mv1.coeff(0)*mv2.coeff(15) - mv1.coeff(2)*mv2.coeff(12) - mv1.coeff(4)*mv2.coeff(19) + mv1.coeff(6)*mv2.coeff(9) + mv1.coeff(7)*mv2.coeff(17) - mv1.coeff(9)*mv2.coeff(6) + mv1.coeff(10)*mv2.coeff(13) + mv1.coeff(12)*mv2.coeff(2) - mv1.coeff(13)*mv2.coeff(10) - mv1.coeff(15)*mv2.coeff(0) - mv1.coeff(17)*mv2.coeff(7) + mv1.coeff(19)*mv2.coeff(4);
		mv3.coeffRef(11) += -mv1.coeff(3)*mv2.coeff(12) + mv1.coeff(8)*mv2.coeff(17) + mv1.coeff(9)*mv2.coeff(18) + mv1.coeff(10)*mv2.coeff(14) + mv1.coeff(11)*mv2.coeff(15) + mv1.coeff(12)*mv2.coeff(3) - mv1.coeff(14)*mv2.coeff(10) - mv1.coeff(15)*mv2.coeff(11) - mv1.coeff(16)*mv2.coeff(19) - mv1.coeff(17)*mv2.coeff(8) - mv1.coeff(18)*mv2.coeff(9) + mv1.coeff(19)*mv2.coeff(16);
		mv3.coeffRef(12) +=  mv1.coeff(1)*mv2.coeff(15) - mv1.coeff(2)*mv2.coeff(14) + mv1.coeff(4)*mv2.coeff(18) - mv1.coeff(5)*mv2.coeff(17) + mv1.coeff(8)*mv2.coeff(9) - mv1.coeff(9)*mv2.coeff(8) - mv1.coeff(10)*mv2.coeff(11) + mv1.coeff(11)*mv2.coeff(10) + mv1.coeff(14)*mv2.coeff(2) - mv1.coeff(15)*mv2.coeff(1) + mv1.coeff(17)*mv2.coeff(5) - mv1.coeff(18)*mv2.coeff(4);
		mv3.coeffRef(13) += -mv1.coeff(3)*mv2.coeff(14) - mv1.coeff(6)*mv2.coeff(17) + mv1.coeff(9)*mv2.coeff(19) - mv1.coeff(10)*mv2.coeff(12) + mv1.coeff(12)*mv2.coeff(10) + mv1.coeff(13)*mv2.coeff(15) + mv1.coeff(14)*mv2.coeff(3) - mv1.coeff(15)*mv2.coeff(13) + mv1.coeff(16)*mv2.coeff(18) + mv1.coeff(17)*mv2.coeff(6) - mv1.coeff(18)*mv2.coeff(16) - mv1.coeff(19)*mv2.coeff(9);
		mv3.coeffRef(14) += -mv1.coeff(3)*mv2.coeff(15) - mv1.coeff(6)*mv2.coeff(18) - mv1.coeff(8)*mv2.coeff(19) - mv1.coeff(11)*mv2.coeff(12) + mv1.coeff(12)*mv2.coeff(11) - mv1.coeff(13)*mv2.coeff(14) + mv1.coeff(14)*mv2.coeff(13) + mv1.coeff(15)*mv2.coeff(3) - mv1.coeff(16)*mv2.coeff(17) + mv1.coeff(17)*mv2.coeff(16) + mv1.coeff(18)*mv2.coeff(6) + mv1.coeff(19)*mv2.coeff(8);
	}


	/// \brief Compute the geometric product between two homogeneous multivectors mv1 (grade 3) and mv2 (grade 3). 
	/// \tparam the type of value that we manipulate, either float or double or something else.
	/// \param mv1 - the first homogeneous multivector of grade 3 represented as an Eigen::VectorXd
	/// \param mv2 - the second homogeneous multivector of grade 3 represented as a Eigen::VectorXd
	/// \param mv3 - the result of mv1 mv2 whose grade is 4
	template<typename T>
	void geometric_3_3_4(const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv1, const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv2, Eigen::Matrix<T, Eigen::Dynamic, 1>& mv3){
		mv3.coeffRef(0) += -mv1.coeff(0)*mv2.coeff(8) + mv1.coeff(1)*mv2.coeff(6) + mv1.coeff(2)*mv2.coeff(16) - mv1.coeff(3)*mv2.coeff(4) - mv1.coeff(4)*mv2.coeff(3) - mv1.coeff(5)*mv2.coeff(13) + mv1.coeff(6)*mv2.coeff(1) + mv1.coeff(7)*mv2.coeff(11) - mv1.coeff(8)*mv2.coeff(0) + mv1.coeff(11)*mv2.coeff(7) - mv1.coeff(13)*mv2.coeff(5) + mv1.coeff(16)*mv2.coeff(2);
		mv3.coeffRef(1) += -mv1.coeff(0)*mv2.coeff(9) - mv1.coeff(1)*mv2.coeff(16) + mv1.coeff(2)*mv2.coeff(6) - mv1.coeff(3)*mv2.coeff(5) + mv1.coeff(4)*mv2.coeff(13) - mv1.coeff(5)*mv2.coeff(3) + mv1.coeff(6)*mv2.coeff(2) - mv1.coeff(7)*mv2.coeff(10) - mv1.coeff(9)*mv2.coeff(0) - mv1.coeff(10)*mv2.coeff(7) + mv1.coeff(13)*mv2.coeff(4) - mv1.coeff(16)*mv2.coeff(1);
		mv3.coeffRef(2) += -mv1.coeff(1)*mv2.coeff(17) - mv1.coeff(2)*mv2.coeff(18) + mv1.coeff(4)*mv2.coeff(14) + mv1.coeff(5)*mv2.coeff(15) - mv1.coeff(8)*mv2.coeff(10) - mv1.coeff(9)*mv2.coeff(11) - mv1.coeff(10)*mv2.coeff(8) - mv1.coeff(11)*mv2.coeff(9) + mv1.coeff(14)*mv2.coeff(4) + mv1.coeff(15)*mv2.coeff(5) - mv1.coeff(17)*mv2.coeff(1) - mv1.coeff(18)*mv2.coeff(2);
		mv3.coeffRef(3) +=  mv1.coeff(0)*mv2.coeff(16) - mv1.coeff(1)*mv2.coeff(9) + mv1.coeff(2)*mv2.coeff(8) - mv1.coeff(3)*mv2.coeff(7) - mv1.coeff(4)*mv2.coeff(11) + mv1.coeff(5)*mv2.coeff(10) - mv1.coeff(7)*mv2.coeff(3) + mv1.coeff(8)*mv2.coeff(2) - mv1.coeff(9)*mv2.coeff(1) + mv1.coeff(10)*mv2.coeff(5) - mv1.coeff(11)*mv2.coeff(4) + mv1.coeff(16)*mv2.coeff(0);
		mv3.coeffRef(4) +=  mv1.coeff(0)*mv2.coeff(17) - mv1.coeff(2)*mv2.coeff(19) - mv1.coeff(4)*mv2.coeff(12) + mv1.coeff(6)*mv2.coeff(10) + mv1.coeff(7)*mv2.coeff(15) - mv1.coeff(9)*mv2.coeff(13) + mv1.coeff(10)*mv2.coeff(6) - mv1.coeff(12)*mv2.coeff(4) - mv1.coeff(13)*mv2.coeff(9) + mv1.coeff(15)*mv2.coeff(7) + mv1.coeff(17)*mv2.coeff(0) - mv1.coeff(19)*mv2.coeff(2);
		mv3.coeffRef(5) +=  mv1.coeff(0)*mv2.coeff(18) + mv1.coeff(1)*mv2.coeff(19) - mv1.coeff(5)*mv2.coeff(12) + mv1.coeff(6)*mv2.coeff(11) - mv1.coeff(7)*mv2.coeff(14) + mv1.coeff(8)*mv2.coeff(13) + mv1.coeff(11)*mv2.coeff(6) - mv1.coeff(12)*mv2.coeff(5) + mv1.coeff(13)*mv2.coeff(8) - mv1.coeff(14)*mv2.coeff(7) + mv1.coeff(18)*mv2.coeff(0) + mv1.coeff(19)*mv2.coeff(1);
		mv3.coeffRef(6) += -mv1.coeff(0)*mv2.coeff(13) + mv1.coeff(1)*mv2.coeff(11) - mv1.coeff(2)*mv2.coeff(10) - mv1.coeff(4)*mv2.coeff(9) + mv1.coeff(5)*mv2.coeff(8) - mv1.coeff(6)*mv2.coeff(7) - mv1.coeff(7)*mv2.coeff(6) + mv1.coeff(8)*mv2.coeff(5) - mv1.coeff(9)*mv2.coeff(4) - mv1.coeff(10)*mv2.coeff(2) + mv1.coeff(11)*mv2.coeff(1) - mv1.coeff(13)*mv2.coeff(0);
		mv3.coeffRef(7) += -mv1.coeff(0)*mv2.coeff(14) + mv1.coeff(1)*mv2.coeff(12) - mv1.coeff(3)*mv2.coeff(10) - mv1.coeff(5)*mv2.coeff(19) + mv1.coeff(7)*mv2.coeff(18) - mv1.coeff(9)*mv2.coeff(16) - mv1.coeff(10)*mv2.coeff(3) + mv1.coeff(12)*mv2.coeff(1) - mv1.coeff(14)*mv2.coeff(0) - mv1.coeff(16)*mv2.coeff(9) + mv1.coeff(18)*mv2.coeff(7) - mv1.coeff(19)*mv2.coeff(5);
		mv3.coeffRef(8) += -mv1.coeff(0)*mv2.coeff(15) + mv1.coeff(2)*mv2.coeff(12) - mv1.coeff(3)*mv2.coeff(11) + mv1.coeff(4)*mv2.coeff(19) - mv1.coeff(7)*mv2.coeff(17) + mv1.coeff(8)*mv2.coeff(16) - mv1.coeff(11)*mv2.coeff(3) + mv1.coeff(12)*mv2.coeff(2) - mv1.coeff(15)*mv2.coeff(0) + mv1.coeff(16)*mv2.coeff(8) - mv1.coeff(17)*mv2.coeff(7) + mv1.coeff(19)*mv2.coeff(4);
		mv3.coeffRef(9) += -mv1.coeff(1)*mv2.coeff(15) + mv1.coeff(2)*mv2.coeff(14) - mv1.coeff(3)*mv2.coeff(13) - mv1.coeff(4)*mv2.coeff(18) + mv1.coeff(5)*mv2.coeff(17) - mv1.coeff(6)*mv2.coeff(16) - mv1.coeff(13)*mv2.coeff(3) + mv1.coeff(14)*mv2.coeff(2) - mv1.coeff(15)*mv2.coeff(1) - mv1.coeff(16)*mv2.coeff(6) + mv1.coeff(17)*mv2.coeff(5) - mv1.coeff(18)*mv2.coeff(4);
		mv3.coeffRef(10) += -mv1.coeff(0)*mv2.coeff(19) + mv1.coeff(1)*mv2.coeff(18) - mv1.coeff(2)*mv2.coeff(17) - mv1.coeff(4)*mv2.coeff(15) + mv1.coeff(5)*mv2.coeff(14) - mv1.coeff(7)*mv2.coeff(12) - mv1.coeff(12)*mv2.coeff(7) + mv1.coeff(14)*mv2.coeff(5) - mv1.coeff(15)*mv2.coeff(4) - mv1.coeff(17)*mv2.coeff(2) + mv1.coeff(18)*mv2.coeff(1) - mv1.coeff(19)*mv2.coeff(0);
		mv3.coeffRef(11) += -mv1.coeff(3)*mv2.coeff(17) + mv1.coeff(6)*mv2.coeff(14) - mv1.coeff(8)*mv2.coeff(12) - mv1.coeff(11)*mv2.coeff(19) - mv1.coeff(12)*mv2.coeff(8) + mv1.coeff(13)*mv2.coeff(18) + mv1.coeff(14)*mv2.coeff(6) - mv1.coeff(15)*mv2.coeff(16) - mv1.coeff(16)*mv2.coeff(15) - mv1.coeff(17)*mv2.coeff(3) + mv1.coeff(18)*mv2.coeff(13) - mv1.coeff(19)*mv2.coeff(11);
		mv3.coeffRef(12) += -mv1.coeff(3)*mv2.coeff(18) + mv1.coeff(6)*mv2.coeff(15) - mv1.coeff(9)*mv2.coeff(12) + mv1.coeff(10)*mv2.coeff(19) - mv1.coeff(12)*mv2.coeff(9) - mv1.coeff(13)*mv2.coeff(17) + mv1.coeff(14)*mv2.coeff(16) + mv1.coeff(15)*mv2.coeff(6) + mv1.coeff(16)*mv2.coeff(14) - mv1.coeff(17)*mv2.coeff(13) - mv1.coeff(18)*mv2.coeff(3) + mv1.coeff(19)*mv2.coeff(10);
		mv3.coeffRef(13) += -mv1.coeff(3)*mv2.coeff(19) + mv1.coeff(8)*mv2.coeff(15) - mv1.coeff(9)*mv2.coeff(14) - mv1.coeff(10)*mv2.coeff(18) + mv1.coeff(11)*mv2.coeff(17) - mv1.coeff(12)*mv2.coeff(16) - mv1.coeff(14)*mv2.coeff(9) + mv1.coeff(15)*mv2.coeff(8) - mv1.coeff(16)*mv2.coeff(12) + mv1.coeff(17)*mv2.coeff(11) - mv1.coeff(18)*mv2.coeff(10) - mv1.coeff(19)*mv2.coeff(3);
		mv3.coeffRef(14) += -mv1.coeff(6)*mv2.coeff(19) + mv1.coeff(8)*mv2.coeff(18) - mv1.coeff(9)*mv2.coeff(17) + mv1.coeff(10)*mv2.coeff(15) - mv1.coeff(11)*mv2.coeff(14) + mv1.coeff(12)*mv2.coeff(13) + mv1.coeff(13)*mv2.coeff(12) - mv1.coeff(14)*mv2.coeff(11) + mv1.coeff(15)*mv2.coeff(10) - mv1.coeff(17)*mv2.coeff(9) + mv1.coeff(18)*mv2.coeff(8) - mv1.coeff(19)*mv2.coeff(6);
	}


	/// \brief Compute the geometric product between two homogeneous multivectors mv1 (grade 3) and mv2 (grade 4). 
	/// \tparam the type of value that we manipulate, either float or double or something else.
	/// \param mv1 - the first homogeneous multivector of grade 3 represented as an Eigen::VectorXd
	/// \param mv2 - the second homogeneous multivector of grade 4 represented as a Eigen::VectorXd
	/// \param mv3 - the result of mv1 mv2 whose grade is 3
	template<typename T>
	void geometric_3_4_3(const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv1, const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv2, Eigen::Matrix<T, Eigen::Dynamic, 1>& mv3){
		mv3.coeffRef(0) += -mv1.coeff(1)*mv2.coeff(7) - mv1.coeff(2)*mv2.coeff(8) + mv1.coeff(4)*mv2.coeff(4) + mv1.coeff(5)*mv2.coeff(5) - mv1.coeff(7)*mv2.coeff(10) - mv1.coeff(8)*mv2.coeff(0) - mv1.coeff(9)*mv2.coeff(1) + mv1.coeff(13)*mv2.coeff(6) - mv1.coeff(16)*mv2.coeff(3);
		mv3.coeffRef(1) +=  mv1.coeff(0)*mv2.coeff(7) - mv1.coeff(2)*mv2.coeff(9) - mv1.coeff(4)*mv2.coeff(2) + mv1.coeff(5)*mv2.coeff(10) + mv1.coeff(6)*mv2.coeff(0) + mv1.coeff(7)*mv2.coeff(5) - mv1.coeff(9)*mv2.coeff(3) - mv1.coeff(11)*mv2.coeff(6) + mv1.coeff(16)*mv2.coeff(1);
		mv3.coeffRef(2) +=  mv1.coeff(0)*mv2.coeff(8) + mv1.coeff(1)*mv2.coeff(9) - mv1.coeff(4)*mv2.coeff(10) - mv1.coeff(5)*mv2.coeff(2) + mv1.coeff(6)*mv2.coeff(1) - mv1.coeff(7)*mv2.coeff(4) + mv1.coeff(8)*mv2.coeff(3) + mv1.coeff(10)*mv2.coeff(6) - mv1.coeff(16)*mv2.coeff(0);
		mv3.coeffRef(3) += -mv1.coeff(4)*mv2.coeff(11) - mv1.coeff(5)*mv2.coeff(12) - mv1.coeff(7)*mv2.coeff(13) + mv1.coeff(10)*mv2.coeff(7) + mv1.coeff(11)*mv2.coeff(8) + mv1.coeff(13)*mv2.coeff(9) - mv1.coeff(17)*mv2.coeff(0) - mv1.coeff(18)*mv2.coeff(1) - mv1.coeff(19)*mv2.coeff(3);
		mv3.coeffRef(4) += -mv1.coeff(0)*mv2.coeff(4) + mv1.coeff(1)*mv2.coeff(2) - mv1.coeff(2)*mv2.coeff(10) - mv1.coeff(3)*mv2.coeff(0) - mv1.coeff(5)*mv2.coeff(9) + mv1.coeff(7)*mv2.coeff(8) - mv1.coeff(9)*mv2.coeff(6) + mv1.coeff(11)*mv2.coeff(3) - mv1.coeff(13)*mv2.coeff(1);
		mv3.coeffRef(5) += -mv1.coeff(0)*mv2.coeff(5) + mv1.coeff(1)*mv2.coeff(10) + mv1.coeff(2)*mv2.coeff(2) - mv1.coeff(3)*mv2.coeff(1) + mv1.coeff(4)*mv2.coeff(9) - mv1.coeff(7)*mv2.coeff(7) + mv1.coeff(8)*mv2.coeff(6) - mv1.coeff(10)*mv2.coeff(3) + mv1.coeff(13)*mv2.coeff(0);
		mv3.coeffRef(6) +=  mv1.coeff(1)*mv2.coeff(11) + mv1.coeff(2)*mv2.coeff(12) - mv1.coeff(7)*mv2.coeff(14) - mv1.coeff(10)*mv2.coeff(4) - mv1.coeff(11)*mv2.coeff(5) + mv1.coeff(14)*mv2.coeff(0) + mv1.coeff(15)*mv2.coeff(1) + mv1.coeff(16)*mv2.coeff(9) - mv1.coeff(19)*mv2.coeff(6);
		mv3.coeffRef(7) += -mv1.coeff(0)*mv2.coeff(10) - mv1.coeff(1)*mv2.coeff(5) + mv1.coeff(2)*mv2.coeff(4) - mv1.coeff(3)*mv2.coeff(3) - mv1.coeff(4)*mv2.coeff(8) + mv1.coeff(5)*mv2.coeff(7) - mv1.coeff(6)*mv2.coeff(6) + mv1.coeff(10)*mv2.coeff(1) - mv1.coeff(11)*mv2.coeff(0);
		mv3.coeffRef(8) += -mv1.coeff(0)*mv2.coeff(11) + mv1.coeff(2)*mv2.coeff(13) + mv1.coeff(5)*mv2.coeff(14) + mv1.coeff(10)*mv2.coeff(2) - mv1.coeff(12)*mv2.coeff(0) - mv1.coeff(13)*mv2.coeff(5) + mv1.coeff(15)*mv2.coeff(3) - mv1.coeff(16)*mv2.coeff(8) + mv1.coeff(18)*mv2.coeff(6);
		mv3.coeffRef(9) += -mv1.coeff(0)*mv2.coeff(12) - mv1.coeff(1)*mv2.coeff(13) - mv1.coeff(4)*mv2.coeff(14) + mv1.coeff(11)*mv2.coeff(2) - mv1.coeff(12)*mv2.coeff(1) + mv1.coeff(13)*mv2.coeff(4) - mv1.coeff(14)*mv2.coeff(3) + mv1.coeff(16)*mv2.coeff(7) - mv1.coeff(17)*mv2.coeff(6);
		mv3.coeffRef(10) +=  mv1.coeff(2)*mv2.coeff(14) - mv1.coeff(3)*mv2.coeff(7) - mv1.coeff(5)*mv2.coeff(13) + mv1.coeff(6)*mv2.coeff(4) + mv1.coeff(7)*mv2.coeff(12) - mv1.coeff(8)*mv2.coeff(2) - mv1.coeff(15)*mv2.coeff(6) + mv1.coeff(18)*mv2.coeff(3) - mv1.coeff(19)*mv2.coeff(1);
		mv3.coeffRef(11) += -mv1.coeff(1)*mv2.coeff(14) - mv1.coeff(3)*mv2.coeff(8) + mv1.coeff(4)*mv2.coeff(13) + mv1.coeff(6)*mv2.coeff(5) - mv1.coeff(7)*mv2.coeff(11) - mv1.coeff(9)*mv2.coeff(2) + mv1.coeff(14)*mv2.coeff(6) - mv1.coeff(17)*mv2.coeff(3) + mv1.coeff(19)*mv2.coeff(0);
		mv3.coeffRef(12) += -mv1.coeff(8)*mv2.coeff(11) - mv1.coeff(9)*mv2.coeff(12) - mv1.coeff(13)*mv2.coeff(14) + mv1.coeff(14)*mv2.coeff(7) + mv1.coeff(15)*mv2.coeff(8) + mv1.coeff(16)*mv2.coeff(13) - mv1.coeff(17)*mv2.coeff(4) - mv1.coeff(18)*mv2.coeff(5) - mv1.coeff(19)*mv2.coeff(10);
		mv3.coeffRef(13) +=  mv1.coeff(0)*mv2.coeff(14) - mv1.coeff(3)*mv2.coeff(9) - mv1.coeff(4)*mv2.coeff(12) + mv1.coeff(5)*mv2.coeff(11) + mv1.coeff(8)*mv2.coeff(5) - mv1.coeff(9)*mv2.coeff(4) - mv1.coeff(12)*mv2.coeff(6) + mv1.coeff(17)*mv2.coeff(1) - mv1.coeff(18)*mv2.coeff(0);
		mv3.coeffRef(14) +=  mv1.coeff(6)*mv2.coeff(11) - mv1.coeff(9)*mv2.coeff(13) + mv1.coeff(11)*mv2.coeff(14) - mv1.coeff(12)*mv2.coeff(7) + mv1.coeff(15)*mv2.coeff(9) - mv1.coeff(16)*mv2.coeff(12) + mv1.coeff(17)*mv2.coeff(2) + mv1.coeff(18)*mv2.coeff(10) - mv1.coeff(19)*mv2.coeff(5);
		mv3.coeffRef(15) +=  mv1.coeff(6)*mv2.coeff(12) + mv1.coeff(8)*mv2.coeff(13) - mv1.coeff(10)*mv2.coeff(14) - mv1.coeff(12)*mv2.coeff(8) - mv1.coeff(14)*mv2.coeff(9) + mv1.coeff(16)*mv2.coeff(11) - mv1.coeff(17)*mv2.coeff(10) + mv1.coeff(18)*mv2.coeff(2) + mv1.coeff(19)*mv2.coeff(4);
		mv3.coeffRef(16) += -mv1.coeff(0)*mv2.coeff(13) + mv1.coeff(1)*mv2.coeff(12) - mv1.coeff(2)*mv2.coeff(11) - mv1.coeff(6)*mv2.coeff(9) + mv1.coeff(8)*mv2.coeff(8) - mv1.coeff(9)*mv2.coeff(7) + mv1.coeff(12)*mv2.coeff(3) - mv1.coeff(14)*mv2.coeff(1) + mv1.coeff(15)*mv2.coeff(0);
		mv3.coeffRef(17) += -mv1.coeff(3)*mv2.coeff(11) - mv1.coeff(9)*mv2.coeff(14) - mv1.coeff(11)*mv2.coeff(13) + mv1.coeff(12)*mv2.coeff(4) + mv1.coeff(13)*mv2.coeff(12) - mv1.coeff(14)*mv2.coeff(2) - mv1.coeff(15)*mv2.coeff(10) + mv1.coeff(18)*mv2.coeff(9) - mv1.coeff(19)*mv2.coeff(8);
		mv3.coeffRef(18) += -mv1.coeff(3)*mv2.coeff(12) + mv1.coeff(8)*mv2.coeff(14) + mv1.coeff(10)*mv2.coeff(13) + mv1.coeff(12)*mv2.coeff(5) - mv1.coeff(13)*mv2.coeff(11) + mv1.coeff(14)*mv2.coeff(10) - mv1.coeff(15)*mv2.coeff(2) - mv1.coeff(17)*mv2.coeff(9) + mv1.coeff(19)*mv2.coeff(7);
		mv3.coeffRef(19) += -mv1.coeff(3)*mv2.coeff(13) - mv1.coeff(6)*mv2.coeff(14) - mv1.coeff(10)*mv2.coeff(12) + mv1.coeff(11)*mv2.coeff(11) - mv1.coeff(12)*mv2.coeff(10) + mv1.coeff(14)*mv2.coeff(5) - mv1.coeff(15)*mv2.coeff(4) + mv1.coeff(17)*mv2.coeff(8) - mv1.coeff(18)*mv2.coeff(7);
	}


	/// \brief Compute the geometric product between two homogeneous multivectors mv1 (grade 3) and mv2 (grade 4). 
	/// \tparam the type of value that we manipulate, either float or double or something else.
	/// \param mv1 - the first homogeneous multivector of grade 3 represented as an Eigen::VectorXd
	/// \param mv2 - the second homogeneous multivector of grade 4 represented as a Eigen::VectorXd
	/// \param mv3 - the result of mv1 mv2 whose grade is 5
	template<typename T>
	void geometric_3_4_5(const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv1, const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv2, Eigen::Matrix<T, Eigen::Dynamic, 1>& mv3){
		mv3.coeffRef(0) +=  mv1.coeff(0)*mv2.coeff(9) - mv1.coeff(1)*mv2.coeff(8) + mv1.coeff(2)*mv2.coeff(7) - mv1.coeff(3)*mv2.coeff(6) + mv1.coeff(4)*mv2.coeff(5) - mv1.coeff(5)*mv2.coeff(4) + mv1.coeff(6)*mv2.coeff(3) + mv1.coeff(7)*mv2.coeff(2) - mv1.coeff(8)*mv2.coeff(1) + mv1.coeff(9)*mv2.coeff(0);
		mv3.coeffRef(1) +=  mv1.coeff(2)*mv2.coeff(14) - mv1.coeff(5)*mv2.coeff(13) + mv1.coeff(7)*mv2.coeff(12) - mv1.coeff(9)*mv2.coeff(10) + mv1.coeff(11)*mv2.coeff(9) - mv1.coeff(13)*mv2.coeff(8) + mv1.coeff(15)*mv2.coeff(6) + mv1.coeff(16)*mv2.coeff(5) - mv1.coeff(18)*mv2.coeff(3) + mv1.coeff(19)*mv2.coeff(1);
		mv3.coeffRef(2) += -mv1.coeff(1)*mv2.coeff(14) + mv1.coeff(4)*mv2.coeff(13) - mv1.coeff(7)*mv2.coeff(11) + mv1.coeff(8)*mv2.coeff(10) - mv1.coeff(10)*mv2.coeff(9) + mv1.coeff(13)*mv2.coeff(7) - mv1.coeff(14)*mv2.coeff(6) - mv1.coeff(16)*mv2.coeff(4) + mv1.coeff(17)*mv2.coeff(3) - mv1.coeff(19)*mv2.coeff(0);
		mv3.coeffRef(3) +=  mv1.coeff(0)*mv2.coeff(14) - mv1.coeff(4)*mv2.coeff(12) + mv1.coeff(5)*mv2.coeff(11) - mv1.coeff(6)*mv2.coeff(10) + mv1.coeff(10)*mv2.coeff(8) - mv1.coeff(11)*mv2.coeff(7) + mv1.coeff(12)*mv2.coeff(6) + mv1.coeff(16)*mv2.coeff(2) - mv1.coeff(17)*mv2.coeff(1) + mv1.coeff(18)*mv2.coeff(0);
		mv3.coeffRef(4) += -mv1.coeff(0)*mv2.coeff(13) + mv1.coeff(1)*mv2.coeff(12) - mv1.coeff(2)*mv2.coeff(11) + mv1.coeff(3)*mv2.coeff(10) - mv1.coeff(10)*mv2.coeff(5) + mv1.coeff(11)*mv2.coeff(4) - mv1.coeff(12)*mv2.coeff(3) - mv1.coeff(13)*mv2.coeff(2) + mv1.coeff(14)*mv2.coeff(1) - mv1.coeff(15)*mv2.coeff(0);
		mv3.coeffRef(5) += -mv1.coeff(3)*mv2.coeff(14) + mv1.coeff(6)*mv2.coeff(13) - mv1.coeff(8)*mv2.coeff(12) + mv1.coeff(9)*mv2.coeff(11) - mv1.coeff(12)*mv2.coeff(9) + mv1.coeff(14)*mv2.coeff(8) - mv1.coeff(15)*mv2.coeff(7) - mv1.coeff(17)*mv2.coeff(5) + mv1.coeff(18)*mv2.coeff(4) - mv1.coeff(19)*mv2.coeff(2);
	}


	/// \brief Compute the geometric product between two homogeneous multivectors mv1 (grade 3) and mv2 (grade 5). 
	/// \tparam the type of value that we manipulate, either float or double or something else.
	/// \param mv1 - the first homogeneous multivector of grade 3 represented as an Eigen::VectorXd
	/// \param mv2 - the second homogeneous multivector of grade 5 represented as a Eigen::VectorXd
	/// \param mv3 - the result of mv1 mv2 whose grade is 4
	template<typename T>
	void geometric_3_5_4(const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv1, const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv2, Eigen::Matrix<T, Eigen::Dynamic, 1>& mv3){
		mv3.coeffRef(0) += -mv1.coeff(2)*mv2.coeff(4) + mv1.coeff(5)*mv2.coeff(3) - mv1.coeff(7)*mv2.coeff(2) + mv1.coeff(9)*mv2.coeff(0);
		mv3.coeffRef(1) +=  mv1.coeff(1)*mv2.coeff(4) - mv1.coeff(4)*mv2.coeff(3) + mv1.coeff(7)*mv2.coeff(1) - mv1.coeff(8)*mv2.coeff(0);
		mv3.coeffRef(2) += -mv1.coeff(7)*mv2.coeff(5) + mv1.coeff(13)*mv2.coeff(4) - mv1.coeff(16)*mv2.coeff(3) + mv1.coeff(19)*mv2.coeff(0);
		mv3.coeffRef(3) += -mv1.coeff(0)*mv2.coeff(4) + mv1.coeff(4)*mv2.coeff(2) - mv1.coeff(5)*mv2.coeff(1) + mv1.coeff(6)*mv2.coeff(0);
		mv3.coeffRef(4) +=  mv1.coeff(5)*mv2.coeff(5) - mv1.coeff(11)*mv2.coeff(4) + mv1.coeff(16)*mv2.coeff(2) - mv1.coeff(18)*mv2.coeff(0);
		mv3.coeffRef(5) += -mv1.coeff(4)*mv2.coeff(5) + mv1.coeff(10)*mv2.coeff(4) - mv1.coeff(16)*mv2.coeff(1) + mv1.coeff(17)*mv2.coeff(0);
		mv3.coeffRef(6) +=  mv1.coeff(0)*mv2.coeff(3) - mv1.coeff(1)*mv2.coeff(2) + mv1.coeff(2)*mv2.coeff(1) - mv1.coeff(3)*mv2.coeff(0);
		mv3.coeffRef(7) += -mv1.coeff(2)*mv2.coeff(5) + mv1.coeff(11)*mv2.coeff(3) - mv1.coeff(13)*mv2.coeff(2) + mv1.coeff(15)*mv2.coeff(0);
		mv3.coeffRef(8) +=  mv1.coeff(1)*mv2.coeff(5) - mv1.coeff(10)*mv2.coeff(3) + mv1.coeff(13)*mv2.coeff(1) - mv1.coeff(14)*mv2.coeff(0);
		mv3.coeffRef(9) += -mv1.coeff(0)*mv2.coeff(5) + mv1.coeff(10)*mv2.coeff(2) - mv1.coeff(11)*mv2.coeff(1) + mv1.coeff(12)*mv2.coeff(0);
		mv3.coeffRef(10) +=  mv1.coeff(3)*mv2.coeff(4) - mv1.coeff(6)*mv2.coeff(3) + mv1.coeff(8)*mv2.coeff(2) - mv1.coeff(9)*mv2.coeff(1);
		mv3.coeffRef(11) +=  mv1.coeff(9)*mv2.coeff(5) - mv1.coeff(15)*mv2.coeff(4) + mv1.coeff(18)*mv2.coeff(3) - mv1.coeff(19)*mv2.coeff(2);
		mv3.coeffRef(12) += -mv1.coeff(8)*mv2.coeff(5) + mv1.coeff(14)*mv2.coeff(4) - mv1.coeff(17)*mv2.coeff(3) + mv1.coeff(19)*mv2.coeff(1);
		mv3.coeffRef(13) +=  mv1.coeff(6)*mv2.coeff(5) - mv1.coeff(12)*mv2.coeff(4) + mv1.coeff(17)*mv2.coeff(2) - mv1.coeff(18)*mv2.coeff(1);
		mv3.coeffRef(14) += -mv1.coeff(3)*mv2.coeff(5) + mv1.coeff(12)*mv2.coeff(3) - mv1.coeff(14)*mv2.coeff(2) + mv1.coeff(15)*mv2.coeff(1);
	}


	/// \brief Compute the geometric product between two homogeneous multivectors mv1 (grade 4) and mv2 (grade 2). 
	/// \tparam the type of value that we manipulate, either float or double or something else.
	/// \param mv1 - the first homogeneous multivector of grade 4 represented as an Eigen::VectorXd
	/// \param mv2 - the second homogeneous multivector of grade 2 represented as a Eigen::VectorXd
	/// \param mv3 - the result of mv1 mv2 whose grade is 4
	template<typename T>
	void geometric_4_2_4(const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv1, const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv2, Eigen::Matrix<T, Eigen::Dynamic, 1>& mv3){
		mv3.coeffRef(0) +=  mv1.coeff(0)*mv2.coeff(4) - mv1.coeff(1)*mv2.coeff(12) - mv1.coeff(2)*mv2.coeff(2) + mv1.coeff(3)*mv2.coeff(10) + mv1.coeff(4)*mv2.coeff(1) - mv1.coeff(6)*mv2.coeff(7) - mv1.coeff(7)*mv2.coeff(0) + mv1.coeff(10)*mv2.coeff(3);
		mv3.coeffRef(1) +=  mv1.coeff(0)*mv2.coeff(12) + mv1.coeff(1)*mv2.coeff(4) - mv1.coeff(2)*mv2.coeff(3) - mv1.coeff(3)*mv2.coeff(9) + mv1.coeff(5)*mv2.coeff(1) + mv1.coeff(6)*mv2.coeff(6) - mv1.coeff(8)*mv2.coeff(0) - mv1.coeff(10)*mv2.coeff(2);
		mv3.coeffRef(2) +=  mv1.coeff(0)*mv2.coeff(13) + mv1.coeff(1)*mv2.coeff(14) - mv1.coeff(4)*mv2.coeff(9) - mv1.coeff(5)*mv2.coeff(10) + mv1.coeff(7)*mv2.coeff(6) + mv1.coeff(8)*mv2.coeff(7) - mv1.coeff(11)*mv2.coeff(2) - mv1.coeff(12)*mv2.coeff(3);
		mv3.coeffRef(3) += -mv1.coeff(0)*mv2.coeff(10) + mv1.coeff(1)*mv2.coeff(9) + mv1.coeff(3)*mv2.coeff(4) - mv1.coeff(4)*mv2.coeff(3) + mv1.coeff(5)*mv2.coeff(2) - mv1.coeff(6)*mv2.coeff(5) - mv1.coeff(9)*mv2.coeff(0) + mv1.coeff(10)*mv2.coeff(1);
		mv3.coeffRef(4) += -mv1.coeff(0)*mv2.coeff(11) + mv1.coeff(2)*mv2.coeff(9) + mv1.coeff(3)*mv2.coeff(14) - mv1.coeff(5)*mv2.coeff(12) - mv1.coeff(7)*mv2.coeff(5) + mv1.coeff(9)*mv2.coeff(7) + mv1.coeff(11)*mv2.coeff(1) - mv1.coeff(13)*mv2.coeff(3);
		mv3.coeffRef(5) += -mv1.coeff(1)*mv2.coeff(11) + mv1.coeff(2)*mv2.coeff(10) - mv1.coeff(3)*mv2.coeff(13) + mv1.coeff(4)*mv2.coeff(12) - mv1.coeff(8)*mv2.coeff(5) - mv1.coeff(9)*mv2.coeff(6) + mv1.coeff(12)*mv2.coeff(1) + mv1.coeff(13)*mv2.coeff(2);
		mv3.coeffRef(6) +=  mv1.coeff(0)*mv2.coeff(7) - mv1.coeff(1)*mv2.coeff(6) + mv1.coeff(3)*mv2.coeff(5) + mv1.coeff(6)*mv2.coeff(4) - mv1.coeff(7)*mv2.coeff(3) + mv1.coeff(8)*mv2.coeff(2) - mv1.coeff(9)*mv2.coeff(1) - mv1.coeff(10)*mv2.coeff(0);
		mv3.coeffRef(7) +=  mv1.coeff(0)*mv2.coeff(8) - mv1.coeff(2)*mv2.coeff(6) + mv1.coeff(4)*mv2.coeff(5) + mv1.coeff(6)*mv2.coeff(14) - mv1.coeff(8)*mv2.coeff(12) + mv1.coeff(9)*mv2.coeff(10) - mv1.coeff(11)*mv2.coeff(0) - mv1.coeff(14)*mv2.coeff(3);
		mv3.coeffRef(8) +=  mv1.coeff(1)*mv2.coeff(8) - mv1.coeff(2)*mv2.coeff(7) + mv1.coeff(5)*mv2.coeff(5) - mv1.coeff(6)*mv2.coeff(13) + mv1.coeff(7)*mv2.coeff(12) - mv1.coeff(9)*mv2.coeff(9) - mv1.coeff(12)*mv2.coeff(0) + mv1.coeff(14)*mv2.coeff(2);
		mv3.coeffRef(9) +=  mv1.coeff(3)*mv2.coeff(8) - mv1.coeff(4)*mv2.coeff(7) + mv1.coeff(5)*mv2.coeff(6) + mv1.coeff(6)*mv2.coeff(11) - mv1.coeff(7)*mv2.coeff(10) + mv1.coeff(8)*mv2.coeff(9) - mv1.coeff(13)*mv2.coeff(0) - mv1.coeff(14)*mv2.coeff(1);
		mv3.coeffRef(10) += -mv1.coeff(0)*mv2.coeff(14) + mv1.coeff(1)*mv2.coeff(13) - mv1.coeff(3)*mv2.coeff(11) + mv1.coeff(6)*mv2.coeff(8) - mv1.coeff(11)*mv2.coeff(3) + mv1.coeff(12)*mv2.coeff(2) - mv1.coeff(13)*mv2.coeff(1) + mv1.coeff(14)*mv2.coeff(0);
		mv3.coeffRef(11) +=  mv1.coeff(2)*mv2.coeff(13) - mv1.coeff(4)*mv2.coeff(11) + mv1.coeff(7)*mv2.coeff(8) + mv1.coeff(10)*mv2.coeff(14) - mv1.coeff(11)*mv2.coeff(4) - mv1.coeff(12)*mv2.coeff(12) + mv1.coeff(13)*mv2.coeff(10) - mv1.coeff(14)*mv2.coeff(7);
		mv3.coeffRef(12) +=  mv1.coeff(2)*mv2.coeff(14) - mv1.coeff(5)*mv2.coeff(11) + mv1.coeff(8)*mv2.coeff(8) - mv1.coeff(10)*mv2.coeff(13) + mv1.coeff(11)*mv2.coeff(12) - mv1.coeff(12)*mv2.coeff(4) - mv1.coeff(13)*mv2.coeff(9) + mv1.coeff(14)*mv2.coeff(6);
		mv3.coeffRef(13) +=  mv1.coeff(4)*mv2.coeff(14) - mv1.coeff(5)*mv2.coeff(13) + mv1.coeff(9)*mv2.coeff(8) + mv1.coeff(10)*mv2.coeff(11) - mv1.coeff(11)*mv2.coeff(10) + mv1.coeff(12)*mv2.coeff(9) - mv1.coeff(13)*mv2.coeff(4) - mv1.coeff(14)*mv2.coeff(5);
		mv3.coeffRef(14) +=  mv1.coeff(7)*mv2.coeff(14) - mv1.coeff(8)*mv2.coeff(13) + mv1.coeff(9)*mv2.coeff(11) - mv1.coeff(10)*mv2.coeff(8) + mv1.coeff(11)*mv2.coeff(7) - mv1.coeff(12)*mv2.coeff(6) + mv1.coeff(13)*mv2.coeff(5) - mv1.coeff(14)*mv2.coeff(4);
	}


	/// \brief Compute the geometric product between two homogeneous multivectors mv1 (grade 4) and mv2 (grade 3). 
	/// \tparam the type of value that we manipulate, either float or double or something else.
	/// \param mv1 - the first homogeneous multivector of grade 4 represented as an Eigen::VectorXd
	/// \param mv2 - the second homogeneous multivector of grade 3 represented as a Eigen::VectorXd
	/// \param mv3 - the result of mv1 mv2 whose grade is 3
	template<typename T>
	void geometric_4_3_3(const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv1, const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv2, Eigen::Matrix<T, Eigen::Dynamic, 1>& mv3){
		mv3.coeffRef(0) += -mv1.coeff(0)*mv2.coeff(8) - mv1.coeff(1)*mv2.coeff(9) - mv1.coeff(3)*mv2.coeff(16) + mv1.coeff(4)*mv2.coeff(4) + mv1.coeff(5)*mv2.coeff(5) + mv1.coeff(6)*mv2.coeff(13) - mv1.coeff(7)*mv2.coeff(1) - mv1.coeff(8)*mv2.coeff(2) - mv1.coeff(10)*mv2.coeff(7);
		mv3.coeffRef(1) +=  mv1.coeff(0)*mv2.coeff(6) + mv1.coeff(1)*mv2.coeff(16) - mv1.coeff(2)*mv2.coeff(4) - mv1.coeff(3)*mv2.coeff(9) + mv1.coeff(5)*mv2.coeff(7) - mv1.coeff(6)*mv2.coeff(11) + mv1.coeff(7)*mv2.coeff(0) - mv1.coeff(9)*mv2.coeff(2) + mv1.coeff(10)*mv2.coeff(5);
		mv3.coeffRef(2) += -mv1.coeff(0)*mv2.coeff(16) + mv1.coeff(1)*mv2.coeff(6) - mv1.coeff(2)*mv2.coeff(5) + mv1.coeff(3)*mv2.coeff(8) - mv1.coeff(4)*mv2.coeff(7) + mv1.coeff(6)*mv2.coeff(10) + mv1.coeff(8)*mv2.coeff(0) + mv1.coeff(9)*mv2.coeff(1) - mv1.coeff(10)*mv2.coeff(4);
		mv3.coeffRef(3) += -mv1.coeff(0)*mv2.coeff(17) - mv1.coeff(1)*mv2.coeff(18) - mv1.coeff(3)*mv2.coeff(19) + mv1.coeff(7)*mv2.coeff(10) + mv1.coeff(8)*mv2.coeff(11) + mv1.coeff(9)*mv2.coeff(13) - mv1.coeff(11)*mv2.coeff(4) - mv1.coeff(12)*mv2.coeff(5) - mv1.coeff(13)*mv2.coeff(7);
		mv3.coeffRef(4) += -mv1.coeff(0)*mv2.coeff(3) - mv1.coeff(1)*mv2.coeff(13) + mv1.coeff(2)*mv2.coeff(1) + mv1.coeff(3)*mv2.coeff(11) - mv1.coeff(4)*mv2.coeff(0) - mv1.coeff(6)*mv2.coeff(9) + mv1.coeff(8)*mv2.coeff(7) - mv1.coeff(9)*mv2.coeff(5) - mv1.coeff(10)*mv2.coeff(2);
		mv3.coeffRef(5) +=  mv1.coeff(0)*mv2.coeff(13) - mv1.coeff(1)*mv2.coeff(3) + mv1.coeff(2)*mv2.coeff(2) - mv1.coeff(3)*mv2.coeff(10) - mv1.coeff(5)*mv2.coeff(0) + mv1.coeff(6)*mv2.coeff(8) - mv1.coeff(7)*mv2.coeff(7) + mv1.coeff(9)*mv2.coeff(4) + mv1.coeff(10)*mv2.coeff(1);
		mv3.coeffRef(6) +=  mv1.coeff(0)*mv2.coeff(14) + mv1.coeff(1)*mv2.coeff(15) - mv1.coeff(4)*mv2.coeff(10) - mv1.coeff(5)*mv2.coeff(11) - mv1.coeff(6)*mv2.coeff(19) + mv1.coeff(9)*mv2.coeff(16) + mv1.coeff(11)*mv2.coeff(1) + mv1.coeff(12)*mv2.coeff(2) - mv1.coeff(14)*mv2.coeff(7);
		mv3.coeffRef(7) += -mv1.coeff(0)*mv2.coeff(11) + mv1.coeff(1)*mv2.coeff(10) - mv1.coeff(3)*mv2.coeff(3) + mv1.coeff(4)*mv2.coeff(2) - mv1.coeff(5)*mv2.coeff(1) - mv1.coeff(6)*mv2.coeff(6) + mv1.coeff(7)*mv2.coeff(5) - mv1.coeff(8)*mv2.coeff(4) - mv1.coeff(10)*mv2.coeff(0);
		mv3.coeffRef(8) += -mv1.coeff(0)*mv2.coeff(12) + mv1.coeff(2)*mv2.coeff(10) + mv1.coeff(3)*mv2.coeff(15) - mv1.coeff(5)*mv2.coeff(13) + mv1.coeff(6)*mv2.coeff(18) - mv1.coeff(8)*mv2.coeff(16) - mv1.coeff(11)*mv2.coeff(0) + mv1.coeff(13)*mv2.coeff(2) + mv1.coeff(14)*mv2.coeff(5);
		mv3.coeffRef(9) += -mv1.coeff(1)*mv2.coeff(12) + mv1.coeff(2)*mv2.coeff(11) - mv1.coeff(3)*mv2.coeff(14) + mv1.coeff(4)*mv2.coeff(13) - mv1.coeff(6)*mv2.coeff(17) + mv1.coeff(7)*mv2.coeff(16) - mv1.coeff(12)*mv2.coeff(0) - mv1.coeff(13)*mv2.coeff(1) - mv1.coeff(14)*mv2.coeff(4);
		mv3.coeffRef(10) += -mv1.coeff(1)*mv2.coeff(19) - mv1.coeff(2)*mv2.coeff(8) + mv1.coeff(3)*mv2.coeff(18) + mv1.coeff(4)*mv2.coeff(6) - mv1.coeff(6)*mv2.coeff(15) - mv1.coeff(7)*mv2.coeff(3) + mv1.coeff(12)*mv2.coeff(7) - mv1.coeff(13)*mv2.coeff(5) + mv1.coeff(14)*mv2.coeff(2);
		mv3.coeffRef(11) +=  mv1.coeff(0)*mv2.coeff(19) - mv1.coeff(2)*mv2.coeff(9) - mv1.coeff(3)*mv2.coeff(17) + mv1.coeff(5)*mv2.coeff(6) + mv1.coeff(6)*mv2.coeff(14) - mv1.coeff(8)*mv2.coeff(3) - mv1.coeff(11)*mv2.coeff(7) + mv1.coeff(13)*mv2.coeff(4) - mv1.coeff(14)*mv2.coeff(1);
		mv3.coeffRef(12) += -mv1.coeff(4)*mv2.coeff(17) - mv1.coeff(5)*mv2.coeff(18) + mv1.coeff(7)*mv2.coeff(14) + mv1.coeff(8)*mv2.coeff(15) - mv1.coeff(10)*mv2.coeff(19) - mv1.coeff(11)*mv2.coeff(8) - mv1.coeff(12)*mv2.coeff(9) + mv1.coeff(13)*mv2.coeff(16) - mv1.coeff(14)*mv2.coeff(13);
		mv3.coeffRef(13) += -mv1.coeff(0)*mv2.coeff(18) + mv1.coeff(1)*mv2.coeff(17) - mv1.coeff(4)*mv2.coeff(9) + mv1.coeff(5)*mv2.coeff(8) - mv1.coeff(6)*mv2.coeff(12) - mv1.coeff(9)*mv2.coeff(3) + mv1.coeff(11)*mv2.coeff(5) - mv1.coeff(12)*mv2.coeff(4) + mv1.coeff(14)*mv2.coeff(0);
		mv3.coeffRef(14) +=  mv1.coeff(2)*mv2.coeff(17) - mv1.coeff(5)*mv2.coeff(19) - mv1.coeff(7)*mv2.coeff(12) + mv1.coeff(9)*mv2.coeff(15) + mv1.coeff(10)*mv2.coeff(18) + mv1.coeff(11)*mv2.coeff(6) - mv1.coeff(12)*mv2.coeff(16) - mv1.coeff(13)*mv2.coeff(9) + mv1.coeff(14)*mv2.coeff(11);
		mv3.coeffRef(15) +=  mv1.coeff(2)*mv2.coeff(18) + mv1.coeff(4)*mv2.coeff(19) - mv1.coeff(8)*mv2.coeff(12) - mv1.coeff(9)*mv2.coeff(14) - mv1.coeff(10)*mv2.coeff(17) + mv1.coeff(11)*mv2.coeff(16) + mv1.coeff(12)*mv2.coeff(6) + mv1.coeff(13)*mv2.coeff(8) - mv1.coeff(14)*mv2.coeff(10);
		mv3.coeffRef(16) +=  mv1.coeff(0)*mv2.coeff(15) - mv1.coeff(1)*mv2.coeff(14) + mv1.coeff(3)*mv2.coeff(12) - mv1.coeff(7)*mv2.coeff(9) + mv1.coeff(8)*mv2.coeff(8) - mv1.coeff(9)*mv2.coeff(6) - mv1.coeff(11)*mv2.coeff(2) + mv1.coeff(12)*mv2.coeff(1) - mv1.coeff(13)*mv2.coeff(0);
		mv3.coeffRef(17) += -mv1.coeff(2)*mv2.coeff(14) + mv1.coeff(4)*mv2.coeff(12) - mv1.coeff(8)*mv2.coeff(19) + mv1.coeff(9)*mv2.coeff(18) - mv1.coeff(10)*mv2.coeff(15) - mv1.coeff(11)*mv2.coeff(3) + mv1.coeff(12)*mv2.coeff(13) - mv1.coeff(13)*mv2.coeff(11) - mv1.coeff(14)*mv2.coeff(9);
		mv3.coeffRef(18) += -mv1.coeff(2)*mv2.coeff(15) + mv1.coeff(5)*mv2.coeff(12) + mv1.coeff(7)*mv2.coeff(19) - mv1.coeff(9)*mv2.coeff(17) + mv1.coeff(10)*mv2.coeff(14) - mv1.coeff(11)*mv2.coeff(13) - mv1.coeff(12)*mv2.coeff(3) + mv1.coeff(13)*mv2.coeff(10) + mv1.coeff(14)*mv2.coeff(8);
		mv3.coeffRef(19) += -mv1.coeff(4)*mv2.coeff(15) + mv1.coeff(5)*mv2.coeff(14) - mv1.coeff(7)*mv2.coeff(18) + mv1.coeff(8)*mv2.coeff(17) - mv1.coeff(10)*mv2.coeff(12) + mv1.coeff(11)*mv2.coeff(11) - mv1.coeff(12)*mv2.coeff(10) - mv1.coeff(13)*mv2.coeff(3) - mv1.coeff(14)*mv2.coeff(6);
	}


	/// \brief Compute the geometric product between two homogeneous multivectors mv1 (grade 4) and mv2 (grade 3). 
	/// \tparam the type of value that we manipulate, either float or double or something else.
	/// \param mv1 - the first homogeneous multivector of grade 4 represented as an Eigen::VectorXd
	/// \param mv2 - the second homogeneous multivector of grade 3 represented as a Eigen::VectorXd
	/// \param mv3 - the result of mv1 mv2 whose grade is 5
	template<typename T>
	void geometric_4_3_5(const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv1, const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv2, Eigen::Matrix<T, Eigen::Dynamic, 1>& mv3){
		mv3.coeffRef(0) += -mv1.coeff(0)*mv2.coeff(9) + mv1.coeff(1)*mv2.coeff(8) - mv1.coeff(2)*mv2.coeff(7) - mv1.coeff(3)*mv2.coeff(6) + mv1.coeff(4)*mv2.coeff(5) - mv1.coeff(5)*mv2.coeff(4) + mv1.coeff(6)*mv2.coeff(3) - mv1.coeff(7)*mv2.coeff(2) + mv1.coeff(8)*mv2.coeff(1) - mv1.coeff(9)*mv2.coeff(0);
		mv3.coeffRef(1) += -mv1.coeff(1)*mv2.coeff(19) + mv1.coeff(3)*mv2.coeff(18) - mv1.coeff(5)*mv2.coeff(16) - mv1.coeff(6)*mv2.coeff(15) + mv1.coeff(8)*mv2.coeff(13) - mv1.coeff(9)*mv2.coeff(11) + mv1.coeff(10)*mv2.coeff(9) - mv1.coeff(12)*mv2.coeff(7) + mv1.coeff(13)*mv2.coeff(5) - mv1.coeff(14)*mv2.coeff(2);
		mv3.coeffRef(2) +=  mv1.coeff(0)*mv2.coeff(19) - mv1.coeff(3)*mv2.coeff(17) + mv1.coeff(4)*mv2.coeff(16) + mv1.coeff(6)*mv2.coeff(14) - mv1.coeff(7)*mv2.coeff(13) + mv1.coeff(9)*mv2.coeff(10) - mv1.coeff(10)*mv2.coeff(8) + mv1.coeff(11)*mv2.coeff(7) - mv1.coeff(13)*mv2.coeff(4) + mv1.coeff(14)*mv2.coeff(1);
		mv3.coeffRef(3) += -mv1.coeff(0)*mv2.coeff(18) + mv1.coeff(1)*mv2.coeff(17) - mv1.coeff(2)*mv2.coeff(16) - mv1.coeff(6)*mv2.coeff(12) + mv1.coeff(7)*mv2.coeff(11) - mv1.coeff(8)*mv2.coeff(10) + mv1.coeff(10)*mv2.coeff(6) - mv1.coeff(11)*mv2.coeff(5) + mv1.coeff(12)*mv2.coeff(4) - mv1.coeff(14)*mv2.coeff(0);
		mv3.coeffRef(4) +=  mv1.coeff(0)*mv2.coeff(15) - mv1.coeff(1)*mv2.coeff(14) + mv1.coeff(2)*mv2.coeff(13) + mv1.coeff(3)*mv2.coeff(12) - mv1.coeff(4)*mv2.coeff(11) + mv1.coeff(5)*mv2.coeff(10) - mv1.coeff(10)*mv2.coeff(3) + mv1.coeff(11)*mv2.coeff(2) - mv1.coeff(12)*mv2.coeff(1) + mv1.coeff(13)*mv2.coeff(0);
		mv3.coeffRef(5) +=  mv1.coeff(2)*mv2.coeff(19) - mv1.coeff(4)*mv2.coeff(18) + mv1.coeff(5)*mv2.coeff(17) + mv1.coeff(7)*mv2.coeff(15) - mv1.coeff(8)*mv2.coeff(14) + mv1.coeff(9)*mv2.coeff(12) - mv1.coeff(11)*mv2.coeff(9) + mv1.coeff(12)*mv2.coeff(8) - mv1.coeff(13)*mv2.coeff(6) + mv1.coeff(14)*mv2.coeff(3);
	}


	/// \brief Compute the geometric product between two homogeneous multivectors mv1 (grade 4) and mv2 (grade 4). 
	/// \tparam the type of value that we manipulate, either float or double or something else.
	/// \param mv1 - the first homogeneous multivector of grade 4 represented as an Eigen::VectorXd
	/// \param mv2 - the second homogeneous multivector of grade 4 represented as a Eigen::VectorXd
	/// \param mv3 - the result of mv1 mv2 whose grade is 2
	template<typename T>
	void geometric_4_4_2(const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv1, const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv2, Eigen::Matrix<T, Eigen::Dynamic, 1>& mv3){
		mv3.coeffRef(0) += -mv1.coeff(0)*mv2.coeff(7) - mv1.coeff(1)*mv2.coeff(8) - mv1.coeff(3)*mv2.coeff(9) + mv1.coeff(6)*mv2.coeff(10) + mv1.coeff(7)*mv2.coeff(0) + mv1.coeff(8)*mv2.coeff(1) + mv1.coeff(9)*mv2.coeff(3) - mv1.coeff(10)*mv2.coeff(6);
		mv3.coeffRef(1) +=  mv1.coeff(0)*mv2.coeff(4) + mv1.coeff(1)*mv2.coeff(5) - mv1.coeff(3)*mv2.coeff(10) - mv1.coeff(4)*mv2.coeff(0) - mv1.coeff(5)*mv2.coeff(1) - mv1.coeff(6)*mv2.coeff(9) + mv1.coeff(9)*mv2.coeff(6) + mv1.coeff(10)*mv2.coeff(3);
		mv3.coeffRef(2) += -mv1.coeff(0)*mv2.coeff(2) + mv1.coeff(1)*mv2.coeff(10) + mv1.coeff(2)*mv2.coeff(0) + mv1.coeff(3)*mv2.coeff(5) - mv1.coeff(5)*mv2.coeff(3) + mv1.coeff(6)*mv2.coeff(8) - mv1.coeff(8)*mv2.coeff(6) - mv1.coeff(10)*mv2.coeff(1);
		mv3.coeffRef(3) += -mv1.coeff(0)*mv2.coeff(10) - mv1.coeff(1)*mv2.coeff(2) + mv1.coeff(2)*mv2.coeff(1) - mv1.coeff(3)*mv2.coeff(4) + mv1.coeff(4)*mv2.coeff(3) - mv1.coeff(6)*mv2.coeff(7) + mv1.coeff(7)*mv2.coeff(6) + mv1.coeff(10)*mv2.coeff(0);
		mv3.coeffRef(4) += -mv1.coeff(0)*mv2.coeff(11) - mv1.coeff(1)*mv2.coeff(12) - mv1.coeff(3)*mv2.coeff(13) - mv1.coeff(6)*mv2.coeff(14) + mv1.coeff(11)*mv2.coeff(0) + mv1.coeff(12)*mv2.coeff(1) + mv1.coeff(13)*mv2.coeff(3) + mv1.coeff(14)*mv2.coeff(6);
		mv3.coeffRef(5) +=  mv1.coeff(3)*mv2.coeff(14) - mv1.coeff(4)*mv2.coeff(7) - mv1.coeff(5)*mv2.coeff(8) - mv1.coeff(6)*mv2.coeff(13) + mv1.coeff(7)*mv2.coeff(4) + mv1.coeff(8)*mv2.coeff(5) + mv1.coeff(13)*mv2.coeff(6) - mv1.coeff(14)*mv2.coeff(3);
		mv3.coeffRef(6) += -mv1.coeff(1)*mv2.coeff(14) + mv1.coeff(2)*mv2.coeff(7) - mv1.coeff(5)*mv2.coeff(9) + mv1.coeff(6)*mv2.coeff(12) - mv1.coeff(7)*mv2.coeff(2) + mv1.coeff(9)*mv2.coeff(5) - mv1.coeff(12)*mv2.coeff(6) + mv1.coeff(14)*mv2.coeff(1);
		mv3.coeffRef(7) +=  mv1.coeff(0)*mv2.coeff(14) + mv1.coeff(2)*mv2.coeff(8) + mv1.coeff(4)*mv2.coeff(9) - mv1.coeff(6)*mv2.coeff(11) - mv1.coeff(8)*mv2.coeff(2) - mv1.coeff(9)*mv2.coeff(4) + mv1.coeff(11)*mv2.coeff(6) - mv1.coeff(14)*mv2.coeff(0);
		mv3.coeffRef(8) += -mv1.coeff(7)*mv2.coeff(11) - mv1.coeff(8)*mv2.coeff(12) - mv1.coeff(9)*mv2.coeff(13) - mv1.coeff(10)*mv2.coeff(14) + mv1.coeff(11)*mv2.coeff(7) + mv1.coeff(12)*mv2.coeff(8) + mv1.coeff(13)*mv2.coeff(9) + mv1.coeff(14)*mv2.coeff(10);
		mv3.coeffRef(9) +=  mv1.coeff(1)*mv2.coeff(13) - mv1.coeff(2)*mv2.coeff(4) - mv1.coeff(3)*mv2.coeff(12) + mv1.coeff(4)*mv2.coeff(2) - mv1.coeff(8)*mv2.coeff(9) + mv1.coeff(9)*mv2.coeff(8) + mv1.coeff(12)*mv2.coeff(3) - mv1.coeff(13)*mv2.coeff(1);
		mv3.coeffRef(10) += -mv1.coeff(0)*mv2.coeff(13) - mv1.coeff(2)*mv2.coeff(5) + mv1.coeff(3)*mv2.coeff(11) + mv1.coeff(5)*mv2.coeff(2) + mv1.coeff(7)*mv2.coeff(9) - mv1.coeff(9)*mv2.coeff(7) - mv1.coeff(11)*mv2.coeff(3) + mv1.coeff(13)*mv2.coeff(0);
		mv3.coeffRef(11) +=  mv1.coeff(4)*mv2.coeff(11) + mv1.coeff(5)*mv2.coeff(12) - mv1.coeff(9)*mv2.coeff(14) + mv1.coeff(10)*mv2.coeff(13) - mv1.coeff(11)*mv2.coeff(4) - mv1.coeff(12)*mv2.coeff(5) - mv1.coeff(13)*mv2.coeff(10) + mv1.coeff(14)*mv2.coeff(9);
		mv3.coeffRef(12) +=  mv1.coeff(0)*mv2.coeff(12) - mv1.coeff(1)*mv2.coeff(11) - mv1.coeff(4)*mv2.coeff(5) + mv1.coeff(5)*mv2.coeff(4) - mv1.coeff(7)*mv2.coeff(8) + mv1.coeff(8)*mv2.coeff(7) + mv1.coeff(11)*mv2.coeff(1) - mv1.coeff(12)*mv2.coeff(0);
		mv3.coeffRef(13) += -mv1.coeff(2)*mv2.coeff(11) + mv1.coeff(5)*mv2.coeff(13) + mv1.coeff(8)*mv2.coeff(14) - mv1.coeff(10)*mv2.coeff(12) + mv1.coeff(11)*mv2.coeff(2) + mv1.coeff(12)*mv2.coeff(10) - mv1.coeff(13)*mv2.coeff(5) - mv1.coeff(14)*mv2.coeff(8);
		mv3.coeffRef(14) += -mv1.coeff(2)*mv2.coeff(12) - mv1.coeff(4)*mv2.coeff(13) - mv1.coeff(7)*mv2.coeff(14) + mv1.coeff(10)*mv2.coeff(11) - mv1.coeff(11)*mv2.coeff(10) + mv1.coeff(12)*mv2.coeff(2) + mv1.coeff(13)*mv2.coeff(4) + mv1.coeff(14)*mv2.coeff(7);
	}


	/// \brief Compute the geometric product between two homogeneous multivectors mv1 (grade 4) and mv2 (grade 4). 
	/// \tparam the type of value that we manipulate, either float or double or something else.
	/// \param mv1 - the first homogeneous multivector of grade 4 represented as an Eigen::VectorXd
	/// \param mv2 - the second homogeneous multivector of grade 4 represented as a Eigen::VectorXd
	/// \param mv3 - the result of mv1 mv2 whose grade is 4
	template<typename T>
	void geometric_4_4_4(const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv1, const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv2, Eigen::Matrix<T, Eigen::Dynamic, 1>& mv3){
		mv3.coeffRef(0) += -mv1.coeff(1)*mv2.coeff(9) + mv1.coeff(3)*mv2.coeff(8) - mv1.coeff(5)*mv2.coeff(6) - mv1.coeff(6)*mv2.coeff(5) + mv1.coeff(8)*mv2.coeff(3) - mv1.coeff(9)*mv2.coeff(1);
		mv3.coeffRef(1) +=  mv1.coeff(0)*mv2.coeff(9) - mv1.coeff(3)*mv2.coeff(7) + mv1.coeff(4)*mv2.coeff(6) + mv1.coeff(6)*mv2.coeff(4) - mv1.coeff(7)*mv2.coeff(3) + mv1.coeff(9)*mv2.coeff(0);
		mv3.coeffRef(2) += -mv1.coeff(3)*mv2.coeff(14) + mv1.coeff(6)*mv2.coeff(13) - mv1.coeff(9)*mv2.coeff(10) - mv1.coeff(10)*mv2.coeff(9) + mv1.coeff(13)*mv2.coeff(6) - mv1.coeff(14)*mv2.coeff(3);
		mv3.coeffRef(3) += -mv1.coeff(0)*mv2.coeff(8) + mv1.coeff(1)*mv2.coeff(7) - mv1.coeff(2)*mv2.coeff(6) - mv1.coeff(6)*mv2.coeff(2) + mv1.coeff(7)*mv2.coeff(1) - mv1.coeff(8)*mv2.coeff(0);
		mv3.coeffRef(4) +=  mv1.coeff(1)*mv2.coeff(14) - mv1.coeff(6)*mv2.coeff(12) + mv1.coeff(8)*mv2.coeff(10) + mv1.coeff(10)*mv2.coeff(8) - mv1.coeff(12)*mv2.coeff(6) + mv1.coeff(14)*mv2.coeff(1);
		mv3.coeffRef(5) += -mv1.coeff(0)*mv2.coeff(14) + mv1.coeff(6)*mv2.coeff(11) - mv1.coeff(7)*mv2.coeff(10) - mv1.coeff(10)*mv2.coeff(7) + mv1.coeff(11)*mv2.coeff(6) - mv1.coeff(14)*mv2.coeff(0);
		mv3.coeffRef(6) +=  mv1.coeff(0)*mv2.coeff(5) - mv1.coeff(1)*mv2.coeff(4) + mv1.coeff(2)*mv2.coeff(3) + mv1.coeff(3)*mv2.coeff(2) - mv1.coeff(4)*mv2.coeff(1) + mv1.coeff(5)*mv2.coeff(0);
		mv3.coeffRef(7) += -mv1.coeff(1)*mv2.coeff(13) + mv1.coeff(3)*mv2.coeff(12) - mv1.coeff(5)*mv2.coeff(10) - mv1.coeff(10)*mv2.coeff(5) + mv1.coeff(12)*mv2.coeff(3) - mv1.coeff(13)*mv2.coeff(1);
		mv3.coeffRef(8) +=  mv1.coeff(0)*mv2.coeff(13) - mv1.coeff(3)*mv2.coeff(11) + mv1.coeff(4)*mv2.coeff(10) + mv1.coeff(10)*mv2.coeff(4) - mv1.coeff(11)*mv2.coeff(3) + mv1.coeff(13)*mv2.coeff(0);
		mv3.coeffRef(9) += -mv1.coeff(0)*mv2.coeff(12) + mv1.coeff(1)*mv2.coeff(11) - mv1.coeff(2)*mv2.coeff(10) - mv1.coeff(10)*mv2.coeff(2) + mv1.coeff(11)*mv2.coeff(1) - mv1.coeff(12)*mv2.coeff(0);
		mv3.coeffRef(10) +=  mv1.coeff(2)*mv2.coeff(9) - mv1.coeff(4)*mv2.coeff(8) + mv1.coeff(5)*mv2.coeff(7) + mv1.coeff(7)*mv2.coeff(5) - mv1.coeff(8)*mv2.coeff(4) + mv1.coeff(9)*mv2.coeff(2);
		mv3.coeffRef(11) +=  mv1.coeff(5)*mv2.coeff(14) - mv1.coeff(8)*mv2.coeff(13) + mv1.coeff(9)*mv2.coeff(12) + mv1.coeff(12)*mv2.coeff(9) - mv1.coeff(13)*mv2.coeff(8) + mv1.coeff(14)*mv2.coeff(5);
		mv3.coeffRef(12) += -mv1.coeff(4)*mv2.coeff(14) + mv1.coeff(7)*mv2.coeff(13) - mv1.coeff(9)*mv2.coeff(11) - mv1.coeff(11)*mv2.coeff(9) + mv1.coeff(13)*mv2.coeff(7) - mv1.coeff(14)*mv2.coeff(4);
		mv3.coeffRef(13) +=  mv1.coeff(2)*mv2.coeff(14) - mv1.coeff(7)*mv2.coeff(12) + mv1.coeff(8)*mv2.coeff(11) + mv1.coeff(11)*mv2.coeff(8) - mv1.coeff(12)*mv2.coeff(7) + mv1.coeff(14)*mv2.coeff(2);
		mv3.coeffRef(14) += -mv1.coeff(2)*mv2.coeff(13) + mv1.coeff(4)*mv2.coeff(12) - mv1.coeff(5)*mv2.coeff(11) - mv1.coeff(11)*mv2.coeff(5) + mv1.coeff(12)*mv2.coeff(4) - mv1.coeff(13)*mv2.coeff(2);
	}


	/// \brief Compute the geometric product between two homogeneous multivectors mv1 (grade 4) and mv2 (grade 5). 
	/// \tparam the type of value that we manipulate, either float or double or something else.
	/// \param mv1 - the first homogeneous multivector of grade 4 represented as an Eigen::VectorXd
	/// \param mv2 - the second homogeneous multivector of grade 5 represented as a Eigen::VectorXd
	/// \param mv3 - the result of mv1 mv2 whose grade is 3
	template<typename T>
	void geometric_4_5_3(const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv1, const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv2, Eigen::Matrix<T, Eigen::Dynamic, 1>& mv3){
		mv3.coeffRef(0) +=  mv1.coeff(3)*mv2.coeff(4) - mv1.coeff(6)*mv2.coeff(3) + mv1.coeff(9)*mv2.coeff(0);
		mv3.coeffRef(1) += -mv1.coeff(1)*mv2.coeff(4) + mv1.coeff(6)*mv2.coeff(2) - mv1.coeff(8)*mv2.coeff(0);
		mv3.coeffRef(2) +=  mv1.coeff(0)*mv2.coeff(4) - mv1.coeff(6)*mv2.coeff(1) + mv1.coeff(7)*mv2.coeff(0);
		mv3.coeffRef(3) +=  mv1.coeff(6)*mv2.coeff(5) - mv1.coeff(10)*mv2.coeff(4) + mv1.coeff(14)*mv2.coeff(0);
		mv3.coeffRef(4) +=  mv1.coeff(1)*mv2.coeff(3) - mv1.coeff(3)*mv2.coeff(2) + mv1.coeff(5)*mv2.coeff(0);
		mv3.coeffRef(5) += -mv1.coeff(0)*mv2.coeff(3) + mv1.coeff(3)*mv2.coeff(1) - mv1.coeff(4)*mv2.coeff(0);
		mv3.coeffRef(6) += -mv1.coeff(3)*mv2.coeff(5) + mv1.coeff(10)*mv2.coeff(3) - mv1.coeff(13)*mv2.coeff(0);
		mv3.coeffRef(7) +=  mv1.coeff(0)*mv2.coeff(2) - mv1.coeff(1)*mv2.coeff(1) + mv1.coeff(2)*mv2.coeff(0);
		mv3.coeffRef(8) +=  mv1.coeff(1)*mv2.coeff(5) - mv1.coeff(10)*mv2.coeff(2) + mv1.coeff(12)*mv2.coeff(0);
		mv3.coeffRef(9) += -mv1.coeff(0)*mv2.coeff(5) + mv1.coeff(10)*mv2.coeff(1) - mv1.coeff(11)*mv2.coeff(0);
		mv3.coeffRef(10) += -mv1.coeff(5)*mv2.coeff(4) + mv1.coeff(8)*mv2.coeff(3) - mv1.coeff(9)*mv2.coeff(2);
		mv3.coeffRef(11) +=  mv1.coeff(4)*mv2.coeff(4) - mv1.coeff(7)*mv2.coeff(3) + mv1.coeff(9)*mv2.coeff(1);
		mv3.coeffRef(12) += -mv1.coeff(9)*mv2.coeff(5) + mv1.coeff(13)*mv2.coeff(4) - mv1.coeff(14)*mv2.coeff(3);
		mv3.coeffRef(13) += -mv1.coeff(2)*mv2.coeff(4) + mv1.coeff(7)*mv2.coeff(2) - mv1.coeff(8)*mv2.coeff(1);
		mv3.coeffRef(14) +=  mv1.coeff(8)*mv2.coeff(5) - mv1.coeff(12)*mv2.coeff(4) + mv1.coeff(14)*mv2.coeff(2);
		mv3.coeffRef(15) += -mv1.coeff(7)*mv2.coeff(5) + mv1.coeff(11)*mv2.coeff(4) - mv1.coeff(14)*mv2.coeff(1);
		mv3.coeffRef(16) +=  mv1.coeff(2)*mv2.coeff(3) - mv1.coeff(4)*mv2.coeff(2) + mv1.coeff(5)*mv2.coeff(1);
		mv3.coeffRef(17) += -mv1.coeff(5)*mv2.coeff(5) + mv1.coeff(12)*mv2.coeff(3) - mv1.coeff(13)*mv2.coeff(2);
		mv3.coeffRef(18) +=  mv1.coeff(4)*mv2.coeff(5) - mv1.coeff(11)*mv2.coeff(3) + mv1.coeff(13)*mv2.coeff(1);
		mv3.coeffRef(19) += -mv1.coeff(2)*mv2.coeff(5) + mv1.coeff(11)*mv2.coeff(2) - mv1.coeff(12)*mv2.coeff(1);
	}


	/// \brief Compute the geometric product between two homogeneous multivectors mv1 (grade 5) and mv2 (grade 2). 
	/// \tparam the type of value that we manipulate, either float or double or something else.
	/// \param mv1 - the first homogeneous multivector of grade 5 represented as an Eigen::VectorXd
	/// \param mv2 - the second homogeneous multivector of grade 2 represented as a Eigen::VectorXd
	/// \param mv3 - the result of mv1 mv2 whose grade is 5
	template<typename T>
	void geometric_5_2_5(const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv1, const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv2, Eigen::Matrix<T, Eigen::Dynamic, 1>& mv3){
		mv3.coeffRef(0) +=  mv1.coeff(0)*mv2.coeff(4) - mv1.coeff(1)*mv2.coeff(3) + mv1.coeff(2)*mv2.coeff(2) - mv1.coeff(3)*mv2.coeff(1) + mv1.coeff(4)*mv2.coeff(0);
		mv3.coeffRef(1) +=  mv1.coeff(0)*mv2.coeff(14) - mv1.coeff(2)*mv2.coeff(12) + mv1.coeff(3)*mv2.coeff(10) - mv1.coeff(4)*mv2.coeff(7) + mv1.coeff(5)*mv2.coeff(3);
		mv3.coeffRef(2) += -mv1.coeff(0)*mv2.coeff(13) + mv1.coeff(1)*mv2.coeff(12) - mv1.coeff(3)*mv2.coeff(9) + mv1.coeff(4)*mv2.coeff(6) - mv1.coeff(5)*mv2.coeff(2);
		mv3.coeffRef(3) +=  mv1.coeff(0)*mv2.coeff(11) - mv1.coeff(1)*mv2.coeff(10) + mv1.coeff(2)*mv2.coeff(9) - mv1.coeff(4)*mv2.coeff(5) + mv1.coeff(5)*mv2.coeff(1);
		mv3.coeffRef(4) += -mv1.coeff(0)*mv2.coeff(8) + mv1.coeff(1)*mv2.coeff(7) - mv1.coeff(2)*mv2.coeff(6) + mv1.coeff(3)*mv2.coeff(5) - mv1.coeff(5)*mv2.coeff(0);
		mv3.coeffRef(5) += -mv1.coeff(1)*mv2.coeff(14) + mv1.coeff(2)*mv2.coeff(13) - mv1.coeff(3)*mv2.coeff(11) + mv1.coeff(4)*mv2.coeff(8) - mv1.coeff(5)*mv2.coeff(4);
	}


	/// \brief Compute the geometric product between two homogeneous multivectors mv1 (grade 5) and mv2 (grade 3). 
	/// \tparam the type of value that we manipulate, either float or double or something else.
	/// \param mv1 - the first homogeneous multivector of grade 5 represented as an Eigen::VectorXd
	/// \param mv2 - the second homogeneous multivector of grade 3 represented as a Eigen::VectorXd
	/// \param mv3 - the result of mv1 mv2 whose grade is 4
	template<typename T>
	void geometric_5_3_4(const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv1, const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv2, Eigen::Matrix<T, Eigen::Dynamic, 1>& mv3){
		mv3.coeffRef(0) += -mv1.coeff(0)*mv2.coeff(9) + mv1.coeff(2)*mv2.coeff(7) - mv1.coeff(3)*mv2.coeff(5) + mv1.coeff(4)*mv2.coeff(2);
		mv3.coeffRef(1) +=  mv1.coeff(0)*mv2.coeff(8) - mv1.coeff(1)*mv2.coeff(7) + mv1.coeff(3)*mv2.coeff(4) - mv1.coeff(4)*mv2.coeff(1);
		mv3.coeffRef(2) += -mv1.coeff(0)*mv2.coeff(19) + mv1.coeff(3)*mv2.coeff(16) - mv1.coeff(4)*mv2.coeff(13) + mv1.coeff(5)*mv2.coeff(7);
		mv3.coeffRef(3) += -mv1.coeff(0)*mv2.coeff(6) + mv1.coeff(1)*mv2.coeff(5) - mv1.coeff(2)*mv2.coeff(4) + mv1.coeff(4)*mv2.coeff(0);
		mv3.coeffRef(4) +=  mv1.coeff(0)*mv2.coeff(18) - mv1.coeff(2)*mv2.coeff(16) + mv1.coeff(4)*mv2.coeff(11) - mv1.coeff(5)*mv2.coeff(5);
		mv3.coeffRef(5) += -mv1.coeff(0)*mv2.coeff(17) + mv1.coeff(1)*mv2.coeff(16) - mv1.coeff(4)*mv2.coeff(10) + mv1.coeff(5)*mv2.coeff(4);
		mv3.coeffRef(6) +=  mv1.coeff(0)*mv2.coeff(3) - mv1.coeff(1)*mv2.coeff(2) + mv1.coeff(2)*mv2.coeff(1) - mv1.coeff(3)*mv2.coeff(0);
		mv3.coeffRef(7) += -mv1.coeff(0)*mv2.coeff(15) + mv1.coeff(2)*mv2.coeff(13) - mv1.coeff(3)*mv2.coeff(11) + mv1.coeff(5)*mv2.coeff(2);
		mv3.coeffRef(8) +=  mv1.coeff(0)*mv2.coeff(14) - mv1.coeff(1)*mv2.coeff(13) + mv1.coeff(3)*mv2.coeff(10) - mv1.coeff(5)*mv2.coeff(1);
		mv3.coeffRef(9) += -mv1.coeff(0)*mv2.coeff(12) + mv1.coeff(1)*mv2.coeff(11) - mv1.coeff(2)*mv2.coeff(10) + mv1.coeff(5)*mv2.coeff(0);
		mv3.coeffRef(10) +=  mv1.coeff(1)*mv2.coeff(9) - mv1.coeff(2)*mv2.coeff(8) + mv1.coeff(3)*mv2.coeff(6) - mv1.coeff(4)*mv2.coeff(3);
		mv3.coeffRef(11) +=  mv1.coeff(2)*mv2.coeff(19) - mv1.coeff(3)*mv2.coeff(18) + mv1.coeff(4)*mv2.coeff(15) - mv1.coeff(5)*mv2.coeff(9);
		mv3.coeffRef(12) += -mv1.coeff(1)*mv2.coeff(19) + mv1.coeff(3)*mv2.coeff(17) - mv1.coeff(4)*mv2.coeff(14) + mv1.coeff(5)*mv2.coeff(8);
		mv3.coeffRef(13) +=  mv1.coeff(1)*mv2.coeff(18) - mv1.coeff(2)*mv2.coeff(17) + mv1.coeff(4)*mv2.coeff(12) - mv1.coeff(5)*mv2.coeff(6);
		mv3.coeffRef(14) += -mv1.coeff(1)*mv2.coeff(15) + mv1.coeff(2)*mv2.coeff(14) - mv1.coeff(3)*mv2.coeff(12) + mv1.coeff(5)*mv2.coeff(3);
	}


	/// \brief Compute the geometric product between two homogeneous multivectors mv1 (grade 5) and mv2 (grade 4). 
	/// \tparam the type of value that we manipulate, either float or double or something else.
	/// \param mv1 - the first homogeneous multivector of grade 5 represented as an Eigen::VectorXd
	/// \param mv2 - the second homogeneous multivector of grade 4 represented as a Eigen::VectorXd
	/// \param mv3 - the result of mv1 mv2 whose grade is 3
	template<typename T>
	void geometric_5_4_3(const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv1, const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv2, Eigen::Matrix<T, Eigen::Dynamic, 1>& mv3){
		mv3.coeffRef(0) += -mv1.coeff(0)*mv2.coeff(9) + mv1.coeff(3)*mv2.coeff(6) - mv1.coeff(4)*mv2.coeff(3);
		mv3.coeffRef(1) +=  mv1.coeff(0)*mv2.coeff(8) - mv1.coeff(2)*mv2.coeff(6) + mv1.coeff(4)*mv2.coeff(1);
		mv3.coeffRef(2) += -mv1.coeff(0)*mv2.coeff(7) + mv1.coeff(1)*mv2.coeff(6) - mv1.coeff(4)*mv2.coeff(0);
		mv3.coeffRef(3) += -mv1.coeff(0)*mv2.coeff(14) + mv1.coeff(4)*mv2.coeff(10) - mv1.coeff(5)*mv2.coeff(6);
		mv3.coeffRef(4) += -mv1.coeff(0)*mv2.coeff(5) + mv1.coeff(2)*mv2.coeff(3) - mv1.coeff(3)*mv2.coeff(1);
		mv3.coeffRef(5) +=  mv1.coeff(0)*mv2.coeff(4) - mv1.coeff(1)*mv2.coeff(3) + mv1.coeff(3)*mv2.coeff(0);
		mv3.coeffRef(6) +=  mv1.coeff(0)*mv2.coeff(13) - mv1.coeff(3)*mv2.coeff(10) + mv1.coeff(5)*mv2.coeff(3);
		mv3.coeffRef(7) += -mv1.coeff(0)*mv2.coeff(2) + mv1.coeff(1)*mv2.coeff(1) - mv1.coeff(2)*mv2.coeff(0);
		mv3.coeffRef(8) += -mv1.coeff(0)*mv2.coeff(12) + mv1.coeff(2)*mv2.coeff(10) - mv1.coeff(5)*mv2.coeff(1);
		mv3.coeffRef(9) +=  mv1.coeff(0)*mv2.coeff(11) - mv1.coeff(1)*mv2.coeff(10) + mv1.coeff(5)*mv2.coeff(0);
		mv3.coeffRef(10) +=  mv1.coeff(2)*mv2.coeff(9) - mv1.coeff(3)*mv2.coeff(8) + mv1.coeff(4)*mv2.coeff(5);
		mv3.coeffRef(11) += -mv1.coeff(1)*mv2.coeff(9) + mv1.coeff(3)*mv2.coeff(7) - mv1.coeff(4)*mv2.coeff(4);
		mv3.coeffRef(12) +=  mv1.coeff(3)*mv2.coeff(14) - mv1.coeff(4)*mv2.coeff(13) + mv1.coeff(5)*mv2.coeff(9);
		mv3.coeffRef(13) +=  mv1.coeff(1)*mv2.coeff(8) - mv1.coeff(2)*mv2.coeff(7) + mv1.coeff(4)*mv2.coeff(2);
		mv3.coeffRef(14) += -mv1.coeff(2)*mv2.coeff(14) + mv1.coeff(4)*mv2.coeff(12) - mv1.coeff(5)*mv2.coeff(8);
		mv3.coeffRef(15) +=  mv1.coeff(1)*mv2.coeff(14) - mv1.coeff(4)*mv2.coeff(11) + mv1.coeff(5)*mv2.coeff(7);
		mv3.coeffRef(16) += -mv1.coeff(1)*mv2.coeff(5) + mv1.coeff(2)*mv2.coeff(4) - mv1.coeff(3)*mv2.coeff(2);
		mv3.coeffRef(17) +=  mv1.coeff(2)*mv2.coeff(13) - mv1.coeff(3)*mv2.coeff(12) + mv1.coeff(5)*mv2.coeff(5);
		mv3.coeffRef(18) += -mv1.coeff(1)*mv2.coeff(13) + mv1.coeff(3)*mv2.coeff(11) - mv1.coeff(5)*mv2.coeff(4);
		mv3.coeffRef(19) +=  mv1.coeff(1)*mv2.coeff(12) - mv1.coeff(2)*mv2.coeff(11) + mv1.coeff(5)*mv2.coeff(2);
	}


	/// \brief Compute the geometric product between two homogeneous multivectors mv1 (grade 5) and mv2 (grade 5). 
	/// \tparam the type of value that we manipulate, either float or double or something else.
	/// \param mv1 - the first homogeneous multivector of grade 5 represented as an Eigen::VectorXd
	/// \param mv2 - the second homogeneous multivector of grade 5 represented as a Eigen::VectorXd
	/// \param mv3 - the result of mv1 mv2 whose grade is 2
	template<typename T>
	void geometric_5_5_2(const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv1, const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv2, Eigen::Matrix<T, Eigen::Dynamic, 1>& mv3){
		mv3.coeffRef(0) +=  mv1.coeff(0)*mv2.coeff(4) - mv1.coeff(4)*mv2.coeff(0);
		mv3.coeffRef(1) += -mv1.coeff(0)*mv2.coeff(3) + mv1.coeff(3)*mv2.coeff(0);
		mv3.coeffRef(2) +=  mv1.coeff(0)*mv2.coeff(2) - mv1.coeff(2)*mv2.coeff(0);
		mv3.coeffRef(3) += -mv1.coeff(0)*mv2.coeff(1) + mv1.coeff(1)*mv2.coeff(0);
		mv3.coeffRef(4) +=  mv1.coeff(0)*mv2.coeff(5) - mv1.coeff(5)*mv2.coeff(0);
		mv3.coeffRef(5) += -mv1.coeff(3)*mv2.coeff(4) + mv1.coeff(4)*mv2.coeff(3);
		mv3.coeffRef(6) +=  mv1.coeff(2)*mv2.coeff(4) - mv1.coeff(4)*mv2.coeff(2);
		mv3.coeffRef(7) += -mv1.coeff(1)*mv2.coeff(4) + mv1.coeff(4)*mv2.coeff(1);
		mv3.coeffRef(8) += -mv1.coeff(4)*mv2.coeff(5) + mv1.coeff(5)*mv2.coeff(4);
		mv3.coeffRef(9) += -mv1.coeff(2)*mv2.coeff(3) + mv1.coeff(3)*mv2.coeff(2);
		mv3.coeffRef(10) +=  mv1.coeff(1)*mv2.coeff(3) - mv1.coeff(3)*mv2.coeff(1);
		mv3.coeffRef(11) +=  mv1.coeff(3)*mv2.coeff(5) - mv1.coeff(5)*mv2.coeff(3);
		mv3.coeffRef(12) += -mv1.coeff(1)*mv2.coeff(2) + mv1.coeff(2)*mv2.coeff(1);
		mv3.coeffRef(13) += -mv1.coeff(2)*mv2.coeff(5) + mv1.coeff(5)*mv2.coeff(2);
		mv3.coeffRef(14) +=  mv1.coeff(1)*mv2.coeff(5) - mv1.coeff(5)*mv2.coeff(1);
	}


	
    template<typename T>
	std::array<std::array<std::array<std::function<void(const Eigen::Matrix<T, Eigen::Dynamic, 1> & , const Eigen::Matrix<T, Eigen::Dynamic, 1> & , Eigen::Matrix<T, Eigen::Dynamic, 1>&)>, 7>, 7>, 7> geometricFunctionsContainer =
	{{
		{{
			{{{},{},{},{},{},{},{}}},
			{{{},{},{},{},{},{},{}}},
			{{{},{},{},{},{},{},{}}},
			{{{},{},{},{},{},{},{}}},
			{{{},{},{},{},{},{},{}}},
			{{{},{},{},{},{},{},{}}},
			{{{},{},{},{},{},{},{}}}
		}},
		{{
			{{{},{},{},{},{},{},{}}},
			{{{},{},{},{},{},{},{}}},
			{{{},{},{},{},{},{},{}}},
			{{{},{},{},{},{},{},{}}},
			{{{},{},{},{},{},{},{}}},
			{{{},{},{},{},{},{},{}}},
			{{{},{},{},{},{},{},{}}}
		}},
		{{
			{{{},{},{},{},{},{},{}}},
			{{{},{},{},{},{},{},{}}},
			{{{},{},geometric_2_2_2<T>,{},{},{},{}}},
			{{{},{},{},geometric_3_3_2<T>,{},{},{}}},
			{{{},{},{},{},geometric_4_4_2<T>,{},{}}},
			{{{},{},{},{},{},geometric_5_5_2<T>,{}}},
			{{{},{},{},{},{},{},{}}}
		}},
		{{
			{{{},{},{},{},{},{},{}}},
			{{{},{},{},{},{},{},{}}},
			{{{},{},{},geometric_2_3_3<T>,{},{},{}}},
			{{{},{},geometric_3_2_3<T>,{},geometric_3_4_3<T>,{},{}}},
			{{{},{},{},geometric_4_3_3<T>,{},geometric_4_5_3<T>,{}}},
			{{{},{},{},{},geometric_5_4_3<T>,{},{}}},
			{{{},{},{},{},{},{},{}}}
		}},
		{{
			{{{},{},{},{},{},{},{}}},
			{{{},{},{},{},{},{},{}}},
			{{{},{},{},{},geometric_2_4_4<T>,{},{}}},
			{{{},{},{},geometric_3_3_4<T>,{},geometric_3_5_4<T>,{}}},
			{{{},{},geometric_4_2_4<T>,{},geometric_4_4_4<T>,{},{}}},
			{{{},{},{},geometric_5_3_4<T>,{},{},{}}},
			{{{},{},{},{},{},{},{}}}
		}},
		{{
			{{{},{},{},{},{},{},{}}},
			{{{},{},{},{},{},{},{}}},
			{{{},{},{},{},{},geometric_2_5_5<T>,{}}},
			{{{},{},{},{},geometric_3_4_5<T>,{},{}}},
			{{{},{},{},geometric_4_3_5<T>,{},{},{}}},
			{{{},{},geometric_5_2_5<T>,{},{},{},{}}},
			{{{},{},{},{},{},{},{}}}
		}},
		{{
			{{{},{},{},{},{},{},{}}},
			{{{},{},{},{},{},{},{}}},
			{{{},{},{},{},{},{},{}}},
			{{{},{},{},{},{},{},{}}},
			{{{},{},{},{},{},{},{}}},
			{{{},{},{},{},{},{},{}}},
			{{{},{},{},{},{},{},{}}}
		}}
	}};

}/// End of Namespace

#endif // C4GA_GEOMETRIC_PRODUCT_EXPLICIT_HPP__