// Copyright (c) 2018 by University Paris-Est Marne-la-Vallee
// OuterExplicit.hpp
// This file is part of the Garamon for c4ga.
// Authors: Stephane Breuils and Vincent Nozick
// Contact: vincent.nozick@u-pem.fr
//
// Licence MIT
// A a copy of the MIT License is given along with this program

/// \file OuterExplicit.hpp
/// \author Stephane Breuils, Vincent Nozick
/// \brief Explicit precomputed per grades outer product.


#ifndef C4GA_OUTER_PRODUCT_EXPLICIT_HPP__
#define C4GA_OUTER_PRODUCT_EXPLICIT_HPP__
#pragma once

#include <Eigen/Core>

#include "c4ga/Mvec.hpp"
#include "c4ga/Outer.hpp"


/*!
 * @namespace c4ga
 */
namespace c4ga {
    template<typename T> class Mvec;

    /// \brief Compute the outer product between two homogeneous multivectors mv1 (grade 0) and mv2 (grade 0). 
	/// \tparam the type of value that we manipulate, either float or double or something.
	/// \param mv1 - the first homogeneous multivector of grade 0 represented as an Eigen::VectorXd
	/// \param mv2 - the second homogeneous multivector of grade 0 represented as a Eigen::VectorXd
	/// \param mv3 - the result of mv1^mv2, which is also a homogeneous multivector of grade 0
	template<typename T>
	void outer_0_0(const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv1, const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv2, Eigen::Matrix<T, Eigen::Dynamic, 1>& mv3){
		mv3.coeffRef(0) +=  mv1.coeff(0)*mv2.coeff(0);
	}


	/// \brief Compute the outer product between two homogeneous multivectors mv1 (grade 0) and mv2 (grade 1). 
	/// \tparam the type of value that we manipulate, either float or double or something.
	/// \param mv1 - the first homogeneous multivector of grade 0 represented as an Eigen::VectorXd
	/// \param mv2 - the second homogeneous multivector of grade 1 represented as a Eigen::VectorXd
	/// \param mv3 - the result of mv1^mv2, which is also a homogeneous multivector of grade 1
	template<typename T>
	void outer_0_1(const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv1, const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv2, Eigen::Matrix<T, Eigen::Dynamic, 1>& mv3){
		mv3 += mv1.coeff(0)*mv2;
	}


	/// \brief Compute the outer product between two homogeneous multivectors mv1 (grade 0) and mv2 (grade 2). 
	/// \tparam the type of value that we manipulate, either float or double or something.
	/// \param mv1 - the first homogeneous multivector of grade 0 represented as an Eigen::VectorXd
	/// \param mv2 - the second homogeneous multivector of grade 2 represented as a Eigen::VectorXd
	/// \param mv3 - the result of mv1^mv2, which is also a homogeneous multivector of grade 2
	template<typename T>
	void outer_0_2(const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv1, const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv2, Eigen::Matrix<T, Eigen::Dynamic, 1>& mv3){
		mv3 += mv1.coeff(0)*mv2;
	}


	/// \brief Compute the outer product between two homogeneous multivectors mv1 (grade 0) and mv2 (grade 3). 
	/// \tparam the type of value that we manipulate, either float or double or something.
	/// \param mv1 - the first homogeneous multivector of grade 0 represented as an Eigen::VectorXd
	/// \param mv2 - the second homogeneous multivector of grade 3 represented as a Eigen::VectorXd
	/// \param mv3 - the result of mv1^mv2, which is also a homogeneous multivector of grade 3
	template<typename T>
	void outer_0_3(const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv1, const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv2, Eigen::Matrix<T, Eigen::Dynamic, 1>& mv3){
		mv3 += mv1.coeff(0)*mv2;
	}


	/// \brief Compute the outer product between two homogeneous multivectors mv1 (grade 0) and mv2 (grade 4). 
	/// \tparam the type of value that we manipulate, either float or double or something.
	/// \param mv1 - the first homogeneous multivector of grade 0 represented as an Eigen::VectorXd
	/// \param mv2 - the second homogeneous multivector of grade 4 represented as a Eigen::VectorXd
	/// \param mv3 - the result of mv1^mv2, which is also a homogeneous multivector of grade 4
	template<typename T>
	void outer_0_4(const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv1, const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv2, Eigen::Matrix<T, Eigen::Dynamic, 1>& mv3){
		mv3 += mv1.coeff(0)*mv2;
	}


	/// \brief Compute the outer product between two homogeneous multivectors mv1 (grade 0) and mv2 (grade 5). 
	/// \tparam the type of value that we manipulate, either float or double or something.
	/// \param mv1 - the first homogeneous multivector of grade 0 represented as an Eigen::VectorXd
	/// \param mv2 - the second homogeneous multivector of grade 5 represented as a Eigen::VectorXd
	/// \param mv3 - the result of mv1^mv2, which is also a homogeneous multivector of grade 5
	template<typename T>
	void outer_0_5(const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv1, const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv2, Eigen::Matrix<T, Eigen::Dynamic, 1>& mv3){
		mv3 += mv1.coeff(0)*mv2;
	}


	/// \brief Compute the outer product between two homogeneous multivectors mv1 (grade 0) and mv2 (grade 6). 
	/// \tparam the type of value that we manipulate, either float or double or something.
	/// \param mv1 - the first homogeneous multivector of grade 0 represented as an Eigen::VectorXd
	/// \param mv2 - the second homogeneous multivector of grade 6 represented as a Eigen::VectorXd
	/// \param mv3 - the result of mv1^mv2, which is also a homogeneous multivector of grade 6
	template<typename T>
	void outer_0_6(const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv1, const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv2, Eigen::Matrix<T, Eigen::Dynamic, 1>& mv3){
		mv3 += mv1.coeff(0)*mv2;
	}


	/// \brief Compute the outer product between two homogeneous multivectors mv1 (grade 1) and mv2 (grade 0). 
	/// \tparam the type of value that we manipulate, either float or double or something.
	/// \param mv1 - the first homogeneous multivector of grade 1 represented as an Eigen::VectorXd
	/// \param mv2 - the second homogeneous multivector of grade 0 represented as a Eigen::VectorXd
	/// \param mv3 - the result of mv1^mv2, which is also a homogeneous multivector of grade 1
	template<typename T>
	void outer_1_0(const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv1, const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv2, Eigen::Matrix<T, Eigen::Dynamic, 1>& mv3){
		mv3 += mv1*mv2.coeff(0);
	}


	/// \brief Compute the outer product between two homogeneous multivectors mv1 (grade 1) and mv2 (grade 1). 
	/// \tparam the type of value that we manipulate, either float or double or something.
	/// \param mv1 - the first homogeneous multivector of grade 1 represented as an Eigen::VectorXd
	/// \param mv2 - the second homogeneous multivector of grade 1 represented as a Eigen::VectorXd
	/// \param mv3 - the result of mv1^mv2, which is also a homogeneous multivector of grade 2
	template<typename T>
	void outer_1_1(const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv1, const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv2, Eigen::Matrix<T, Eigen::Dynamic, 1>& mv3){
		mv3.coeffRef(0) +=  mv1.coeff(0)*mv2.coeff(1) - mv1.coeff(1)*mv2.coeff(0);
		mv3.coeffRef(1) +=  mv1.coeff(0)*mv2.coeff(2) - mv1.coeff(2)*mv2.coeff(0);
		mv3.coeffRef(2) +=  mv1.coeff(0)*mv2.coeff(3) - mv1.coeff(3)*mv2.coeff(0);
		mv3.coeffRef(3) +=  mv1.coeff(0)*mv2.coeff(4) - mv1.coeff(4)*mv2.coeff(0);
		mv3.coeffRef(4) +=  mv1.coeff(0)*mv2.coeff(5) - mv1.coeff(5)*mv2.coeff(0);
		mv3.coeffRef(5) +=  mv1.coeff(1)*mv2.coeff(2) - mv1.coeff(2)*mv2.coeff(1);
		mv3.coeffRef(6) +=  mv1.coeff(1)*mv2.coeff(3) - mv1.coeff(3)*mv2.coeff(1);
		mv3.coeffRef(7) +=  mv1.coeff(1)*mv2.coeff(4) - mv1.coeff(4)*mv2.coeff(1);
		mv3.coeffRef(8) +=  mv1.coeff(1)*mv2.coeff(5) - mv1.coeff(5)*mv2.coeff(1);
		mv3.coeffRef(9) +=  mv1.coeff(2)*mv2.coeff(3) - mv1.coeff(3)*mv2.coeff(2);
		mv3.coeffRef(10) +=  mv1.coeff(2)*mv2.coeff(4) - mv1.coeff(4)*mv2.coeff(2);
		mv3.coeffRef(11) +=  mv1.coeff(2)*mv2.coeff(5) - mv1.coeff(5)*mv2.coeff(2);
		mv3.coeffRef(12) +=  mv1.coeff(3)*mv2.coeff(4) - mv1.coeff(4)*mv2.coeff(3);
		mv3.coeffRef(13) +=  mv1.coeff(3)*mv2.coeff(5) - mv1.coeff(5)*mv2.coeff(3);
		mv3.coeffRef(14) +=  mv1.coeff(4)*mv2.coeff(5) - mv1.coeff(5)*mv2.coeff(4);
	}


	/// \brief Compute the outer product between two homogeneous multivectors mv1 (grade 1) and mv2 (grade 2). 
	/// \tparam the type of value that we manipulate, either float or double or something.
	/// \param mv1 - the first homogeneous multivector of grade 1 represented as an Eigen::VectorXd
	/// \param mv2 - the second homogeneous multivector of grade 2 represented as a Eigen::VectorXd
	/// \param mv3 - the result of mv1^mv2, which is also a homogeneous multivector of grade 3
	template<typename T>
	void outer_1_2(const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv1, const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv2, Eigen::Matrix<T, Eigen::Dynamic, 1>& mv3){
		mv3.coeffRef(0) +=  mv1.coeff(0)*mv2.coeff(5) - mv1.coeff(1)*mv2.coeff(1) + mv1.coeff(2)*mv2.coeff(0);
		mv3.coeffRef(1) +=  mv1.coeff(0)*mv2.coeff(6) - mv1.coeff(1)*mv2.coeff(2) + mv1.coeff(3)*mv2.coeff(0);
		mv3.coeffRef(2) +=  mv1.coeff(0)*mv2.coeff(7) - mv1.coeff(1)*mv2.coeff(3) + mv1.coeff(4)*mv2.coeff(0);
		mv3.coeffRef(3) +=  mv1.coeff(0)*mv2.coeff(8) - mv1.coeff(1)*mv2.coeff(4) + mv1.coeff(5)*mv2.coeff(0);
		mv3.coeffRef(4) +=  mv1.coeff(0)*mv2.coeff(9) - mv1.coeff(2)*mv2.coeff(2) + mv1.coeff(3)*mv2.coeff(1);
		mv3.coeffRef(5) +=  mv1.coeff(0)*mv2.coeff(10) - mv1.coeff(2)*mv2.coeff(3) + mv1.coeff(4)*mv2.coeff(1);
		mv3.coeffRef(6) +=  mv1.coeff(0)*mv2.coeff(11) - mv1.coeff(2)*mv2.coeff(4) + mv1.coeff(5)*mv2.coeff(1);
		mv3.coeffRef(7) +=  mv1.coeff(0)*mv2.coeff(12) - mv1.coeff(3)*mv2.coeff(3) + mv1.coeff(4)*mv2.coeff(2);
		mv3.coeffRef(8) +=  mv1.coeff(0)*mv2.coeff(13) - mv1.coeff(3)*mv2.coeff(4) + mv1.coeff(5)*mv2.coeff(2);
		mv3.coeffRef(9) +=  mv1.coeff(0)*mv2.coeff(14) - mv1.coeff(4)*mv2.coeff(4) + mv1.coeff(5)*mv2.coeff(3);
		mv3.coeffRef(10) +=  mv1.coeff(1)*mv2.coeff(9) - mv1.coeff(2)*mv2.coeff(6) + mv1.coeff(3)*mv2.coeff(5);
		mv3.coeffRef(11) +=  mv1.coeff(1)*mv2.coeff(10) - mv1.coeff(2)*mv2.coeff(7) + mv1.coeff(4)*mv2.coeff(5);
		mv3.coeffRef(12) +=  mv1.coeff(1)*mv2.coeff(11) - mv1.coeff(2)*mv2.coeff(8) + mv1.coeff(5)*mv2.coeff(5);
		mv3.coeffRef(13) +=  mv1.coeff(1)*mv2.coeff(12) - mv1.coeff(3)*mv2.coeff(7) + mv1.coeff(4)*mv2.coeff(6);
		mv3.coeffRef(14) +=  mv1.coeff(1)*mv2.coeff(13) - mv1.coeff(3)*mv2.coeff(8) + mv1.coeff(5)*mv2.coeff(6);
		mv3.coeffRef(15) +=  mv1.coeff(1)*mv2.coeff(14) - mv1.coeff(4)*mv2.coeff(8) + mv1.coeff(5)*mv2.coeff(7);
		mv3.coeffRef(16) +=  mv1.coeff(2)*mv2.coeff(12) - mv1.coeff(3)*mv2.coeff(10) + mv1.coeff(4)*mv2.coeff(9);
		mv3.coeffRef(17) +=  mv1.coeff(2)*mv2.coeff(13) - mv1.coeff(3)*mv2.coeff(11) + mv1.coeff(5)*mv2.coeff(9);
		mv3.coeffRef(18) +=  mv1.coeff(2)*mv2.coeff(14) - mv1.coeff(4)*mv2.coeff(11) + mv1.coeff(5)*mv2.coeff(10);
		mv3.coeffRef(19) +=  mv1.coeff(3)*mv2.coeff(14) - mv1.coeff(4)*mv2.coeff(13) + mv1.coeff(5)*mv2.coeff(12);
	}


	/// \brief Compute the outer product between two homogeneous multivectors mv1 (grade 1) and mv2 (grade 3). 
	/// \tparam the type of value that we manipulate, either float or double or something.
	/// \param mv1 - the first homogeneous multivector of grade 1 represented as an Eigen::VectorXd
	/// \param mv2 - the second homogeneous multivector of grade 3 represented as a Eigen::VectorXd
	/// \param mv3 - the result of mv1^mv2, which is also a homogeneous multivector of grade 4
	template<typename T>
	void outer_1_3(const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv1, const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv2, Eigen::Matrix<T, Eigen::Dynamic, 1>& mv3){
		mv3.coeffRef(0) +=  mv1.coeff(0)*mv2.coeff(10) - mv1.coeff(1)*mv2.coeff(4) + mv1.coeff(2)*mv2.coeff(1) - mv1.coeff(3)*mv2.coeff(0);
		mv3.coeffRef(1) +=  mv1.coeff(0)*mv2.coeff(11) - mv1.coeff(1)*mv2.coeff(5) + mv1.coeff(2)*mv2.coeff(2) - mv1.coeff(4)*mv2.coeff(0);
		mv3.coeffRef(2) +=  mv1.coeff(0)*mv2.coeff(12) - mv1.coeff(1)*mv2.coeff(6) + mv1.coeff(2)*mv2.coeff(3) - mv1.coeff(5)*mv2.coeff(0);
		mv3.coeffRef(3) +=  mv1.coeff(0)*mv2.coeff(13) - mv1.coeff(1)*mv2.coeff(7) + mv1.coeff(3)*mv2.coeff(2) - mv1.coeff(4)*mv2.coeff(1);
		mv3.coeffRef(4) +=  mv1.coeff(0)*mv2.coeff(14) - mv1.coeff(1)*mv2.coeff(8) + mv1.coeff(3)*mv2.coeff(3) - mv1.coeff(5)*mv2.coeff(1);
		mv3.coeffRef(5) +=  mv1.coeff(0)*mv2.coeff(15) - mv1.coeff(1)*mv2.coeff(9) + mv1.coeff(4)*mv2.coeff(3) - mv1.coeff(5)*mv2.coeff(2);
		mv3.coeffRef(6) +=  mv1.coeff(0)*mv2.coeff(16) - mv1.coeff(2)*mv2.coeff(7) + mv1.coeff(3)*mv2.coeff(5) - mv1.coeff(4)*mv2.coeff(4);
		mv3.coeffRef(7) +=  mv1.coeff(0)*mv2.coeff(17) - mv1.coeff(2)*mv2.coeff(8) + mv1.coeff(3)*mv2.coeff(6) - mv1.coeff(5)*mv2.coeff(4);
		mv3.coeffRef(8) +=  mv1.coeff(0)*mv2.coeff(18) - mv1.coeff(2)*mv2.coeff(9) + mv1.coeff(4)*mv2.coeff(6) - mv1.coeff(5)*mv2.coeff(5);
		mv3.coeffRef(9) +=  mv1.coeff(0)*mv2.coeff(19) - mv1.coeff(3)*mv2.coeff(9) + mv1.coeff(4)*mv2.coeff(8) - mv1.coeff(5)*mv2.coeff(7);
		mv3.coeffRef(10) +=  mv1.coeff(1)*mv2.coeff(16) - mv1.coeff(2)*mv2.coeff(13) + mv1.coeff(3)*mv2.coeff(11) - mv1.coeff(4)*mv2.coeff(10);
		mv3.coeffRef(11) +=  mv1.coeff(1)*mv2.coeff(17) - mv1.coeff(2)*mv2.coeff(14) + mv1.coeff(3)*mv2.coeff(12) - mv1.coeff(5)*mv2.coeff(10);
		mv3.coeffRef(12) +=  mv1.coeff(1)*mv2.coeff(18) - mv1.coeff(2)*mv2.coeff(15) + mv1.coeff(4)*mv2.coeff(12) - mv1.coeff(5)*mv2.coeff(11);
		mv3.coeffRef(13) +=  mv1.coeff(1)*mv2.coeff(19) - mv1.coeff(3)*mv2.coeff(15) + mv1.coeff(4)*mv2.coeff(14) - mv1.coeff(5)*mv2.coeff(13);
		mv3.coeffRef(14) +=  mv1.coeff(2)*mv2.coeff(19) - mv1.coeff(3)*mv2.coeff(18) + mv1.coeff(4)*mv2.coeff(17) - mv1.coeff(5)*mv2.coeff(16);
	}


	/// \brief Compute the outer product between two homogeneous multivectors mv1 (grade 1) and mv2 (grade 4). 
	/// \tparam the type of value that we manipulate, either float or double or something.
	/// \param mv1 - the first homogeneous multivector of grade 1 represented as an Eigen::VectorXd
	/// \param mv2 - the second homogeneous multivector of grade 4 represented as a Eigen::VectorXd
	/// \param mv3 - the result of mv1^mv2, which is also a homogeneous multivector of grade 5
	template<typename T>
	void outer_1_4(const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv1, const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv2, Eigen::Matrix<T, Eigen::Dynamic, 1>& mv3){
		mv3.coeffRef(0) +=  mv1.coeff(0)*mv2.coeff(10) - mv1.coeff(1)*mv2.coeff(6) + mv1.coeff(2)*mv2.coeff(3) - mv1.coeff(3)*mv2.coeff(1) + mv1.coeff(4)*mv2.coeff(0);
		mv3.coeffRef(1) +=  mv1.coeff(0)*mv2.coeff(11) - mv1.coeff(1)*mv2.coeff(7) + mv1.coeff(2)*mv2.coeff(4) - mv1.coeff(3)*mv2.coeff(2) + mv1.coeff(5)*mv2.coeff(0);
		mv3.coeffRef(2) +=  mv1.coeff(0)*mv2.coeff(12) - mv1.coeff(1)*mv2.coeff(8) + mv1.coeff(2)*mv2.coeff(5) - mv1.coeff(4)*mv2.coeff(2) + mv1.coeff(5)*mv2.coeff(1);
		mv3.coeffRef(3) +=  mv1.coeff(0)*mv2.coeff(13) - mv1.coeff(1)*mv2.coeff(9) + mv1.coeff(3)*mv2.coeff(5) - mv1.coeff(4)*mv2.coeff(4) + mv1.coeff(5)*mv2.coeff(3);
		mv3.coeffRef(4) +=  mv1.coeff(0)*mv2.coeff(14) - mv1.coeff(2)*mv2.coeff(9) + mv1.coeff(3)*mv2.coeff(8) - mv1.coeff(4)*mv2.coeff(7) + mv1.coeff(5)*mv2.coeff(6);
		mv3.coeffRef(5) +=  mv1.coeff(1)*mv2.coeff(14) - mv1.coeff(2)*mv2.coeff(13) + mv1.coeff(3)*mv2.coeff(12) - mv1.coeff(4)*mv2.coeff(11) + mv1.coeff(5)*mv2.coeff(10);
	}


	/// \brief Compute the outer product between two homogeneous multivectors mv1 (grade 1) and mv2 (grade 5). 
	/// \tparam the type of value that we manipulate, either float or double or something.
	/// \param mv1 - the first homogeneous multivector of grade 1 represented as an Eigen::VectorXd
	/// \param mv2 - the second homogeneous multivector of grade 5 represented as a Eigen::VectorXd
	/// \param mv3 - the result of mv1^mv2, which is also a homogeneous multivector of grade 6
	template<typename T>
	void outer_1_5(const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv1, const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv2, Eigen::Matrix<T, Eigen::Dynamic, 1>& mv3){
		mv3.coeffRef(0) +=  mv1.coeff(0)*mv2.coeff(5) - mv1.coeff(1)*mv2.coeff(4) + mv1.coeff(2)*mv2.coeff(3) - mv1.coeff(3)*mv2.coeff(2) + mv1.coeff(4)*mv2.coeff(1) - mv1.coeff(5)*mv2.coeff(0);
	}


	/// \brief Compute the outer product between two homogeneous multivectors mv1 (grade 2) and mv2 (grade 0). 
	/// \tparam the type of value that we manipulate, either float or double or something.
	/// \param mv1 - the first homogeneous multivector of grade 2 represented as an Eigen::VectorXd
	/// \param mv2 - the second homogeneous multivector of grade 0 represented as a Eigen::VectorXd
	/// \param mv3 - the result of mv1^mv2, which is also a homogeneous multivector of grade 2
	template<typename T>
	void outer_2_0(const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv1, const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv2, Eigen::Matrix<T, Eigen::Dynamic, 1>& mv3){
		mv3 += mv1*mv2.coeff(0);
	}


	/// \brief Compute the outer product between two homogeneous multivectors mv1 (grade 2) and mv2 (grade 1). 
	/// \tparam the type of value that we manipulate, either float or double or something.
	/// \param mv1 - the first homogeneous multivector of grade 2 represented as an Eigen::VectorXd
	/// \param mv2 - the second homogeneous multivector of grade 1 represented as a Eigen::VectorXd
	/// \param mv3 - the result of mv1^mv2, which is also a homogeneous multivector of grade 3
	template<typename T>
	void outer_2_1(const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv1, const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv2, Eigen::Matrix<T, Eigen::Dynamic, 1>& mv3){
		mv3.coeffRef(0) +=  mv1.coeff(0)*mv2.coeff(2) - mv1.coeff(1)*mv2.coeff(1) + mv1.coeff(5)*mv2.coeff(0);
		mv3.coeffRef(1) +=  mv1.coeff(0)*mv2.coeff(3) - mv1.coeff(2)*mv2.coeff(1) + mv1.coeff(6)*mv2.coeff(0);
		mv3.coeffRef(2) +=  mv1.coeff(0)*mv2.coeff(4) - mv1.coeff(3)*mv2.coeff(1) + mv1.coeff(7)*mv2.coeff(0);
		mv3.coeffRef(3) +=  mv1.coeff(0)*mv2.coeff(5) - mv1.coeff(4)*mv2.coeff(1) + mv1.coeff(8)*mv2.coeff(0);
		mv3.coeffRef(4) +=  mv1.coeff(1)*mv2.coeff(3) - mv1.coeff(2)*mv2.coeff(2) + mv1.coeff(9)*mv2.coeff(0);
		mv3.coeffRef(5) +=  mv1.coeff(1)*mv2.coeff(4) - mv1.coeff(3)*mv2.coeff(2) + mv1.coeff(10)*mv2.coeff(0);
		mv3.coeffRef(6) +=  mv1.coeff(1)*mv2.coeff(5) - mv1.coeff(4)*mv2.coeff(2) + mv1.coeff(11)*mv2.coeff(0);
		mv3.coeffRef(7) +=  mv1.coeff(2)*mv2.coeff(4) - mv1.coeff(3)*mv2.coeff(3) + mv1.coeff(12)*mv2.coeff(0);
		mv3.coeffRef(8) +=  mv1.coeff(2)*mv2.coeff(5) - mv1.coeff(4)*mv2.coeff(3) + mv1.coeff(13)*mv2.coeff(0);
		mv3.coeffRef(9) +=  mv1.coeff(3)*mv2.coeff(5) - mv1.coeff(4)*mv2.coeff(4) + mv1.coeff(14)*mv2.coeff(0);
		mv3.coeffRef(10) +=  mv1.coeff(5)*mv2.coeff(3) - mv1.coeff(6)*mv2.coeff(2) + mv1.coeff(9)*mv2.coeff(1);
		mv3.coeffRef(11) +=  mv1.coeff(5)*mv2.coeff(4) - mv1.coeff(7)*mv2.coeff(2) + mv1.coeff(10)*mv2.coeff(1);
		mv3.coeffRef(12) +=  mv1.coeff(5)*mv2.coeff(5) - mv1.coeff(8)*mv2.coeff(2) + mv1.coeff(11)*mv2.coeff(1);
		mv3.coeffRef(13) +=  mv1.coeff(6)*mv2.coeff(4) - mv1.coeff(7)*mv2.coeff(3) + mv1.coeff(12)*mv2.coeff(1);
		mv3.coeffRef(14) +=  mv1.coeff(6)*mv2.coeff(5) - mv1.coeff(8)*mv2.coeff(3) + mv1.coeff(13)*mv2.coeff(1);
		mv3.coeffRef(15) +=  mv1.coeff(7)*mv2.coeff(5) - mv1.coeff(8)*mv2.coeff(4) + mv1.coeff(14)*mv2.coeff(1);
		mv3.coeffRef(16) +=  mv1.coeff(9)*mv2.coeff(4) - mv1.coeff(10)*mv2.coeff(3) + mv1.coeff(12)*mv2.coeff(2);
		mv3.coeffRef(17) +=  mv1.coeff(9)*mv2.coeff(5) - mv1.coeff(11)*mv2.coeff(3) + mv1.coeff(13)*mv2.coeff(2);
		mv3.coeffRef(18) +=  mv1.coeff(10)*mv2.coeff(5) - mv1.coeff(11)*mv2.coeff(4) + mv1.coeff(14)*mv2.coeff(2);
		mv3.coeffRef(19) +=  mv1.coeff(12)*mv2.coeff(5) - mv1.coeff(13)*mv2.coeff(4) + mv1.coeff(14)*mv2.coeff(3);
	}


	/// \brief Compute the outer product between two homogeneous multivectors mv1 (grade 2) and mv2 (grade 2). 
	/// \tparam the type of value that we manipulate, either float or double or something.
	/// \param mv1 - the first homogeneous multivector of grade 2 represented as an Eigen::VectorXd
	/// \param mv2 - the second homogeneous multivector of grade 2 represented as a Eigen::VectorXd
	/// \param mv3 - the result of mv1^mv2, which is also a homogeneous multivector of grade 4
	template<typename T>
	void outer_2_2(const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv1, const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv2, Eigen::Matrix<T, Eigen::Dynamic, 1>& mv3){
		mv3.coeffRef(0) +=  mv1.coeff(0)*mv2.coeff(9) - mv1.coeff(1)*mv2.coeff(6) + mv1.coeff(2)*mv2.coeff(5) + mv1.coeff(5)*mv2.coeff(2) - mv1.coeff(6)*mv2.coeff(1) + mv1.coeff(9)*mv2.coeff(0);
		mv3.coeffRef(1) +=  mv1.coeff(0)*mv2.coeff(10) - mv1.coeff(1)*mv2.coeff(7) + mv1.coeff(3)*mv2.coeff(5) + mv1.coeff(5)*mv2.coeff(3) - mv1.coeff(7)*mv2.coeff(1) + mv1.coeff(10)*mv2.coeff(0);
		mv3.coeffRef(2) +=  mv1.coeff(0)*mv2.coeff(11) - mv1.coeff(1)*mv2.coeff(8) + mv1.coeff(4)*mv2.coeff(5) + mv1.coeff(5)*mv2.coeff(4) - mv1.coeff(8)*mv2.coeff(1) + mv1.coeff(11)*mv2.coeff(0);
		mv3.coeffRef(3) +=  mv1.coeff(0)*mv2.coeff(12) - mv1.coeff(2)*mv2.coeff(7) + mv1.coeff(3)*mv2.coeff(6) + mv1.coeff(6)*mv2.coeff(3) - mv1.coeff(7)*mv2.coeff(2) + mv1.coeff(12)*mv2.coeff(0);
		mv3.coeffRef(4) +=  mv1.coeff(0)*mv2.coeff(13) - mv1.coeff(2)*mv2.coeff(8) + mv1.coeff(4)*mv2.coeff(6) + mv1.coeff(6)*mv2.coeff(4) - mv1.coeff(8)*mv2.coeff(2) + mv1.coeff(13)*mv2.coeff(0);
		mv3.coeffRef(5) +=  mv1.coeff(0)*mv2.coeff(14) - mv1.coeff(3)*mv2.coeff(8) + mv1.coeff(4)*mv2.coeff(7) + mv1.coeff(7)*mv2.coeff(4) - mv1.coeff(8)*mv2.coeff(3) + mv1.coeff(14)*mv2.coeff(0);
		mv3.coeffRef(6) +=  mv1.coeff(1)*mv2.coeff(12) - mv1.coeff(2)*mv2.coeff(10) + mv1.coeff(3)*mv2.coeff(9) + mv1.coeff(9)*mv2.coeff(3) - mv1.coeff(10)*mv2.coeff(2) + mv1.coeff(12)*mv2.coeff(1);
		mv3.coeffRef(7) +=  mv1.coeff(1)*mv2.coeff(13) - mv1.coeff(2)*mv2.coeff(11) + mv1.coeff(4)*mv2.coeff(9) + mv1.coeff(9)*mv2.coeff(4) - mv1.coeff(11)*mv2.coeff(2) + mv1.coeff(13)*mv2.coeff(1);
		mv3.coeffRef(8) +=  mv1.coeff(1)*mv2.coeff(14) - mv1.coeff(3)*mv2.coeff(11) + mv1.coeff(4)*mv2.coeff(10) + mv1.coeff(10)*mv2.coeff(4) - mv1.coeff(11)*mv2.coeff(3) + mv1.coeff(14)*mv2.coeff(1);
		mv3.coeffRef(9) +=  mv1.coeff(2)*mv2.coeff(14) - mv1.coeff(3)*mv2.coeff(13) + mv1.coeff(4)*mv2.coeff(12) + mv1.coeff(12)*mv2.coeff(4) - mv1.coeff(13)*mv2.coeff(3) + mv1.coeff(14)*mv2.coeff(2);
		mv3.coeffRef(10) +=  mv1.coeff(5)*mv2.coeff(12) - mv1.coeff(6)*mv2.coeff(10) + mv1.coeff(7)*mv2.coeff(9) + mv1.coeff(9)*mv2.coeff(7) - mv1.coeff(10)*mv2.coeff(6) + mv1.coeff(12)*mv2.coeff(5);
		mv3.coeffRef(11) +=  mv1.coeff(5)*mv2.coeff(13) - mv1.coeff(6)*mv2.coeff(11) + mv1.coeff(8)*mv2.coeff(9) + mv1.coeff(9)*mv2.coeff(8) - mv1.coeff(11)*mv2.coeff(6) + mv1.coeff(13)*mv2.coeff(5);
		mv3.coeffRef(12) +=  mv1.coeff(5)*mv2.coeff(14) - mv1.coeff(7)*mv2.coeff(11) + mv1.coeff(8)*mv2.coeff(10) + mv1.coeff(10)*mv2.coeff(8) - mv1.coeff(11)*mv2.coeff(7) + mv1.coeff(14)*mv2.coeff(5);
		mv3.coeffRef(13) +=  mv1.coeff(6)*mv2.coeff(14) - mv1.coeff(7)*mv2.coeff(13) + mv1.coeff(8)*mv2.coeff(12) + mv1.coeff(12)*mv2.coeff(8) - mv1.coeff(13)*mv2.coeff(7) + mv1.coeff(14)*mv2.coeff(6);
		mv3.coeffRef(14) +=  mv1.coeff(9)*mv2.coeff(14) - mv1.coeff(10)*mv2.coeff(13) + mv1.coeff(11)*mv2.coeff(12) + mv1.coeff(12)*mv2.coeff(11) - mv1.coeff(13)*mv2.coeff(10) + mv1.coeff(14)*mv2.coeff(9);
	}


	/// \brief Compute the outer product between two homogeneous multivectors mv1 (grade 2) and mv2 (grade 3). 
	/// \tparam the type of value that we manipulate, either float or double or something.
	/// \param mv1 - the first homogeneous multivector of grade 2 represented as an Eigen::VectorXd
	/// \param mv2 - the second homogeneous multivector of grade 3 represented as a Eigen::VectorXd
	/// \param mv3 - the result of mv1^mv2, which is also a homogeneous multivector of grade 5
	template<typename T>
	void outer_2_3(const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv1, const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv2, Eigen::Matrix<T, Eigen::Dynamic, 1>& mv3){
		mv3.coeffRef(0) +=  mv1.coeff(0)*mv2.coeff(16) - mv1.coeff(1)*mv2.coeff(13) + mv1.coeff(2)*mv2.coeff(11) - mv1.coeff(3)*mv2.coeff(10) + mv1.coeff(5)*mv2.coeff(7) - mv1.coeff(6)*mv2.coeff(5) + mv1.coeff(7)*mv2.coeff(4) + mv1.coeff(9)*mv2.coeff(2) - mv1.coeff(10)*mv2.coeff(1) + mv1.coeff(12)*mv2.coeff(0);
		mv3.coeffRef(1) +=  mv1.coeff(0)*mv2.coeff(17) - mv1.coeff(1)*mv2.coeff(14) + mv1.coeff(2)*mv2.coeff(12) - mv1.coeff(4)*mv2.coeff(10) + mv1.coeff(5)*mv2.coeff(8) - mv1.coeff(6)*mv2.coeff(6) + mv1.coeff(8)*mv2.coeff(4) + mv1.coeff(9)*mv2.coeff(3) - mv1.coeff(11)*mv2.coeff(1) + mv1.coeff(13)*mv2.coeff(0);
		mv3.coeffRef(2) +=  mv1.coeff(0)*mv2.coeff(18) - mv1.coeff(1)*mv2.coeff(15) + mv1.coeff(3)*mv2.coeff(12) - mv1.coeff(4)*mv2.coeff(11) + mv1.coeff(5)*mv2.coeff(9) - mv1.coeff(7)*mv2.coeff(6) + mv1.coeff(8)*mv2.coeff(5) + mv1.coeff(10)*mv2.coeff(3) - mv1.coeff(11)*mv2.coeff(2) + mv1.coeff(14)*mv2.coeff(0);
		mv3.coeffRef(3) +=  mv1.coeff(0)*mv2.coeff(19) - mv1.coeff(2)*mv2.coeff(15) + mv1.coeff(3)*mv2.coeff(14) - mv1.coeff(4)*mv2.coeff(13) + mv1.coeff(6)*mv2.coeff(9) - mv1.coeff(7)*mv2.coeff(8) + mv1.coeff(8)*mv2.coeff(7) + mv1.coeff(12)*mv2.coeff(3) - mv1.coeff(13)*mv2.coeff(2) + mv1.coeff(14)*mv2.coeff(1);
		mv3.coeffRef(4) +=  mv1.coeff(1)*mv2.coeff(19) - mv1.coeff(2)*mv2.coeff(18) + mv1.coeff(3)*mv2.coeff(17) - mv1.coeff(4)*mv2.coeff(16) + mv1.coeff(9)*mv2.coeff(9) - mv1.coeff(10)*mv2.coeff(8) + mv1.coeff(11)*mv2.coeff(7) + mv1.coeff(12)*mv2.coeff(6) - mv1.coeff(13)*mv2.coeff(5) + mv1.coeff(14)*mv2.coeff(4);
		mv3.coeffRef(5) +=  mv1.coeff(5)*mv2.coeff(19) - mv1.coeff(6)*mv2.coeff(18) + mv1.coeff(7)*mv2.coeff(17) - mv1.coeff(8)*mv2.coeff(16) + mv1.coeff(9)*mv2.coeff(15) - mv1.coeff(10)*mv2.coeff(14) + mv1.coeff(11)*mv2.coeff(13) + mv1.coeff(12)*mv2.coeff(12) - mv1.coeff(13)*mv2.coeff(11) + mv1.coeff(14)*mv2.coeff(10);
	}


	/// \brief Compute the outer product between two homogeneous multivectors mv1 (grade 2) and mv2 (grade 4). 
	/// \tparam the type of value that we manipulate, either float or double or something.
	/// \param mv1 - the first homogeneous multivector of grade 2 represented as an Eigen::VectorXd
	/// \param mv2 - the second homogeneous multivector of grade 4 represented as a Eigen::VectorXd
	/// \param mv3 - the result of mv1^mv2, which is also a homogeneous multivector of grade 6
	template<typename T>
	void outer_2_4(const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv1, const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv2, Eigen::Matrix<T, Eigen::Dynamic, 1>& mv3){
		mv3.coeffRef(0) +=  mv1.coeff(0)*mv2.coeff(14) - mv1.coeff(1)*mv2.coeff(13) + mv1.coeff(2)*mv2.coeff(12) - mv1.coeff(3)*mv2.coeff(11) + mv1.coeff(4)*mv2.coeff(10) + mv1.coeff(5)*mv2.coeff(9) - mv1.coeff(6)*mv2.coeff(8) + mv1.coeff(7)*mv2.coeff(7) - mv1.coeff(8)*mv2.coeff(6) + mv1.coeff(9)*mv2.coeff(5) - mv1.coeff(10)*mv2.coeff(4) + mv1.coeff(11)*mv2.coeff(3) + mv1.coeff(12)*mv2.coeff(2) - mv1.coeff(13)*mv2.coeff(1) + mv1.coeff(14)*mv2.coeff(0);
	}


	/// \brief Compute the outer product between two homogeneous multivectors mv1 (grade 3) and mv2 (grade 0). 
	/// \tparam the type of value that we manipulate, either float or double or something.
	/// \param mv1 - the first homogeneous multivector of grade 3 represented as an Eigen::VectorXd
	/// \param mv2 - the second homogeneous multivector of grade 0 represented as a Eigen::VectorXd
	/// \param mv3 - the result of mv1^mv2, which is also a homogeneous multivector of grade 3
	template<typename T>
	void outer_3_0(const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv1, const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv2, Eigen::Matrix<T, Eigen::Dynamic, 1>& mv3){
		mv3 += mv1*mv2.coeff(0);
	}


	/// \brief Compute the outer product between two homogeneous multivectors mv1 (grade 3) and mv2 (grade 1). 
	/// \tparam the type of value that we manipulate, either float or double or something.
	/// \param mv1 - the first homogeneous multivector of grade 3 represented as an Eigen::VectorXd
	/// \param mv2 - the second homogeneous multivector of grade 1 represented as a Eigen::VectorXd
	/// \param mv3 - the result of mv1^mv2, which is also a homogeneous multivector of grade 4
	template<typename T>
	void outer_3_1(const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv1, const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv2, Eigen::Matrix<T, Eigen::Dynamic, 1>& mv3){
		mv3.coeffRef(0) +=  mv1.coeff(0)*mv2.coeff(3) - mv1.coeff(1)*mv2.coeff(2) + mv1.coeff(4)*mv2.coeff(1) - mv1.coeff(10)*mv2.coeff(0);
		mv3.coeffRef(1) +=  mv1.coeff(0)*mv2.coeff(4) - mv1.coeff(2)*mv2.coeff(2) + mv1.coeff(5)*mv2.coeff(1) - mv1.coeff(11)*mv2.coeff(0);
		mv3.coeffRef(2) +=  mv1.coeff(0)*mv2.coeff(5) - mv1.coeff(3)*mv2.coeff(2) + mv1.coeff(6)*mv2.coeff(1) - mv1.coeff(12)*mv2.coeff(0);
		mv3.coeffRef(3) +=  mv1.coeff(1)*mv2.coeff(4) - mv1.coeff(2)*mv2.coeff(3) + mv1.coeff(7)*mv2.coeff(1) - mv1.coeff(13)*mv2.coeff(0);
		mv3.coeffRef(4) +=  mv1.coeff(1)*mv2.coeff(5) - mv1.coeff(3)*mv2.coeff(3) + mv1.coeff(8)*mv2.coeff(1) - mv1.coeff(14)*mv2.coeff(0);
		mv3.coeffRef(5) +=  mv1.coeff(2)*mv2.coeff(5) - mv1.coeff(3)*mv2.coeff(4) + mv1.coeff(9)*mv2.coeff(1) - mv1.coeff(15)*mv2.coeff(0);
		mv3.coeffRef(6) +=  mv1.coeff(4)*mv2.coeff(4) - mv1.coeff(5)*mv2.coeff(3) + mv1.coeff(7)*mv2.coeff(2) - mv1.coeff(16)*mv2.coeff(0);
		mv3.coeffRef(7) +=  mv1.coeff(4)*mv2.coeff(5) - mv1.coeff(6)*mv2.coeff(3) + mv1.coeff(8)*mv2.coeff(2) - mv1.coeff(17)*mv2.coeff(0);
		mv3.coeffRef(8) +=  mv1.coeff(5)*mv2.coeff(5) - mv1.coeff(6)*mv2.coeff(4) + mv1.coeff(9)*mv2.coeff(2) - mv1.coeff(18)*mv2.coeff(0);
		mv3.coeffRef(9) +=  mv1.coeff(7)*mv2.coeff(5) - mv1.coeff(8)*mv2.coeff(4) + mv1.coeff(9)*mv2.coeff(3) - mv1.coeff(19)*mv2.coeff(0);
		mv3.coeffRef(10) +=  mv1.coeff(10)*mv2.coeff(4) - mv1.coeff(11)*mv2.coeff(3) + mv1.coeff(13)*mv2.coeff(2) - mv1.coeff(16)*mv2.coeff(1);
		mv3.coeffRef(11) +=  mv1.coeff(10)*mv2.coeff(5) - mv1.coeff(12)*mv2.coeff(3) + mv1.coeff(14)*mv2.coeff(2) - mv1.coeff(17)*mv2.coeff(1);
		mv3.coeffRef(12) +=  mv1.coeff(11)*mv2.coeff(5) - mv1.coeff(12)*mv2.coeff(4) + mv1.coeff(15)*mv2.coeff(2) - mv1.coeff(18)*mv2.coeff(1);
		mv3.coeffRef(13) +=  mv1.coeff(13)*mv2.coeff(5) - mv1.coeff(14)*mv2.coeff(4) + mv1.coeff(15)*mv2.coeff(3) - mv1.coeff(19)*mv2.coeff(1);
		mv3.coeffRef(14) +=  mv1.coeff(16)*mv2.coeff(5) - mv1.coeff(17)*mv2.coeff(4) + mv1.coeff(18)*mv2.coeff(3) - mv1.coeff(19)*mv2.coeff(2);
	}


	/// \brief Compute the outer product between two homogeneous multivectors mv1 (grade 3) and mv2 (grade 2). 
	/// \tparam the type of value that we manipulate, either float or double or something.
	/// \param mv1 - the first homogeneous multivector of grade 3 represented as an Eigen::VectorXd
	/// \param mv2 - the second homogeneous multivector of grade 2 represented as a Eigen::VectorXd
	/// \param mv3 - the result of mv1^mv2, which is also a homogeneous multivector of grade 5
	template<typename T>
	void outer_3_2(const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv1, const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv2, Eigen::Matrix<T, Eigen::Dynamic, 1>& mv3){
		mv3.coeffRef(0) +=  mv1.coeff(0)*mv2.coeff(12) - mv1.coeff(1)*mv2.coeff(10) + mv1.coeff(2)*mv2.coeff(9) + mv1.coeff(4)*mv2.coeff(7) - mv1.coeff(5)*mv2.coeff(6) + mv1.coeff(7)*mv2.coeff(5) - mv1.coeff(10)*mv2.coeff(3) + mv1.coeff(11)*mv2.coeff(2) - mv1.coeff(13)*mv2.coeff(1) + mv1.coeff(16)*mv2.coeff(0);
		mv3.coeffRef(1) +=  mv1.coeff(0)*mv2.coeff(13) - mv1.coeff(1)*mv2.coeff(11) + mv1.coeff(3)*mv2.coeff(9) + mv1.coeff(4)*mv2.coeff(8) - mv1.coeff(6)*mv2.coeff(6) + mv1.coeff(8)*mv2.coeff(5) - mv1.coeff(10)*mv2.coeff(4) + mv1.coeff(12)*mv2.coeff(2) - mv1.coeff(14)*mv2.coeff(1) + mv1.coeff(17)*mv2.coeff(0);
		mv3.coeffRef(2) +=  mv1.coeff(0)*mv2.coeff(14) - mv1.coeff(2)*mv2.coeff(11) + mv1.coeff(3)*mv2.coeff(10) + mv1.coeff(5)*mv2.coeff(8) - mv1.coeff(6)*mv2.coeff(7) + mv1.coeff(9)*mv2.coeff(5) - mv1.coeff(11)*mv2.coeff(4) + mv1.coeff(12)*mv2.coeff(3) - mv1.coeff(15)*mv2.coeff(1) + mv1.coeff(18)*mv2.coeff(0);
		mv3.coeffRef(3) +=  mv1.coeff(1)*mv2.coeff(14) - mv1.coeff(2)*mv2.coeff(13) + mv1.coeff(3)*mv2.coeff(12) + mv1.coeff(7)*mv2.coeff(8) - mv1.coeff(8)*mv2.coeff(7) + mv1.coeff(9)*mv2.coeff(6) - mv1.coeff(13)*mv2.coeff(4) + mv1.coeff(14)*mv2.coeff(3) - mv1.coeff(15)*mv2.coeff(2) + mv1.coeff(19)*mv2.coeff(0);
		mv3.coeffRef(4) +=  mv1.coeff(4)*mv2.coeff(14) - mv1.coeff(5)*mv2.coeff(13) + mv1.coeff(6)*mv2.coeff(12) + mv1.coeff(7)*mv2.coeff(11) - mv1.coeff(8)*mv2.coeff(10) + mv1.coeff(9)*mv2.coeff(9) - mv1.coeff(16)*mv2.coeff(4) + mv1.coeff(17)*mv2.coeff(3) - mv1.coeff(18)*mv2.coeff(2) + mv1.coeff(19)*mv2.coeff(1);
		mv3.coeffRef(5) +=  mv1.coeff(10)*mv2.coeff(14) - mv1.coeff(11)*mv2.coeff(13) + mv1.coeff(12)*mv2.coeff(12) + mv1.coeff(13)*mv2.coeff(11) - mv1.coeff(14)*mv2.coeff(10) + mv1.coeff(15)*mv2.coeff(9) - mv1.coeff(16)*mv2.coeff(8) + mv1.coeff(17)*mv2.coeff(7) - mv1.coeff(18)*mv2.coeff(6) + mv1.coeff(19)*mv2.coeff(5);
	}


	/// \brief Compute the outer product between two homogeneous multivectors mv1 (grade 3) and mv2 (grade 3). 
	/// \tparam the type of value that we manipulate, either float or double or something.
	/// \param mv1 - the first homogeneous multivector of grade 3 represented as an Eigen::VectorXd
	/// \param mv2 - the second homogeneous multivector of grade 3 represented as a Eigen::VectorXd
	/// \param mv3 - the result of mv1^mv2, which is also a homogeneous multivector of grade 6
	template<typename T>
	void outer_3_3(const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv1, const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv2, Eigen::Matrix<T, Eigen::Dynamic, 1>& mv3){
		mv3.coeffRef(0) +=  mv1.coeff(0)*mv2.coeff(19) - mv1.coeff(1)*mv2.coeff(18) + mv1.coeff(2)*mv2.coeff(17) - mv1.coeff(3)*mv2.coeff(16) + mv1.coeff(4)*mv2.coeff(15) - mv1.coeff(5)*mv2.coeff(14) + mv1.coeff(6)*mv2.coeff(13) + mv1.coeff(7)*mv2.coeff(12) - mv1.coeff(8)*mv2.coeff(11) + mv1.coeff(9)*mv2.coeff(10) - mv1.coeff(10)*mv2.coeff(9) + mv1.coeff(11)*mv2.coeff(8) - mv1.coeff(12)*mv2.coeff(7) - mv1.coeff(13)*mv2.coeff(6) + mv1.coeff(14)*mv2.coeff(5) - mv1.coeff(15)*mv2.coeff(4) + mv1.coeff(16)*mv2.coeff(3) - mv1.coeff(17)*mv2.coeff(2) + mv1.coeff(18)*mv2.coeff(1) - mv1.coeff(19)*mv2.coeff(0);
	}


	/// \brief Compute the outer product between two homogeneous multivectors mv1 (grade 4) and mv2 (grade 0). 
	/// \tparam the type of value that we manipulate, either float or double or something.
	/// \param mv1 - the first homogeneous multivector of grade 4 represented as an Eigen::VectorXd
	/// \param mv2 - the second homogeneous multivector of grade 0 represented as a Eigen::VectorXd
	/// \param mv3 - the result of mv1^mv2, which is also a homogeneous multivector of grade 4
	template<typename T>
	void outer_4_0(const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv1, const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv2, Eigen::Matrix<T, Eigen::Dynamic, 1>& mv3){
		mv3 += mv1*mv2.coeff(0);
	}


	/// \brief Compute the outer product between two homogeneous multivectors mv1 (grade 4) and mv2 (grade 1). 
	/// \tparam the type of value that we manipulate, either float or double or something.
	/// \param mv1 - the first homogeneous multivector of grade 4 represented as an Eigen::VectorXd
	/// \param mv2 - the second homogeneous multivector of grade 1 represented as a Eigen::VectorXd
	/// \param mv3 - the result of mv1^mv2, which is also a homogeneous multivector of grade 5
	template<typename T>
	void outer_4_1(const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv1, const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv2, Eigen::Matrix<T, Eigen::Dynamic, 1>& mv3){
		mv3.coeffRef(0) +=  mv1.coeff(0)*mv2.coeff(4) - mv1.coeff(1)*mv2.coeff(3) + mv1.coeff(3)*mv2.coeff(2) - mv1.coeff(6)*mv2.coeff(1) + mv1.coeff(10)*mv2.coeff(0);
		mv3.coeffRef(1) +=  mv1.coeff(0)*mv2.coeff(5) - mv1.coeff(2)*mv2.coeff(3) + mv1.coeff(4)*mv2.coeff(2) - mv1.coeff(7)*mv2.coeff(1) + mv1.coeff(11)*mv2.coeff(0);
		mv3.coeffRef(2) +=  mv1.coeff(1)*mv2.coeff(5) - mv1.coeff(2)*mv2.coeff(4) + mv1.coeff(5)*mv2.coeff(2) - mv1.coeff(8)*mv2.coeff(1) + mv1.coeff(12)*mv2.coeff(0);
		mv3.coeffRef(3) +=  mv1.coeff(3)*mv2.coeff(5) - mv1.coeff(4)*mv2.coeff(4) + mv1.coeff(5)*mv2.coeff(3) - mv1.coeff(9)*mv2.coeff(1) + mv1.coeff(13)*mv2.coeff(0);
		mv3.coeffRef(4) +=  mv1.coeff(6)*mv2.coeff(5) - mv1.coeff(7)*mv2.coeff(4) + mv1.coeff(8)*mv2.coeff(3) - mv1.coeff(9)*mv2.coeff(2) + mv1.coeff(14)*mv2.coeff(0);
		mv3.coeffRef(5) +=  mv1.coeff(10)*mv2.coeff(5) - mv1.coeff(11)*mv2.coeff(4) + mv1.coeff(12)*mv2.coeff(3) - mv1.coeff(13)*mv2.coeff(2) + mv1.coeff(14)*mv2.coeff(1);
	}


	/// \brief Compute the outer product between two homogeneous multivectors mv1 (grade 4) and mv2 (grade 2). 
	/// \tparam the type of value that we manipulate, either float or double or something.
	/// \param mv1 - the first homogeneous multivector of grade 4 represented as an Eigen::VectorXd
	/// \param mv2 - the second homogeneous multivector of grade 2 represented as a Eigen::VectorXd
	/// \param mv3 - the result of mv1^mv2, which is also a homogeneous multivector of grade 6
	template<typename T>
	void outer_4_2(const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv1, const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv2, Eigen::Matrix<T, Eigen::Dynamic, 1>& mv3){
		mv3.coeffRef(0) +=  mv1.coeff(0)*mv2.coeff(14) - mv1.coeff(1)*mv2.coeff(13) + mv1.coeff(2)*mv2.coeff(12) + mv1.coeff(3)*mv2.coeff(11) - mv1.coeff(4)*mv2.coeff(10) + mv1.coeff(5)*mv2.coeff(9) - mv1.coeff(6)*mv2.coeff(8) + mv1.coeff(7)*mv2.coeff(7) - mv1.coeff(8)*mv2.coeff(6) + mv1.coeff(9)*mv2.coeff(5) + mv1.coeff(10)*mv2.coeff(4) - mv1.coeff(11)*mv2.coeff(3) + mv1.coeff(12)*mv2.coeff(2) - mv1.coeff(13)*mv2.coeff(1) + mv1.coeff(14)*mv2.coeff(0);
	}


	/// \brief Compute the outer product between two homogeneous multivectors mv1 (grade 5) and mv2 (grade 0). 
	/// \tparam the type of value that we manipulate, either float or double or something.
	/// \param mv1 - the first homogeneous multivector of grade 5 represented as an Eigen::VectorXd
	/// \param mv2 - the second homogeneous multivector of grade 0 represented as a Eigen::VectorXd
	/// \param mv3 - the result of mv1^mv2, which is also a homogeneous multivector of grade 5
	template<typename T>
	void outer_5_0(const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv1, const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv2, Eigen::Matrix<T, Eigen::Dynamic, 1>& mv3){
		mv3 += mv1*mv2.coeff(0);
	}


	/// \brief Compute the outer product between two homogeneous multivectors mv1 (grade 5) and mv2 (grade 1). 
	/// \tparam the type of value that we manipulate, either float or double or something.
	/// \param mv1 - the first homogeneous multivector of grade 5 represented as an Eigen::VectorXd
	/// \param mv2 - the second homogeneous multivector of grade 1 represented as a Eigen::VectorXd
	/// \param mv3 - the result of mv1^mv2, which is also a homogeneous multivector of grade 6
	template<typename T>
	void outer_5_1(const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv1, const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv2, Eigen::Matrix<T, Eigen::Dynamic, 1>& mv3){
		mv3.coeffRef(0) +=  mv1.coeff(0)*mv2.coeff(5) - mv1.coeff(1)*mv2.coeff(4) + mv1.coeff(2)*mv2.coeff(3) - mv1.coeff(3)*mv2.coeff(2) + mv1.coeff(4)*mv2.coeff(1) - mv1.coeff(5)*mv2.coeff(0);
	}


	/// \brief Compute the outer product between two homogeneous multivectors mv1 (grade 6) and mv2 (grade 0). 
	/// \tparam the type of value that we manipulate, either float or double or something.
	/// \param mv1 - the first homogeneous multivector of grade 6 represented as an Eigen::VectorXd
	/// \param mv2 - the second homogeneous multivector of grade 0 represented as a Eigen::VectorXd
	/// \param mv3 - the result of mv1^mv2, which is also a homogeneous multivector of grade 6
	template<typename T>
	void outer_6_0(const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv1, const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv2, Eigen::Matrix<T, Eigen::Dynamic, 1>& mv3){
		mv3 += mv1*mv2.coeff(0);
	}


	

    template<typename T>
	std::array<std::array<std::function<void(const Eigen::Matrix<T, Eigen::Dynamic, 1> & , const Eigen::Matrix<T, Eigen::Dynamic, 1> & , Eigen::Matrix<T, Eigen::Dynamic, 1>&)>, 7>, 7> outerFunctionsContainer = {{
		{{outer_0_0<T>,outer_0_1<T>,outer_0_2<T>,outer_0_3<T>,outer_0_4<T>,outer_0_5<T>,outer_0_6<T>}},
		{{outer_1_0<T>,outer_1_1<T>,outer_1_2<T>,outer_1_3<T>,outer_1_4<T>,outer_1_5<T>,{}}},
		{{outer_2_0<T>,outer_2_1<T>,outer_2_2<T>,outer_2_3<T>,outer_2_4<T>,{},{}}},
		{{outer_3_0<T>,outer_3_1<T>,outer_3_2<T>,outer_3_3<T>,{},{},{}}},
		{{outer_4_0<T>,outer_4_1<T>,outer_4_2<T>,{},{},{},{}}},
		{{outer_5_0<T>,outer_5_1<T>,{},{},{},{},{}}},
		{{outer_6_0<T>,{},{},{},{},{},{}}}
	}};

}/// End of Namespace

#endif // C4GA_OUTER_PRODUCT_EXPLICIT_HPP__