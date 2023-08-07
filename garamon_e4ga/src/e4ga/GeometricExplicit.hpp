// Copyright (c) 2018 by University Paris-Est Marne-la-Vallee
// GeometricExplicit.hpp
// This file is part of the Garamon for e4ga.
// Authors: Stephane Breuils and Vincent Nozick
// Conctat: vincent.nozick@u-pem.fr
//
// Licence MIT
// A a copy of the MIT License is given along with this program

/// \file GeometricExplicit.hpp
/// \author Stephane Breuils, Vincent Nozick
/// \brief Explicit precomputed per grades geometric products of e4ga.


#ifndef E4GA_GEOMETRIC_PRODUCT_EXPLICIT_HPP__
#define E4GA_GEOMETRIC_PRODUCT_EXPLICIT_HPP__
#pragma once

#include <Eigen/Core>

#include "e4ga/Mvec.hpp"
#include "e4ga/Constants.hpp"


/*!
 * @namespace e4ga
 */
namespace e4ga {
    template<typename T> class Mvec;

    /// \brief Compute the geometric product between two homogeneous multivectors mv1 (grade 2) and mv2 (grade 2). 
	/// \tparam the type of value that we manipulate, either float or double or something else.
	/// \param mv1 - the first homogeneous multivector of grade 2 represented as an Eigen::VectorXd
	/// \param mv2 - the second homogeneous multivector of grade 2 represented as a Eigen::VectorXd
	/// \param mv3 - the result of mv1 mv2 whose grade is 2
	template<typename T>
	void geometric_2_2_2(const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv1, const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv2, Eigen::Matrix<T, Eigen::Dynamic, 1>& mv3){
		mv3.coeffRef(0) += -mv1.coeff(1)*mv2.coeff(3) - mv1.coeff(2)*mv2.coeff(4) + mv1.coeff(3)*mv2.coeff(1) + mv1.coeff(4)*mv2.coeff(2);
		mv3.coeffRef(1) +=  mv1.coeff(0)*mv2.coeff(3) - mv1.coeff(2)*mv2.coeff(5) - mv1.coeff(3)*mv2.coeff(0) + mv1.coeff(5)*mv2.coeff(2);
		mv3.coeffRef(2) +=  mv1.coeff(0)*mv2.coeff(4) + mv1.coeff(1)*mv2.coeff(5) - mv1.coeff(4)*mv2.coeff(0) - mv1.coeff(5)*mv2.coeff(1);
		mv3.coeffRef(3) += -mv1.coeff(0)*mv2.coeff(1) + mv1.coeff(1)*mv2.coeff(0) - mv1.coeff(4)*mv2.coeff(5) + mv1.coeff(5)*mv2.coeff(4);
		mv3.coeffRef(4) += -mv1.coeff(0)*mv2.coeff(2) + mv1.coeff(2)*mv2.coeff(0) + mv1.coeff(3)*mv2.coeff(5) - mv1.coeff(5)*mv2.coeff(3);
		mv3.coeffRef(5) += -mv1.coeff(1)*mv2.coeff(2) + mv1.coeff(2)*mv2.coeff(1) - mv1.coeff(3)*mv2.coeff(4) + mv1.coeff(4)*mv2.coeff(3);
	}


	/// \brief Compute the geometric product between two homogeneous multivectors mv1 (grade 2) and mv2 (grade 3). 
	/// \tparam the type of value that we manipulate, either float or double or something else.
	/// \param mv1 - the first homogeneous multivector of grade 2 represented as an Eigen::VectorXd
	/// \param mv2 - the second homogeneous multivector of grade 3 represented as a Eigen::VectorXd
	/// \param mv3 - the result of mv1 mv2 whose grade is 3
	template<typename T>
	void geometric_2_3_3(const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv1, const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv2, Eigen::Matrix<T, Eigen::Dynamic, 1>& mv3){
		mv3.coeffRef(0) +=  mv1.coeff(2)*mv2.coeff(3) - mv1.coeff(4)*mv2.coeff(2) + mv1.coeff(5)*mv2.coeff(1);
		mv3.coeffRef(1) += -mv1.coeff(1)*mv2.coeff(3) + mv1.coeff(3)*mv2.coeff(2) - mv1.coeff(5)*mv2.coeff(0);
		mv3.coeffRef(2) +=  mv1.coeff(0)*mv2.coeff(3) - mv1.coeff(3)*mv2.coeff(1) + mv1.coeff(4)*mv2.coeff(0);
		mv3.coeffRef(3) += -mv1.coeff(0)*mv2.coeff(2) + mv1.coeff(1)*mv2.coeff(1) - mv1.coeff(2)*mv2.coeff(0);
	}


	/// \brief Compute the geometric product between two homogeneous multivectors mv1 (grade 3) and mv2 (grade 2). 
	/// \tparam the type of value that we manipulate, either float or double or something else.
	/// \param mv1 - the first homogeneous multivector of grade 3 represented as an Eigen::VectorXd
	/// \param mv2 - the second homogeneous multivector of grade 2 represented as a Eigen::VectorXd
	/// \param mv3 - the result of mv1 mv2 whose grade is 3
	template<typename T>
	void geometric_3_2_3(const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv1, const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv2, Eigen::Matrix<T, Eigen::Dynamic, 1>& mv3){
		mv3.coeffRef(0) += -mv1.coeff(1)*mv2.coeff(5) + mv1.coeff(2)*mv2.coeff(4) - mv1.coeff(3)*mv2.coeff(2);
		mv3.coeffRef(1) +=  mv1.coeff(0)*mv2.coeff(5) - mv1.coeff(2)*mv2.coeff(3) + mv1.coeff(3)*mv2.coeff(1);
		mv3.coeffRef(2) += -mv1.coeff(0)*mv2.coeff(4) + mv1.coeff(1)*mv2.coeff(3) - mv1.coeff(3)*mv2.coeff(0);
		mv3.coeffRef(3) +=  mv1.coeff(0)*mv2.coeff(2) - mv1.coeff(1)*mv2.coeff(1) + mv1.coeff(2)*mv2.coeff(0);
	}


	/// \brief Compute the geometric product between two homogeneous multivectors mv1 (grade 3) and mv2 (grade 3). 
	/// \tparam the type of value that we manipulate, either float or double or something else.
	/// \param mv1 - the first homogeneous multivector of grade 3 represented as an Eigen::VectorXd
	/// \param mv2 - the second homogeneous multivector of grade 3 represented as a Eigen::VectorXd
	/// \param mv3 - the result of mv1 mv2 whose grade is 2
	template<typename T>
	void geometric_3_3_2(const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv1, const Eigen::Matrix<T, Eigen::Dynamic, 1>& mv2, Eigen::Matrix<T, Eigen::Dynamic, 1>& mv3){
		mv3.coeffRef(0) += -mv1.coeff(2)*mv2.coeff(3) + mv1.coeff(3)*mv2.coeff(2);
		mv3.coeffRef(1) +=  mv1.coeff(1)*mv2.coeff(3) - mv1.coeff(3)*mv2.coeff(1);
		mv3.coeffRef(2) += -mv1.coeff(0)*mv2.coeff(3) + mv1.coeff(3)*mv2.coeff(0);
		mv3.coeffRef(3) += -mv1.coeff(1)*mv2.coeff(2) + mv1.coeff(2)*mv2.coeff(1);
		mv3.coeffRef(4) +=  mv1.coeff(0)*mv2.coeff(2) - mv1.coeff(2)*mv2.coeff(0);
		mv3.coeffRef(5) += -mv1.coeff(0)*mv2.coeff(1) + mv1.coeff(1)*mv2.coeff(0);
	}


	
    template<typename T>
	std::array<std::array<std::array<std::function<void(const Eigen::Matrix<T, Eigen::Dynamic, 1> & , const Eigen::Matrix<T, Eigen::Dynamic, 1> & , Eigen::Matrix<T, Eigen::Dynamic, 1>&)>, 5>, 5>, 5> geometricFunctionsContainer =
	{{
		{{
			{{{},{},{},{},{}}},
			{{{},{},{},{},{}}},
			{{{},{},{},{},{}}},
			{{{},{},{},{},{}}},
			{{{},{},{},{},{}}}
		}},
		{{
			{{{},{},{},{},{}}},
			{{{},{},{},{},{}}},
			{{{},{},{},{},{}}},
			{{{},{},{},{},{}}},
			{{{},{},{},{},{}}}
		}},
		{{
			{{{},{},{},{},{}}},
			{{{},{},{},{},{}}},
			{{{},{},geometric_2_2_2<T>,{},{}}},
			{{{},{},{},geometric_3_3_2<T>,{}}},
			{{{},{},{},{},{}}}
		}},
		{{
			{{{},{},{},{},{}}},
			{{{},{},{},{},{}}},
			{{{},{},{},geometric_2_3_3<T>,{}}},
			{{{},{},geometric_3_2_3<T>,{},{}}},
			{{{},{},{},{},{}}}
		}},
		{{
			{{{},{},{},{},{}}},
			{{{},{},{},{},{}}},
			{{{},{},{},{},{}}},
			{{{},{},{},{},{}}},
			{{{},{},{},{},{}}}
		}}
	}};

}/// End of Namespace

#endif // E4GA_GEOMETRIC_PRODUCT_EXPLICIT_HPP__