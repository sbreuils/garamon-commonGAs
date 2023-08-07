// Copyright (c) 2018 by University Paris-Est Marne-la-Vallee
// Constants.hpp
// This file is part of the Garamon for e2ga.
// Authors: Stephane Breuils and Vincent Nozick
// Contact: vincent.nozick@u-pem.fr
//
// Licence MIT
// A a copy of the MIT License is given along with this program

/// \file Constants.hpp
/// \author Stephane Breuils, Vincent Nozick
/// \brief Constant values and data related to the specified geometric algebra (e2ga)


// Doxygen
/// \version 0.1
/// \mainpage
/// \tableofcontents
/// \section instroduction_sec What for?
/// Garamon is a C++ library to represent an manipulate the geometric algebra objects.
/// \section install_bigsec How to install
/// \subsection dependencies_sec Dependecies
/// \li cmake (at least version 3.10)
/// \li Eigen (at least version 3.3.4)
/// \li Doxygen (if you want the documentation)
/// \subsection install_sec Install with cmake (Linux / MacOs)
/// \li go to garamon dir
/// \li mkdir build
/// \li cd build
/// \li cmake ..
/// \li make
/// \li (optional) sudo make install
/// \li (optional for documentation) make html
/// \li The documentation is located in [path to build]/doc/doc-doxygen/html/index.html


#ifndef E2GA_CONSTANTS_HPP__
#define E2GA_CONSTANTS_HPP__
#pragma once


#include <array>
#include <vector>
#include <utility>
#include <string>
#include <Eigen/Sparse>

#include "e2ga/Utility.hpp"
#include "e2ga/BasisTransformations.hpp"
#include "e2ga/DualCoefficients.hpp"


/*!
 * @namespace e2ga
 */
namespace e2ga{

    constexpr unsigned int algebraDimension = 2; /*!< dimension of the algebra (number of  basis vectors of grade 1) */

    constexpr unsigned int perGradeStartingIndex[3] = {0,1,3};  /*!< array specifying the index of each first element of grade k in the full multivector */

    constexpr unsigned int binomialArray[3] = {1,2,1};  /*!< array of the (dimension + 1) first binomial coefficients */

    constexpr unsigned int xorIndexToGrade[] = {0,1,1,2}; /*!< given a Xor index in a multivector, this array indicates the corresponding grade*/ 

    constexpr unsigned int xorIndexToHomogeneousIndex[] = {0,0,1,0}; /*!< given a Xor index in a multivector, this array indicates the corresponding index in the whole homogeneous vector*/

    std::array<std::vector<unsigned int>, 3> dualPermutations = {{ { 0}, {{ 1, 0}}, {0} }};

    const std::array<Eigen::Matrix<double, Eigen::Dynamic,1>, 3> dualCoefficients = loadFastDualArray<double>(); /*!< array containing the coefficients needed to compute the dual */ 
    
    template<typename T>
    std::array<T, 4> recursiveDualCoefficients = {{ 1.000000, -1.000000, 1.000000, 1.000000}}; /*!< array containing the coefficients needed to compute the recursive product like (primal^dual) */
    

    constexpr double pseudoScalarInverse = -1.000000; /*!< compute the inverse of the pseudo scalar */

    const int signReversePerGrade[3] = {1,1,-1}; /*!< array of signs to avoid the computation of (-1)^k*(k-1)/2 during the reverse operation */

    const std::vector<std::string> basisVectors = {"1", "2"}; /*!< name of the basis vectors (of grade 1) */

    const std::string metric =
"\
	e1	e2	\n\
e1	1	0	\n\
e2	0	1	\n\
"; /*!< metric / quadratic form of the algebra (inner product between basis vectors) */


    template<class T>
    const T zero = 0;

    const unsigned int scalar = 0;
    const unsigned int E1 = 1;
    const unsigned int E2 = 2;
    const unsigned int E12 = 3;
    /*!< defines the constants for the cga */



    template<typename T>
    Eigen::Matrix<T, Eigen::Dynamic, 1> diagonalMetric = []()->Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>{Eigen::Matrix<T, Eigen::Dynamic, 1> tmp(2);tmp<<1.000000,1.000000; return tmp;}();   /*!< defines the diagonal metric (stored as a vector) */


}  // namespace


#endif // E2GA_CONSTANTS_HPP__
