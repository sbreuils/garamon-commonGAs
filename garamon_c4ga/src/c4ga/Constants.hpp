// Copyright (c) 2018 by University Paris-Est Marne-la-Vallee
// Constants.hpp
// This file is part of the Garamon for c4ga.
// Authors: Stephane Breuils and Vincent Nozick
// Contact: vincent.nozick@u-pem.fr
//
// Licence MIT
// A a copy of the MIT License is given along with this program

/// \file Constants.hpp
/// \author Stephane Breuils, Vincent Nozick
/// \brief Constant values and data related to the specified geometric algebra (c4ga)


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


#ifndef C4GA_CONSTANTS_HPP__
#define C4GA_CONSTANTS_HPP__
#pragma once


#include <array>
#include <vector>
#include <utility>
#include <string>
#include <Eigen/Sparse>

#include "c4ga/Utility.hpp"
#include "c4ga/BasisTransformations.hpp"
#include "c4ga/DualCoefficients.hpp"


/*!
 * @namespace c4ga
 */
namespace c4ga{

    constexpr unsigned int algebraDimension = 6; /*!< dimension of the algebra (number of  basis vectors of grade 1) */

    constexpr unsigned int perGradeStartingIndex[7] = {0,1,7,22,42,57,63};  /*!< array specifying the index of each first element of grade k in the full multivector */

    constexpr unsigned int binomialArray[7] = {1,6,15,20,15,6,1};  /*!< array of the (dimension + 1) first binomial coefficients */

    constexpr unsigned int xorIndexToGrade[] = {0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6}; /*!< given a Xor index in a multivector, this array indicates the corresponding grade*/ 

    constexpr unsigned int xorIndexToHomogeneousIndex[] = {0,0,1,0,2,1,5,0,3,2,6,1,9,4,10,0,4,3,7,2,10,5,11,1,12,7,13,3,16,6,10,0,5,4,8,3,11,6,12,2,13,8,14,4,17,7,11,1,14,9,15,5,18,8,12,2,19,9,13,3,14,4,5,0}; /*!< given a Xor index in a multivector, this array indicates the corresponding index in the whole homogeneous vector*/

    const std::array<std::vector<unsigned int>, 7> dualPermutations = {{ { 0}, {{ 0,4,3,2,1,5}}, {{ 6,3,1,0,10,9,8,7,14,5,4,13,2,12,11}}, {{ 7,5,4,16,2,1,13,0,11,10,9,8,19,6,18,17,3,15,14,12}}, {{ 3,2,12,1,10,9,0,7,6,5,4,14,13,11,8}}, {{ 0,4,3,2,1,5}}, {0} }}; /*!< array referring to some permutations required to compute the dual. */

    const std::array<Eigen::Matrix<double, Eigen::Dynamic,1>, 7> dualCoefficients = loadFastDualArray<double>(); /*!< array containing some basis change coefficients required to compute the dual */
    
    template<typename T>
    std::array<T, 64> recursiveDualCoefficients = {{ 1.000000,1.000000,1.000000,-1.000000,-1.000000,1.000000,-1.000000,-1.000000,1.000000,-1.000000,1.000000,1.000000,-1.000000,-1.000000,1.000000,1.000000,-1.000000,-1.000000,-1.000000,-1.000000,1.000000,1.000000,-1.000000,-1.000000,-1.000000,-1.000000,1.000000,-1.000000,-1.000000,1.000000,1.000000,1.000000,-1.000000,1.000000,1.000000,-1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,-1.000000,-1.000000,1.000000,1.000000,-1.000000,1.000000,-1.000000,-1.000000,1.000000,1.000000,-1.000000,-1.000000,-1.000000,-1.000000,1.000000,1.000000,-1.000000,-1.000000,-1.000000,1.000000,1.000000,-1.000000,1.000000}}; /*!< array containing the coefficients needed to compute the recursive product like (primal^dual) */

    constexpr double pseudoScalarInverse = 1.000000; /*!< compute the inverse of the pseudo scalar */

    const int signReversePerGrade[7] = {1,1,-1,-1,1,1,-1}; /*!< array of signs to avoid the computation of (-1)^k*(k-1)/2 during the reverse operation */

    const std::vector<std::string> basisVectors = {"0", "1", "2", "3", "4", "i"}; /*!< name of the basis vectors (of grade 1) */

    const std::string metric =
"\
	e0	e1	e2	e3	e4	ei	\n\
e0	0	0	0	0	0	-1	\n\
e1	0	1	0	0	0	0	\n\
e2	0	0	1	0	0	0	\n\
e3	0	0	0	1	0	0	\n\
e4	0	0	0	0	1	0	\n\
ei	-1	0	0	0	0	0	\n\
"; /*!< metric / quadratic form of the algebra (inner product between basis vectors) */


    template<class T>
    const T zero = 0;

    const unsigned int scalar = 0;
    const unsigned int E0 = 1;
    const unsigned int E1 = 2;
    const unsigned int E2 = 4;
    const unsigned int E3 = 8;
    const unsigned int E4 = 16;
    const unsigned int Ei = 32;
    const unsigned int E01 = 3;
    const unsigned int E02 = 5;
    const unsigned int E03 = 9;
    const unsigned int E04 = 17;
    const unsigned int E0i = 33;
    const unsigned int E12 = 6;
    const unsigned int E13 = 10;
    const unsigned int E14 = 18;
    const unsigned int E1i = 34;
    const unsigned int E23 = 12;
    const unsigned int E24 = 20;
    const unsigned int E2i = 36;
    const unsigned int E34 = 24;
    const unsigned int E3i = 40;
    const unsigned int E4i = 48;
    const unsigned int E012 = 7;
    const unsigned int E013 = 11;
    const unsigned int E014 = 19;
    const unsigned int E01i = 35;
    const unsigned int E023 = 13;
    const unsigned int E024 = 21;
    const unsigned int E02i = 37;
    const unsigned int E034 = 25;
    const unsigned int E03i = 41;
    const unsigned int E04i = 49;
    const unsigned int E123 = 14;
    const unsigned int E124 = 22;
    const unsigned int E12i = 38;
    const unsigned int E134 = 26;
    const unsigned int E13i = 42;
    const unsigned int E14i = 50;
    const unsigned int E234 = 28;
    const unsigned int E23i = 44;
    const unsigned int E24i = 52;
    const unsigned int E34i = 56;
    const unsigned int E0123 = 15;
    const unsigned int E0124 = 23;
    const unsigned int E012i = 39;
    const unsigned int E0134 = 27;
    const unsigned int E013i = 43;
    const unsigned int E014i = 51;
    const unsigned int E0234 = 29;
    const unsigned int E023i = 45;
    const unsigned int E024i = 53;
    const unsigned int E034i = 57;
    const unsigned int E1234 = 30;
    const unsigned int E123i = 46;
    const unsigned int E124i = 54;
    const unsigned int E134i = 58;
    const unsigned int E234i = 60;
    const unsigned int E01234 = 31;
    const unsigned int E0123i = 47;
    const unsigned int E0124i = 55;
    const unsigned int E0134i = 59;
    const unsigned int E0234i = 61;
    const unsigned int E1234i = 62;
    const unsigned int E01234i = 63;
    /*!< defines the constants for the cga */

    template<typename T>
    const std::array<Eigen::SparseMatrix<T, Eigen::ColMajor>,7> transformationMatrices = loadMatrices<T>(); /*!< set transformation matrices to transform a k-vector from the orhogonal basis to the original basis */
    template<typename T>
    const std::array<Eigen::SparseMatrix<T, Eigen::ColMajor>,7> transformationMatricesInverse = loadMatricesInverse<T>(); /*!< set transformation matrices to transform a k-vector from the original basis to the orhogonal basis */


    template<typename T>
    Eigen::Matrix<T, Eigen::Dynamic, 1> diagonalMetric = []()->Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>{Eigen::Matrix<T, Eigen::Dynamic, 1> tmp(6);tmp<<2.000000,-2.000000,1.000000,1.000000,1.000000,1.000000; return tmp;}();   /*!< defines the diagonal metric (stored as a vector) */


}  // namespace


#endif // C4GA_CONSTANTS_HPP__
