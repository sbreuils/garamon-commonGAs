// Copyright (c) 2018 by University Paris-Est Marne-la-Vallee
// Inner.hpp
// This file is part of the Garamon for c4ga.
// Authors: Stephane Breuils and Vincent Nozick
// Contact: vincent.nozick@u-pem.fr
//
// Licence MIT
// A a copy of the MIT License is given along with this program

// \file BasisTransformations.hpp
// \author V. Nozick, S. Breuils
// \brief this files generates and load the elements of the transformation matrices into array of Eigen matrices 


#ifndef C4GA_BASISTRANSFORMATIONS_HPP__
#define C4GA_BASISTRANSFORMATIONS_HPP__
#pragma once


#include <Eigen/Sparse>
#include <vector>
#include <array>
#include <algorithm>
#include <iterator>
#include <sstream>


/*!
 * @namespace c4ga
 */
namespace c4ga{

	/// Decode the string encodedPerGradeMatrixComponents to a vector of T. This vector will contain the components of the transformation matrices
	template<typename T>
	const std::vector<T> decodeStringToVecOfT(std::string encodedPerGradeMatrixComponents) {
	    std::vector<T> resultDecodedGrade;
	    std::istringstream streamStringOverFloat(encodedPerGradeMatrixComponents);
	    std::copy(std::istream_iterator<float>(streamStringOverFloat),
		std::istream_iterator<float>(),
		std::back_inserter(resultDecodedGrade));
	    return resultDecodedGrade;
	}

	/// contains and load the components of the grade 0 direct transformation matrix 
 	template<typename T>
	Eigen::SparseMatrix<T, Eigen::ColMajor>  loadgrade0Matrix() {
		const std::string grade0MatrixComponents = " 0 0 1.000000 ";
		std::vector<T> gradeVectorComponents = decodeStringToVecOfT<T>(grade0MatrixComponents);
		Eigen::SparseMatrix<T, Eigen::ColMajor> perGradeBasisTransformMatrix		 = Eigen::SparseMatrix<T, Eigen::ColMajor>(1,1);
		perGradeBasisTransformMatrix.reserve(1);
		perGradeBasisTransformMatrix.insert((int)gradeVectorComponents[0],(int)gradeVectorComponents[1]) = gradeVectorComponents[2];
		return perGradeBasisTransformMatrix;
	}

	/// contains and load the components of the grade 1 direct transformation matrix 
 	template<typename T>
	Eigen::SparseMatrix<T, Eigen::ColMajor>  loadgrade1Matrix() {
		const std::string grade1MatrixComponents = " 0 0 1.000000 5 0 -1.000000 0 1 1.000000 5 1 1.000000 2 2 1.000000 3 3 1.000000 4 4 1.000000 1 5 1.000000 ";
		std::vector<T> gradeVectorComponents = decodeStringToVecOfT<T>(grade1MatrixComponents);
		Eigen::SparseMatrix<T, Eigen::ColMajor> perGradeBasisTransformMatrix		 = Eigen::SparseMatrix<T, Eigen::ColMajor>(6,6);
		perGradeBasisTransformMatrix.reserve(8);
		for(unsigned int i=0;i<8;++i)
			perGradeBasisTransformMatrix.insert((int)gradeVectorComponents[3*i],(int)gradeVectorComponents[(3*i)+1]) = gradeVectorComponents[(3*i)+2];
		return perGradeBasisTransformMatrix;
	}

	/// contains and load the components of the grade 2 direct transformation matrix 
 	template<typename T>
	Eigen::SparseMatrix<T, Eigen::ColMajor>  loadgrade2Matrix() {
		const std::string grade2MatrixComponents = " 4 0 2.000000 1 1 1.000000 11 1 1.000000 2 2 1.000000 13 2 1.000000 3 3 1.000000 14 3 1.000000 0 4 1.000000 8 4 1.000000 1 5 1.000000 11 5 -1.000000 2 6 1.000000 13 6 -1.000000 3 7 1.000000 14 7 -1.000000 0 8 1.000000 8 8 -1.000000 9 9 1.000000 10 10 1.000000 5 11 -1.000000 12 12 1.000000 6 13 -1.000000 7 14 -1.000000 ";
		std::vector<T> gradeVectorComponents = decodeStringToVecOfT<T>(grade2MatrixComponents);
		Eigen::SparseMatrix<T, Eigen::ColMajor> perGradeBasisTransformMatrix		 = Eigen::SparseMatrix<T, Eigen::ColMajor>(15,15);
		perGradeBasisTransformMatrix.reserve(23);
		for(unsigned int i=0;i<23;++i)
			perGradeBasisTransformMatrix.insert((int)gradeVectorComponents[3*i],(int)gradeVectorComponents[(3*i)+1]) = gradeVectorComponents[(3*i)+2];
		return perGradeBasisTransformMatrix;
	}

	/// contains and load the components of the grade 3 direct transformation matrix 
 	template<typename T>
	Eigen::SparseMatrix<T, Eigen::ColMajor>  loadgrade3Matrix() {
		const std::string grade3MatrixComponents = " 6 0 -2.000000 8 1 -2.000000 9 2 -2.000000 3 3 -2.000000 4 4 1.000000 17 4 -1.000000 5 5 1.000000 18 5 -1.000000 0 6 -1.000000 12 6 1.000000 7 7 1.000000 19 7 -1.000000 1 8 -1.000000 14 8 1.000000 2 9 -1.000000 15 9 1.000000 4 10 1.000000 17 10 1.000000 5 11 1.000000 18 11 1.000000 0 12 -1.000000 12 12 -1.000000 7 13 1.000000 19 13 1.000000 1 14 -1.000000 14 14 -1.000000 2 15 -1.000000 15 15 -1.000000 16 16 1.000000 10 17 1.000000 11 18 1.000000 13 19 1.000000 ";
		std::vector<T> gradeVectorComponents = decodeStringToVecOfT<T>(grade3MatrixComponents);
		Eigen::SparseMatrix<T, Eigen::ColMajor> perGradeBasisTransformMatrix		 = Eigen::SparseMatrix<T, Eigen::ColMajor>(20,20);
		perGradeBasisTransformMatrix.reserve(32);
		for(unsigned int i=0;i<32;++i)
			perGradeBasisTransformMatrix.insert((int)gradeVectorComponents[3*i],(int)gradeVectorComponents[(3*i)+1]) = gradeVectorComponents[(3*i)+2];
		return perGradeBasisTransformMatrix;
	}

	/// contains and load the components of the grade 4 direct transformation matrix 
 	template<typename T>
	Eigen::SparseMatrix<T, Eigen::ColMajor>  loadgrade4Matrix() {
		const std::string grade4MatrixComponents = " 7 0 2.000000 8 1 2.000000 2 2 -2.000000 9 3 2.000000 4 4 -2.000000 5 5 -2.000000 6 6 1.000000 14 6 1.000000 0 7 1.000000 11 7 1.000000 1 8 1.000000 12 8 1.000000 3 9 1.000000 13 9 1.000000 6 10 1.000000 14 10 -1.000000 0 11 1.000000 11 11 -1.000000 1 12 1.000000 12 12 -1.000000 3 13 1.000000 13 13 -1.000000 10 14 -1.000000 ";
		std::vector<T> gradeVectorComponents = decodeStringToVecOfT<T>(grade4MatrixComponents);
		Eigen::SparseMatrix<T, Eigen::ColMajor> perGradeBasisTransformMatrix		 = Eigen::SparseMatrix<T, Eigen::ColMajor>(15,15);
		perGradeBasisTransformMatrix.reserve(23);
		for(unsigned int i=0;i<23;++i)
			perGradeBasisTransformMatrix.insert((int)gradeVectorComponents[3*i],(int)gradeVectorComponents[(3*i)+1]) = gradeVectorComponents[(3*i)+2];
		return perGradeBasisTransformMatrix;
	}

	/// contains and load the components of the grade 5 direct transformation matrix 
 	template<typename T>
	Eigen::SparseMatrix<T, Eigen::ColMajor>  loadgrade5Matrix() {
		const std::string grade5MatrixComponents = " 4 0 -2.000000 1 1 -2.000000 2 2 -2.000000 3 3 -2.000000 0 4 -1.000000 5 4 1.000000 0 5 -1.000000 5 5 -1.000000 ";
		std::vector<T> gradeVectorComponents = decodeStringToVecOfT<T>(grade5MatrixComponents);
		Eigen::SparseMatrix<T, Eigen::ColMajor> perGradeBasisTransformMatrix		 = Eigen::SparseMatrix<T, Eigen::ColMajor>(6,6);
		perGradeBasisTransformMatrix.reserve(8);
		for(unsigned int i=0;i<8;++i)
			perGradeBasisTransformMatrix.insert((int)gradeVectorComponents[3*i],(int)gradeVectorComponents[(3*i)+1]) = gradeVectorComponents[(3*i)+2];
		return perGradeBasisTransformMatrix;
	}

	/// contains and load the components of the grade 6 direct transformation matrix 
 	template<typename T>
	Eigen::SparseMatrix<T, Eigen::ColMajor>  loadgrade6Matrix() {
		const std::string grade6MatrixComponents = " 0 0 -2.000000 ";
		std::vector<T> gradeVectorComponents = decodeStringToVecOfT<T>(grade6MatrixComponents);
		Eigen::SparseMatrix<T, Eigen::ColMajor> perGradeBasisTransformMatrix		 = Eigen::SparseMatrix<T, Eigen::ColMajor>(1,1);
		perGradeBasisTransformMatrix.reserve(1);
		perGradeBasisTransformMatrix.insert((int)gradeVectorComponents[0],(int)gradeVectorComponents[1]) = gradeVectorComponents[2];
		return perGradeBasisTransformMatrix;
	}




	/// contains and load the components of the grade 0 inverse transformation matrix 
 	template<typename T>
	Eigen::SparseMatrix<T, Eigen::ColMajor>  loadgrade0InverseMatrix() {
		const std::string grade0MatrixComponents = " 0 0 1.000000 ";
		std::vector<T> gradeVectorComponents = decodeStringToVecOfT<T>(grade0MatrixComponents);
		Eigen::SparseMatrix<T, Eigen::ColMajor> perGradeBasisTransformMatrix		 = Eigen::SparseMatrix<T, Eigen::ColMajor>(1,1);
		perGradeBasisTransformMatrix.reserve(1);
		perGradeBasisTransformMatrix.insert((int)gradeVectorComponents[0],(int)gradeVectorComponents[1]) = gradeVectorComponents[2];
		return perGradeBasisTransformMatrix;
	}

	/// contains and load the components of the grade 1 inverse transformation matrix 
 	template<typename T>
	Eigen::SparseMatrix<T, Eigen::ColMajor>  loadgrade1InverseMatrix() {
		const std::string grade1MatrixComponents = " 0 0 0.500000 1 0 0.500000 5 1 1.000000 2 2 1.000000 3 3 1.000000 4 4 1.000000 0 5 -0.500000 1 5 0.500000 ";
		std::vector<T> gradeVectorComponents = decodeStringToVecOfT<T>(grade1MatrixComponents);
		Eigen::SparseMatrix<T, Eigen::ColMajor> perGradeBasisTransformMatrix		 = Eigen::SparseMatrix<T, Eigen::ColMajor>(6,6);
		perGradeBasisTransformMatrix.reserve(8);
		for(unsigned int i=0;i<8;++i)
			perGradeBasisTransformMatrix.insert((int)gradeVectorComponents[3*i],(int)gradeVectorComponents[(3*i)+1]) = gradeVectorComponents[(3*i)+2];
		return perGradeBasisTransformMatrix;
	}

	/// contains and load the components of the grade 2 inverse transformation matrix 
 	template<typename T>
	Eigen::SparseMatrix<T, Eigen::ColMajor>  loadgrade2InverseMatrix() {
		const std::string grade2MatrixComponents = " 4 0 0.500000 8 0 0.500000 1 1 0.500000 5 1 0.500000 2 2 0.500000 6 2 0.500000 3 3 0.500000 7 3 0.500000 0 4 0.500000 11 5 -1.000000 13 6 -1.000000 14 7 -1.000000 4 8 0.500000 8 8 -0.500000 9 9 1.000000 10 10 1.000000 1 11 0.500000 5 11 -0.500000 12 12 1.000000 2 13 0.500000 6 13 -0.500000 3 14 0.500000 7 14 -0.500000 ";
		std::vector<T> gradeVectorComponents = decodeStringToVecOfT<T>(grade2MatrixComponents);
		Eigen::SparseMatrix<T, Eigen::ColMajor> perGradeBasisTransformMatrix		 = Eigen::SparseMatrix<T, Eigen::ColMajor>(15,15);
		perGradeBasisTransformMatrix.reserve(23);
		for(unsigned int i=0;i<23;++i)
			perGradeBasisTransformMatrix.insert((int)gradeVectorComponents[3*i],(int)gradeVectorComponents[(3*i)+1]) = gradeVectorComponents[(3*i)+2];
		return perGradeBasisTransformMatrix;
	}

	/// contains and load the components of the grade 3 inverse transformation matrix 
 	template<typename T>
	Eigen::SparseMatrix<T, Eigen::ColMajor>  loadgrade3InverseMatrix() {
		const std::string grade3MatrixComponents = " 6 0 -0.500000 12 0 -0.500000 8 1 -0.500000 14 1 -0.500000 9 2 -0.500000 15 2 -0.500000 3 3 -0.500000 4 4 0.500000 10 4 0.500000 5 5 0.500000 11 5 0.500000 0 6 -0.500000 7 7 0.500000 13 7 0.500000 1 8 -0.500000 2 9 -0.500000 17 10 1.000000 18 11 1.000000 6 12 0.500000 12 12 -0.500000 19 13 1.000000 8 14 0.500000 14 14 -0.500000 9 15 0.500000 15 15 -0.500000 16 16 1.000000 4 17 -0.500000 10 17 0.500000 5 18 -0.500000 11 18 0.500000 7 19 -0.500000 13 19 0.500000 ";
		std::vector<T> gradeVectorComponents = decodeStringToVecOfT<T>(grade3MatrixComponents);
		Eigen::SparseMatrix<T, Eigen::ColMajor> perGradeBasisTransformMatrix		 = Eigen::SparseMatrix<T, Eigen::ColMajor>(20,20);
		perGradeBasisTransformMatrix.reserve(32);
		for(unsigned int i=0;i<32;++i)
			perGradeBasisTransformMatrix.insert((int)gradeVectorComponents[3*i],(int)gradeVectorComponents[(3*i)+1]) = gradeVectorComponents[(3*i)+2];
		return perGradeBasisTransformMatrix;
	}

	/// contains and load the components of the grade 4 inverse transformation matrix 
 	template<typename T>
	Eigen::SparseMatrix<T, Eigen::ColMajor>  loadgrade4InverseMatrix() {
		const std::string grade4MatrixComponents = " 7 0 0.500000 11 0 0.500000 8 1 0.500000 12 1 0.500000 2 2 -0.500000 9 3 0.500000 13 3 0.500000 4 4 -0.500000 5 5 -0.500000 6 6 0.500000 10 6 0.500000 0 7 0.500000 1 8 0.500000 3 9 0.500000 14 10 -1.000000 7 11 0.500000 11 11 -0.500000 8 12 0.500000 12 12 -0.500000 9 13 0.500000 13 13 -0.500000 6 14 0.500000 10 14 -0.500000 ";
		std::vector<T> gradeVectorComponents = decodeStringToVecOfT<T>(grade4MatrixComponents);
		Eigen::SparseMatrix<T, Eigen::ColMajor> perGradeBasisTransformMatrix		 = Eigen::SparseMatrix<T, Eigen::ColMajor>(15,15);
		perGradeBasisTransformMatrix.reserve(23);
		for(unsigned int i=0;i<23;++i)
			perGradeBasisTransformMatrix.insert((int)gradeVectorComponents[3*i],(int)gradeVectorComponents[(3*i)+1]) = gradeVectorComponents[(3*i)+2];
		return perGradeBasisTransformMatrix;
	}

	/// contains and load the components of the grade 5 inverse transformation matrix 
 	template<typename T>
	Eigen::SparseMatrix<T, Eigen::ColMajor>  loadgrade5InverseMatrix() {
		const std::string grade5MatrixComponents = " 4 0 -0.500000 5 0 -0.500000 1 1 -0.500000 2 2 -0.500000 3 3 -0.500000 0 4 -0.500000 4 5 0.500000 5 5 -0.500000 ";
		std::vector<T> gradeVectorComponents = decodeStringToVecOfT<T>(grade5MatrixComponents);
		Eigen::SparseMatrix<T, Eigen::ColMajor> perGradeBasisTransformMatrix		 = Eigen::SparseMatrix<T, Eigen::ColMajor>(6,6);
		perGradeBasisTransformMatrix.reserve(8);
		for(unsigned int i=0;i<8;++i)
			perGradeBasisTransformMatrix.insert((int)gradeVectorComponents[3*i],(int)gradeVectorComponents[(3*i)+1]) = gradeVectorComponents[(3*i)+2];
		return perGradeBasisTransformMatrix;
	}

	/// contains and load the components of the grade 6 inverse transformation matrix 
 	template<typename T>
	Eigen::SparseMatrix<T, Eigen::ColMajor>  loadgrade6InverseMatrix() {
		const std::string grade6MatrixComponents = " 0 0 -0.500000 ";
		std::vector<T> gradeVectorComponents = decodeStringToVecOfT<T>(grade6MatrixComponents);
		Eigen::SparseMatrix<T, Eigen::ColMajor> perGradeBasisTransformMatrix		 = Eigen::SparseMatrix<T, Eigen::ColMajor>(1,1);
		perGradeBasisTransformMatrix.reserve(1);
		perGradeBasisTransformMatrix.insert((int)gradeVectorComponents[0],(int)gradeVectorComponents[1]) = gradeVectorComponents[2];
		return perGradeBasisTransformMatrix;
	}




	/// initialize all the direct transformation matrices using array of eigen sparse matrices
	template<typename T>
	const std::array<Eigen::SparseMatrix<T, Eigen::ColMajor>,7> loadMatrices() {
		std::array<Eigen::SparseMatrix<T, Eigen::ColMajor>,7>  transformationMatrices;
		transformationMatrices[0] = loadgrade0Matrix<T>();
		transformationMatrices[1] = loadgrade1Matrix<T>();
		transformationMatrices[2] = loadgrade2Matrix<T>();
		transformationMatrices[3] = loadgrade3Matrix<T>();
		transformationMatrices[4] = loadgrade4Matrix<T>();
		transformationMatrices[5] = loadgrade5Matrix<T>();
		transformationMatrices[6] = loadgrade6Matrix<T>();
		return transformationMatrices;
	}


	/// initialize all the inverse transformation matrices using array of eigen sparse matrices
	template<typename T>
	const std::array<Eigen::SparseMatrix<T, Eigen::ColMajor>,7> loadMatricesInverse() {
		std::array<Eigen::SparseMatrix<T, Eigen::ColMajor>,7>  transformationMatrices;
		transformationMatrices[0] = loadgrade0InverseMatrix<T>();
		transformationMatrices[1] = loadgrade1InverseMatrix<T>();
		transformationMatrices[2] = loadgrade2InverseMatrix<T>();
		transformationMatrices[3] = loadgrade3InverseMatrix<T>();
		transformationMatrices[4] = loadgrade4InverseMatrix<T>();
		transformationMatrices[5] = loadgrade5InverseMatrix<T>();
		transformationMatrices[6] = loadgrade6InverseMatrix<T>();
		return transformationMatrices;
	}


}/// End of Namespace

#endif // C4GA_BASISTRANSFORMATIONS_HPP__
