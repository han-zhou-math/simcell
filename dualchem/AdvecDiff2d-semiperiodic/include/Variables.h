/*=============================================================================
*   
*   Filename : Variables.h
*   Creator : Han Zhou
*   Date : 07/22/21
*   Description : 
*
=============================================================================*/
 
#ifndef _VARIABLES_H
#define _VARIABLES_H

#include "Array.H"
#include "MyVector.H"
#include "MyMatrix.H"
#include "MyTensor.H"
//#include "AbstractMatrix.H"
//#include "SparseTensor.H"

//=============================================================================

typedef Array< double, 2 > Array2d;
typedef Array< int, 2 > Array2i;
typedef Array< bool, 2 > Array2b;

typedef Array< double, 3 > Array3d;
typedef Array< int, 3 > Array3i;
typedef Array< bool, 3 > Array3b;

typedef Array< double, 4 > Array4d;
typedef Array< int, 4 > Array4i;
typedef Array< bool, 4 > Array4b;

typedef Array< double, 6 > Array6d;
typedef Array< int, 6 > Array6i;
typedef Array< bool, 6 > Array6b;

typedef Array< double, 7 > Array7d;
typedef Array< int, 7 > Array7i;
typedef Array< bool, 7 > Array7b;

typedef Array< double, 9 > Array9d;
typedef Array< int, 9 > Array9i;
typedef Array< bool, 9 > Array9b;

typedef Array< double, 10 > Array10d;
typedef Array< int, 10 > Array10i;
typedef Array< bool, 10 > Array10b;

typedef Array< double, 12 > Array12d;
typedef Array< int, 12 > Array12i;
typedef Array< bool, 12 > Array12b;

//=============================================================================

typedef Vector< double > VectorXd;
typedef Vector< float > VectorXf;
typedef Vector< int > VectorXi;
typedef Vector< bool > VectorXb;

typedef Vector< Array2d > VectorX2d;
typedef Vector< Array2i > VectorX2i;
typedef Vector< Array2b > VectorX2b;

typedef Vector< Array3d > VectorX3d;
typedef Vector< Array3i > VectorX3i;
typedef Vector< Array3b > VectorX3b;

typedef Vector< Array4d > VectorX4d;
typedef Vector< Array4i > VectorX4i;
typedef Vector< Array4b > VectorX4b;

typedef Vector< Array6d > VectorX6d;
typedef Vector< Array6i > VectorX6i;
typedef Vector< Array6b > VectorX6b;

typedef Vector< Array7d > VectorX7d;
typedef Vector< Array7i > VectorX7i;
typedef Vector< Array7b > VectorX7b;

typedef Vector< Array9d > VectorX9d;
typedef Vector< Array9i > VectorX9i;
typedef Vector< Array9b > VectorX9b;

typedef Vector< Array10d > VectorX10d;
typedef Vector< Array10i > VectorX10i;
typedef Vector< Array10b > VectorX10b;

typedef Vector< Array12d > VectorX12d;
typedef Vector< Array12i > VectorX12i;
typedef Vector< Array12b > VectorX12b;

//=============================================================================

typedef Matrix< double > MatrixXd;
typedef Matrix< float > MatrixXf;
typedef Matrix< int > MatrixXi;
typedef Matrix< bool > MatrixXb;

typedef Matrix< Array2d > MatrixX2d;
typedef Matrix< Array2i > MatrixX2i;
typedef Matrix< Array2b > MatrixX2b;

typedef Matrix< Array3d > MatrixX3d;
typedef Matrix< Array3i > MatrixX3i;
typedef Matrix< Array3b > MatrixX3b;

typedef Matrix< Array6d > MatrixX6d;
typedef Matrix< Array6i > MatrixX6i;
typedef Matrix< Array6b > MatrixX6b;

//=============================================================================

typedef Tensor< double > TensorXd;
typedef Tensor< int > TensorXi;
typedef Tensor< bool > TensorXb;

typedef Tensor< Array3d > TensorX3d;
typedef Tensor< Array3i > TensorX3i;
typedef Tensor< Array3b > TensorX3b;

//=============================================================================
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#endif
