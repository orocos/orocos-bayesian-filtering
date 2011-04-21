// $Id$
// Copyright (C) 2003 Klaas Gadeyne <first dot last at gmail dot com>
//                    Wim Meeussen  <wim dot meeussen at mech dot kuleuven dot ac dot be>
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation; either version 2.1 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
//

#include "nonlinearanalyticconditionalgaussian_ginac.h"
#include <cmath>
#include <cassert>
#include "../wrappers/rng/rng.h" // Wrapper around several rng
                                 // libraries

namespace BFL
{
  // constructor
  NonLinearAnalyticConditionalGaussian_Ginac::NonLinearAnalyticConditionalGaussian_Ginac
  (const GiNaC::matrix& func,
   const vector<GiNaC::symbol>& u,
   const vector<GiNaC::symbol>& x,
   const Gaussian& additiveNoise,
   const vector<GiNaC::symbol>& cond )
    :AnalyticConditionalGaussianAdditiveNoise(additiveNoise,3),
     func_sym        (func),
     cond_sym        (cond),
     u_sym           (u),
     x_sym           (x),
     cond_size       (cond_sym.size()),
     u_size          (u_sym.size()),
     x_size          (x_sym.size()),
     func_size       (func_sym.rows()),
     dfunc_dcond     (cond_size),
     dfunc_dx        (x_size)
  {
    // test for consistent input
    assert (func_sym.cols() == 1);
    assert (additiveNoise.DimensionGet() == cond_size);

    // derive func to cond
    for (unsigned int i=0; i < cond_size; i++)
      dfunc_dcond[i] = func_sym.diff(cond_sym[i]);

    // derive func to x
    for (unsigned int i=0; i < x_size; i++)
      dfunc_dx[i] = func_sym.diff(x_sym[i]);
  }

  // constructor
  NonLinearAnalyticConditionalGaussian_Ginac::NonLinearAnalyticConditionalGaussian_Ginac
  (const GiNaC::matrix& func,
   const vector<GiNaC::symbol>& u,
   const vector<GiNaC::symbol>& x,
   const Gaussian& additiveNoise)
    : AnalyticConditionalGaussianAdditiveNoise(additiveNoise,2),
      func_sym        (func),
      u_sym           (u),
      x_sym           (x),
      cond_size       (0),
      u_size          (u_sym.size()),
      x_size          (x_sym.size()),
      func_size       (func_sym.rows()),
      dfunc_dx        (x_size)
  {
    // test for consistent input
    assert (func_sym.cols() == 1);

    // derive func to x
    for (unsigned int i=0; i < x_size; i++)
      dfunc_dx[i] = func_sym.diff(x_sym[i]);
  }
  // copy constructor
  NonLinearAnalyticConditionalGaussian_Ginac::NonLinearAnalyticConditionalGaussian_Ginac
  (const NonLinearAnalyticConditionalGaussian_Ginac& g)
    : AnalyticConditionalGaussianAdditiveNoise( Gaussian(g.AdditiveNoiseMuGet(),g.AdditiveNoiseSigmaGet()) ,2),
      func_sym        (g.func_sym),
      cond_sym        (g.cond_sym),
      u_sym           (g.u_sym),
      x_sym           (g.x_sym),
      cond_size       (cond_sym.size()),
      u_size          (u_sym.size()),
      x_size          (x_sym.size()),
      func_size       (func_sym.rows()),
      dfunc_dcond     (cond_size),
      dfunc_dx        (x_size)
  {
    // test for consistent input
    assert (func_sym.cols() == 1);
    if (cond_size!=0) assert (g.AdditiveNoiseSigmaGet().rows() == cond_size);

    // derive func to cond
    for (unsigned int i=0; i<cond_size; i++)
      dfunc_dcond[i] = func_sym.diff(cond_sym[i]);

    // derive func to x
    for (unsigned int i=0; i < x_size; i++)
      dfunc_dx[i] = func_sym.diff(x_sym[i]);
  }

  // destructor
  NonLinearAnalyticConditionalGaussian_Ginac::~NonLinearAnalyticConditionalGaussian_Ginac()
  {}


  std::ostream& operator<< (std::ostream& os, NonLinearAnalyticConditionalGaussian_Ginac& p)
  {
    os << "function:    " << p.func_size << endl;
    os << "input:       " << p.u_size << endl;
    os << "State:       " << p.x_size << endl;
    os << "Conditional: " << ((p.cond_size) !=0) << endl;
    return os;
  }


  MatrixWrapper::ColumnVector
  NonLinearAnalyticConditionalGaussian_Ginac::ExpectedValueGet() const
  {
    MatrixWrapper::ColumnVector u_num   (u_size);
    MatrixWrapper::ColumnVector x_num   (x_size);
    MatrixWrapper::ColumnVector func_num(func_size);
    GiNaC::ex substitute (func_size);
    MatrixWrapper::ColumnVector expected(func_size);

    u_num = ConditionalArgumentGet(1);
    x_num = ConditionalArgumentGet(0);

    // use Mu of additive noise
    if (cond_size!=0)
      for (unsigned int i=0; i<u_size; i++)
	for (unsigned int j=0; j<cond_size; j++)
	  if (u_sym[i] == cond_sym[j])
	      u_num(i+1) += (this->AdditiveNoiseMuGet())(j+1);


    // evaluate func
    for (unsigned int i=0; i<func_size; i++)
      {
	// temp variable to substitute in
	substitute = func_sym(i,0);

	// substitute all u_sym with u_num
	for (unsigned int j=0; j<u_size; j++)
	  substitute = substitute.subs( u_sym[j]==u_num(j+1) );

	// substitute all x_sym with x_num
	for (unsigned int j=0; j<x_size; j++)
	  substitute = substitute.subs( x_sym[j]==x_num(j+1) );

	// build matrix func_num
	func_num(i+1) = GiNaC::ex_to<GiNaC::numeric>( substitute.evalf() ).to_double();
      }
    expected = func_num;

    if (cond_size==0)
      expected += AdditiveNoiseMuGet();

    return expected;
  }


  MatrixWrapper::SymmetricMatrix
  NonLinearAnalyticConditionalGaussian_Ginac::CovarianceGet() const
  {
    if (cond_size!=0)
      {
	MatrixWrapper::ColumnVector u_num   (u_size);
	MatrixWrapper::ColumnVector x_num   (x_size);
	GiNaC::ex substitute (func_size);
	MatrixWrapper::Matrix D             (func_size,cond_size);


	u_num = ConditionalArgumentGet(1);
	x_num = ConditionalArgumentGet(0);

	for (unsigned int i=0; i<cond_size; i++)
	  {
	    // temp variable to substitute in
	    substitute = dfunc_dcond[i];

	    // substitute all u_sym with u_num
	    for (unsigned int j=0; j<u_size; j++)
	      substitute = substitute.subs( u_sym[j]==u_num(j+1) );

	    // substitute all x_sym with x_num
	    for (unsigned int j=0; j<x_size; j++)
	      substitute = substitute.subs( x_sym[j]==x_num(j+1) );

	    // convert substitute back to matrix
	    GiNaC::matrix substitute_matrix = GiNaC::ex_to<GiNaC::matrix>(substitute);

	    // build matrix D
	    for (unsigned int j=0; j<func_size; j++)
	      D(j+1,i+1) = GiNaC::ex_to<GiNaC::numeric>( substitute_matrix(j,0).evalf() ).to_double();
	  }
	//cout << "D: " << D << endl;
	//cout << "CondCov:\n" << (Matrix)cond_covariance << endl;
	MatrixWrapper::Matrix temp = D * (MatrixWrapper::Matrix)AdditiveNoiseSigmaGet() * D.transpose();

	// convert func_covariance_matrix to symmetric matrix
	MatrixWrapper::SymmetricMatrix additiveNoise(temp.rows());
	temp.convertToSymmetricMatrix(additiveNoise);
	return additiveNoise;

      }
    else
      {
	return AdditiveNoiseSigmaGet();
      }
  }


  MatrixWrapper::Matrix
  NonLinearAnalyticConditionalGaussian_Ginac::dfGet(unsigned int i) const
  {
    // Check if i = 0, since this is the old df_dxGet method!
    assert(i == 0);

    // evaluate function
    MatrixWrapper::ColumnVector u_num   (u_size);
    MatrixWrapper::ColumnVector x_num   (x_size);
    GiNaC::ex substitute (func_size);
    MatrixWrapper::Matrix F             (func_size, x_size);

    u_num = ConditionalArgumentGet(1);
    x_num = ConditionalArgumentGet(0);

    // numeric evaluation of derivative: dfunc_dx = F
    for (unsigned int i=0; i<x_size; i++)
      {
	// temp variable to substitute in
	substitute = dfunc_dx[i];

	// substitute all u_sym with u_num
	for (unsigned int j=0; j<u_size; j++)
	  substitute = substitute.subs( u_sym[j]==u_num(j+1) );

	// substitute all x_sym with x_num
	for (unsigned int j=0; j<x_size; j++)
	  substitute = substitute.subs( x_sym[j]==x_num(j+1) );

	// convert substitute to matrix. Now all elements in matrix are accessible
	GiNaC::matrix substitute_matrix = GiNaC::ex_to<GiNaC::matrix>(substitute);

	// build matrix F
	for (unsigned int j=0; j<func_size; j++)
	  F(j+1,i+1) = GiNaC::ex_to<GiNaC::numeric>( substitute_matrix(j,0).evalf() ).to_double();
      }

    return F;
  }



  GiNaC::matrix
  NonLinearAnalyticConditionalGaussian_Ginac::FunctionGet()
  {
    return func_sym;
  }


  vector<GiNaC::symbol>
  NonLinearAnalyticConditionalGaussian_Ginac::InputGet()
  {
    return u_sym;
  }

  vector<GiNaC::symbol>
  NonLinearAnalyticConditionalGaussian_Ginac::StateGet()
  {
    return x_sym;
  }

  vector<GiNaC::symbol>
  NonLinearAnalyticConditionalGaussian_Ginac::ConditionalGet()
  {
    if ( cond_size==0 )
      return vector<GiNaC::symbol>(0);
    else
      return cond_sym;
  }

} // End namespace BFL
