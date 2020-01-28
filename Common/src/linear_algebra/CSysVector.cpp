/*!
 * \file vector_structure.cpp
 * \brief Main classes required for solving linear systems of equations
 * \author F. Palacios, J. Hicken
 * \version 7.0.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation 
 * (http://su2foundation.org)
 *
 * Copyright 2012-2019, SU2 Contributors (cf. AUTHORS.md)
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

#include "../../include/linear_algebra/CSysVector.hpp"

template<class ScalarType>
CSysVector<ScalarType>::CSysVector(void) {

  nElm = 0; nElmDomain = 0;
  nBlk = 0; nBlkDomain = 0;
  nVar = 0;
  
  vec_val = NULL;
  nElm = 0;
  nElmDomain = 0;
  nVar = 0;
  nBlk = 0;
  nBlkDomain = 0;

}

template<class ScalarType>
CSysVector<ScalarType>::CSysVector(const unsigned long & size, const ScalarType & val) {

  nElm = size; nElmDomain = size;
  nBlk = nElm; nBlkDomain = nElmDomain;
  nVar = 1;

  /*--- Check for invalid size, then allocate memory and initialize values ---*/
  if ( (nElm >= ULONG_MAX) ) {
    char buf[100];
    SPRINTF(buf, "Invalid input: size = %lu", size );
    SU2_MPI::Error(string(buf), CURRENT_FUNCTION);
  }

  vec_val = new ScalarType[nElm];
  for (unsigned int i = 0; i < nElm; i++)
    vec_val[i] = val;

#ifdef HAVE_MPI
  unsigned long nElmLocal = (unsigned long)nElm;
  SU2_MPI::Allreduce(&nElmLocal, &nElmGlobal, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif

}

template<class ScalarType>
CSysVector<ScalarType>::CSysVector(const unsigned long & numBlk, const unsigned long & numBlkDomain, const unsigned short & numVar,
                       const ScalarType & val) {

  nElm = numBlk*numVar; nElmDomain = numBlkDomain*numVar;
  nBlk = numBlk; nBlkDomain = numBlkDomain;
  nVar = numVar;

  /*--- Check for invalid input, then allocate memory and initialize values ---*/
  if ( nElm >= ULONG_MAX ) {
    char buf[100];
    SPRINTF(buf, "invalid input: numBlk, numVar = %lu, %u", numBlk, numVar );
    SU2_MPI::Error(string(buf), CURRENT_FUNCTION);
  }

  vec_val = new ScalarType[nElm];
  for (unsigned int i = 0; i < nElm; i++)
    vec_val[i] = val;

#ifdef HAVE_MPI
  unsigned long nElmLocal = (unsigned long)nElm;
  SU2_MPI::Allreduce(&nElmLocal, &nElmGlobal, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif

}

template<class ScalarType>
CSysVector<ScalarType>::CSysVector(const CSysVector<ScalarType> & u) {

  /*--- Copy size information, allocate memory, and initialize values ---*/
  nElm = u.nElm; nElmDomain = u.nElmDomain;
  nBlk = u.nBlk; nBlkDomain = u.nBlkDomain;
  nVar = u.nVar;

  vec_val = new ScalarType[nElm];
  for (unsigned long i = 0; i < nElm; i++)
    vec_val[i] = u.vec_val[i];

#ifdef HAVE_MPI
  nElmGlobal = u.nElmGlobal;
#endif

}

template<class ScalarType>
CSysVector<ScalarType>::CSysVector(const unsigned long & size, const ScalarType* u_array) {

  nElm = size; nElmDomain = size;
  nBlk = nElm; nBlkDomain = nElmDomain;
  nVar = 1;

  /*--- Check for invalid size, then allocate memory and initialize values ---*/
  if ( nElm >= ULONG_MAX ) {
    char buf[100];
    SPRINTF(buf, "Invalid input: size = %lu", size );
    SU2_MPI::Error(string(buf), CURRENT_FUNCTION);
  }

  vec_val = new ScalarType[nElm];
  for (unsigned long i = 0; i < nElm; i++)
    vec_val[i] = u_array[i];

#ifdef HAVE_MPI
  unsigned long nElmLocal = (unsigned long)nElm;
  SU2_MPI::Allreduce(&nElmLocal, &nElmGlobal, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif

}

template<class ScalarType>
CSysVector<ScalarType>::CSysVector(const unsigned long & numBlk, const unsigned long & numBlkDomain, const unsigned short & numVar,
                       const ScalarType* u_array) {

  nElm = numBlk*numVar; nElmDomain = numBlkDomain*numVar;
  nBlk = numBlk; nBlkDomain = numBlkDomain;
  nVar = numVar;

  /*--- check for invalid input, then allocate memory and initialize values ---*/
  if ( nElm >= ULONG_MAX ) {
    char buf[100];
    SPRINTF(buf, "invalid input: numBlk, numVar = %lu, %u", numBlk, numVar );
    SU2_MPI::Error(string(buf), CURRENT_FUNCTION);
  }

  vec_val = new ScalarType[nElm];
  for (unsigned long i = 0; i < nElm; i++)
    vec_val[i] = u_array[i];

#ifdef HAVE_MPI
  unsigned long nElmLocal = (unsigned long)nElm;
  SU2_MPI::Allreduce(&nElmLocal, &nElmGlobal, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif

}

template<class ScalarType>
CSysVector<ScalarType>::~CSysVector() {
  delete [] vec_val;

  nElm = 0; nElmDomain = 0;
  nBlk = 0; nBlkDomain = 0;
  nVar = 0;

}

template<class ScalarType>
void CSysVector<ScalarType>::Initialize(const unsigned long & numBlk, const unsigned long & numBlkDomain, const unsigned short & numVar, const ScalarType & val) {

  nElm = numBlk*numVar; nElmDomain = numBlkDomain*numVar;
  nBlk = numBlk; nBlkDomain = numBlkDomain;
  nVar = numVar;

  /*--- Check for invalid input, then allocate memory and initialize values ---*/
  if ( nElm >= ULONG_MAX ) {
    char buf[100];
    SPRINTF(buf, "invalid input: numBlk, numVar = %lu, %u", numBlk, numVar );
    SU2_MPI::Error(string(buf), CURRENT_FUNCTION);
  }

  vec_val = new ScalarType[nElm];
  for (unsigned long i = 0; i < nElm; i++)
    vec_val[i] = val;

#ifdef HAVE_MPI
  unsigned long nElmLocal = (unsigned long)nElm;
  SU2_MPI::Allreduce(&nElmLocal, &nElmGlobal, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif

}

template<class ScalarType>
void CSysVector<ScalarType>::Equals_AX(const ScalarType & a, CSysVector<ScalarType> & x) {
  /*--- check that *this and x are compatible ---*/
  if (nElm != x.nElm) {
    cerr << "CSysVector::Equals_AX(): " << "sizes do not match";
    throw(-1);
  }
  for (unsigned long i = 0; i < nElm; i++)
    vec_val[i] = a * x.vec_val[i];
}

template<class ScalarType>
void CSysVector<ScalarType>::Plus_AX(const ScalarType & a, CSysVector<ScalarType> & x) {
  /*--- check that *this and x are compatible ---*/
  if (nElm != x.nElm) {
    SU2_MPI::Error("Sizes do not match", CURRENT_FUNCTION);
  }
  for (unsigned long i = 0; i < nElm; i++)
    vec_val[i] += a * x.vec_val[i];
}

template<class ScalarType>
void CSysVector<ScalarType>::Equals_AX_Plus_BY(const ScalarType & a, CSysVector<ScalarType> & x, const ScalarType & b, CSysVector<ScalarType> & y) {
  /*--- check that *this, x and y are compatible ---*/
  if ((nElm != x.nElm) || (nElm != y.nElm)) {
    SU2_MPI::Error("Sizes do not match", CURRENT_FUNCTION);
  }
  for (unsigned long i = 0; i < nElm; i++)
    vec_val[i] = a * x.vec_val[i] + b * y.vec_val[i];
}

template<class ScalarType>
CSysVector<ScalarType> & CSysVector<ScalarType>::operator=(const CSysVector<ScalarType> & u) {

  /*--- check if self-assignment, otherwise perform deep copy ---*/
  if (this == &u) return *this;

  /*--- determine if (re-)allocation is needed ---*/
  if (nElm != u.nElm && vec_val != NULL) {delete [] vec_val; vec_val = NULL;}
  if (vec_val == NULL) vec_val = new ScalarType[u.nElm];

  /*--- copy ---*/
  nElm = u.nElm;
  nElmDomain = u.nElmDomain;
  nBlk = u.nBlk;
  nBlkDomain = u.nBlkDomain;
  nVar = u.nVar;

  for (unsigned long i = 0; i < nElm; i++)
    vec_val[i] = u.vec_val[i];

#ifdef HAVE_MPI
  nElmGlobal = u.nElmGlobal;
#endif

  return *this;
}

template<class ScalarType>
CSysVector<ScalarType> & CSysVector<ScalarType>::operator=(const ScalarType & val) {
  for (unsigned long i = 0; i < nElm; i++)
    vec_val[i] = val;
  return *this;
}

template<class ScalarType>
CSysVector<ScalarType> CSysVector<ScalarType>::operator+(const CSysVector<ScalarType> & u) const {

  /*--- Use copy constructor and compound addition-assignment ---*/
  CSysVector<ScalarType> sum(*this);
  sum += u;
  return sum;
}

template<class ScalarType>
CSysVector<ScalarType> & CSysVector<ScalarType>::operator+=(const CSysVector<ScalarType> & u) {

  /*--- Check for consistent sizes, then add elements ---*/
  if (nElm != u.nElm) {
    SU2_MPI::Error("Sizes do not match", CURRENT_FUNCTION);
  }
  for (unsigned long i = 0; i < nElm; i++)
    vec_val[i] += u.vec_val[i];
  return *this;
}

template<class ScalarType>
CSysVector<ScalarType> CSysVector<ScalarType>::operator-(const CSysVector<ScalarType> & u) const {

  /*--- Use copy constructor and compound subtraction-assignment ---*/
  CSysVector<ScalarType> diff(*this);
  diff -= u;
  return diff;
}

template<class ScalarType>
CSysVector<ScalarType> & CSysVector<ScalarType>::operator-=(const CSysVector<ScalarType> & u) {

  /*--- Check for consistent sizes, then subtract elements ---*/
  if (nElm != u.nElm) {
    SU2_MPI::Error("Sizes do not match", CURRENT_FUNCTION);
  }
  for (unsigned long i = 0; i < nElm; i++)
    vec_val[i] -= u.vec_val[i];
  return *this;
}

template<class ScalarType>
CSysVector<ScalarType> CSysVector<ScalarType>::operator*(const ScalarType & val) const {

  /*--- use copy constructor and compound scalar
   multiplication-assignment ---*/
  CSysVector<ScalarType> prod(*this);
  prod *= val;
  return prod;
}

template<class ScalarType>
CSysVector<ScalarType> operator*(const ScalarType & val, const CSysVector<ScalarType> & u) {

  /*--- use copy constructor and compound scalar
   multiplication-assignment ---*/
  CSysVector<ScalarType> prod(u);
  prod *= val;
  return prod;
}

template<class ScalarType>
CSysVector<ScalarType> & CSysVector<ScalarType>::operator*=(const ScalarType & val) {

  for (unsigned long i = 0; i < nElm; i++)
    vec_val[i] *= val;
  return *this;
}

template<class ScalarType>
CSysVector<ScalarType> CSysVector<ScalarType>::operator/(const ScalarType & val) const {

  /*--- use copy constructor and compound scalar
   division-assignment ---*/
  CSysVector quotient(*this);
  quotient /= val;
  return quotient;
}

template<class ScalarType>
CSysVector<ScalarType> & CSysVector<ScalarType>::operator/=(const ScalarType & val) {

  for (unsigned long i = 0; i < nElm; i++)
    vec_val[i] /= val;
  return *this;
}

template<class ScalarType>
ScalarType CSysVector<ScalarType>::norm() const {

  /*--- just call dotProd on this*, then sqrt ---*/
  ScalarType val = dotProd(*this, *this);
  if (val < 0.0) {
    SU2_MPI::Error("Inner product of CSysVector is negative", CURRENT_FUNCTION);
  }
  return sqrt(val);
}

template<class ScalarType>
void CSysVector<ScalarType>::CopyToArray(ScalarType* u_array) {

  for (unsigned long i = 0; i < nElm; i++)
    u_array[i] = vec_val[i];
}

template<class ScalarType>
void CSysVector<ScalarType>::AddBlock(unsigned long val_ipoint, ScalarType *val_residual) {
  unsigned short iVar;

  for (iVar = 0; iVar < nVar; iVar++)
    vec_val[val_ipoint*nVar+iVar] += val_residual[iVar];
}

template<class ScalarType>
void CSysVector<ScalarType>::SubtractBlock(unsigned long val_ipoint, ScalarType *val_residual) {
  unsigned short iVar;

  for (iVar = 0; iVar < nVar; iVar++)
    vec_val[val_ipoint*nVar+iVar] -= val_residual[iVar];
}

template<class ScalarType>
void CSysVector<ScalarType>::SetBlock(unsigned long val_ipoint, ScalarType *val_residual) {
  unsigned short iVar;

  for (iVar = 0; iVar < nVar; iVar++)
    vec_val[val_ipoint*nVar+iVar] = val_residual[iVar];
}

template<class ScalarType>
void CSysVector<ScalarType>::SetBlock(unsigned long val_ipoint, unsigned short val_var, ScalarType val_residual) {

  vec_val[val_ipoint*nVar+val_var] = val_residual;
}

template<class ScalarType>
void CSysVector<ScalarType>::SetBlock_Zero(unsigned long val_ipoint) {
  unsigned short iVar;

  for (iVar = 0; iVar < nVar; iVar++)
    vec_val[val_ipoint*nVar+iVar] = 0.0;
}

template<class ScalarType>
void CSysVector<ScalarType>::SetBlock_Zero(unsigned long val_ipoint, unsigned short val_var) {
    vec_val[val_ipoint*nVar+val_var] = 0.0;
}

template<class ScalarType>
ScalarType CSysVector<ScalarType>::GetBlock(unsigned long val_ipoint, unsigned short val_var) {
  return vec_val[val_ipoint*nVar + val_var];
}

template<class ScalarType>
ScalarType *CSysVector<ScalarType>::GetBlock(unsigned long val_ipoint) {
  return &vec_val[val_ipoint*nVar];
}

template<class ScalarType>
template<class T>
void CSysVector<ScalarType>::PassiveCopy(const CSysVector<T>& other) {

  /*--- This is a method and not the overload of an operator to make sure who
   calls it knows the consequence to the derivative information (lost) ---*/

  /*--- check if self-assignment, otherwise perform deep copy ---*/
  if ((const void*)this == (const void*)&other) return;

  /*--- determine if (re-)allocation is needed ---*/
  if (nElm != other.GetLocSize() && vec_val != NULL) {
    delete [] vec_val;
    vec_val = NULL;
  }

  /*--- copy ---*/
  nElm = other.GetLocSize();
  nElmDomain = other.GetNElmDomain();
  nBlk = other.GetNBlk();
  nBlkDomain = other.GetNBlkDomain();
  nVar = other.GetNVar();

  if (vec_val == NULL)
    vec_val = new ScalarType[nElm];

  for (unsigned long i = 0; i < nElm; i++)
    vec_val[i] = SU2_TYPE::GetValue(other[i]);

#ifdef HAVE_MPI
  nElmGlobal = other.GetSize();
#endif
}

template<class ScalarType>
ScalarType dotProd(const CSysVector<ScalarType> & u, const CSysVector<ScalarType> & v) {

  /*--- check for consistent sizes ---*/
  if (u.nElm != v.nElm) {
    SU2_MPI::Error("Sizes do not match", CURRENT_FUNCTION);
  }

  /*--- find local inner product and, if a parallel run, sum over all
   processors (we use nElemDomain instead of nElem) ---*/
  ScalarType loc_prod = 0.0;
  for (unsigned long i = 0; i < u.nElmDomain; i++)
    loc_prod += u.vec_val[i]*v.vec_val[i];
  ScalarType prod = 0.0;

#ifdef HAVE_MPI
  SelectMPIWrapper<ScalarType>::W::Allreduce(&loc_prod, &prod, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
  prod = loc_prod;
#endif

  return prod;
}

/*--- Explicit instantiations ---*/
template class CSysVector<su2double>;
template CSysVector<su2double> operator*(const su2double&, const CSysVector<su2double>&);
template void CSysVector<su2double>::PassiveCopy(const CSysVector<su2double>&);
template su2double dotProd<su2double>(const CSysVector<su2double> & u, const CSysVector<su2double> & v);

template class CSysVector<unsigned long>;

#ifdef CODI_REVERSE_TYPE
template class CSysVector<passivedouble>;
template CSysVector<passivedouble> operator*(const passivedouble&, const CSysVector<passivedouble>&);
template void CSysVector<su2double>::PassiveCopy(const CSysVector<passivedouble>&);
template void CSysVector<passivedouble>::PassiveCopy(const CSysVector<su2double>&);
template passivedouble dotProd<passivedouble>(const CSysVector<passivedouble> & u, const CSysVector<passivedouble> & v);
#endif
