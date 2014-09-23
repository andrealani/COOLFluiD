// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "SAMGLSS/SAMGLSSMatrix.hh"
#include <fstream>
#include "Common/CFLog.hh"
#include "Framework/BlockAccumulator.hh"
#include "Framework/LSSVector.hh"
#include "Common/BadValueException.hh"

namespace COOLFluiD {
  namespace SAMGLSS {

//////////////////////////////////////////////////////////////////////////////

void SAMGLSSMatrix::createJAStructure(
  std::vector< std::vector< CFuint > >& nz )
{
  CFAUTOTRACE;

  cf_assert( m_nnu && m_nna && m_nsys && (m_ia!=CFNULL) && (m_ja!=CFNULL) );
  m_preconstruct = true;

  // construct ja off-diagonal entries (r/c are 0-based row/column indices)
  for (int r=0;r<m_nnu;++r) {  // row index (skipping diagonal entry)
    int j = m_ia[r];
    for (CFuint n=0;n<nz[r].size();++n) {          // for each neighbor
      const int c = static_cast< int >(nz[r][n]);  // column index, 0-based
      if (c!=r)
        m_ja[j++] = c+1;
    }  // for each neigbor (index n)
  }  // for each row (r)
}

//////////////////////////////////////////////////////////////////////////////

void SAMGLSSMatrix::createSeqAIJ(
  const CFint m, const CFint n, const CFint nz, const CFint* nnz,
  const char* name)
{
  CFAUTOTRACE;
  createSeqBAIJ(1,m,n,nz,nnz,name);
}

//////////////////////////////////////////////////////////////////////////////

void SAMGLSSMatrix::createSeqBAIJ(
  const CFuint blockSize, const CFint m, const CFint n, const CFint nz,
  const CFint* nnz, const char* name)
{
  CFAUTOTRACE;

  m_name = std::string(name);  // why not?
  m_preconstruct = false;
  m_nnu = m;
  m_nnu_nrhalo = n;
  m_nsys = static_cast< int >(blockSize);

  // matrix row indices construction
  m_ia = new int[m_nnu+1];
  m_ia[0] = 1;
  for (int r=0;r<m;++r)
    m_ia[r+1] = m_ia[r]+(nz? nz:nnz[r]);

  // number of non-zero elements
  m_nna = m_ia[m_nnu]-1;

  /* fill matrix diagonal indices, where r/c are absolute row/column indices,
     0-based, first adding diagonal then clear the row */
  m_ja = new int[m_nna];
  for (int r=0;r<m_nnu;++r) {
    int j = m_ia[r]-1;
    m_ja[j++] = r+1;
    for (;j<m_ia[r+1]-1;++j)
      m_ja[j]=0;
  }

  // allocate memory for matrix values and reset
  m_a = new double[m_nna];
  resetToZeroEntries();
}

//////////////////////////////////////////////////////////////////////////////

#ifdef CF_HAVE_MPI
void SAMGLSSMatrix::createParAIJ(
  MPI_Comm comm,
  const CFint m, const CFint n, const CFint M, const CFint N,
  const CFint dnz, const int* dnnz, const CFint onz, const int* onnz,
  const char* name)
{
  CFAUTOTRACE;
  createParBAIJ(comm,1,m,n,M,N,dnz,dnnz,onz,onnz,name);
}
#endif

//////////////////////////////////////////////////////////////////////////////

#ifdef CF_HAVE_MPI
void SAMGLSSMatrix::createParBAIJ(
  MPI_Comm comm, const CFuint blockSize,
  const CFint m, const CFint n, const CFint M, const CFint N,
  const CFint dnz, const int* dnnz, const CFint onz, const int* onnz,
  const char* name)
{
  CFAUTOTRACE;
  const int nz = (dnz+onz)*static_cast< int >(blockSize);
  int *nnz = new int[m];
  for (int i=0;i<m;++i)
    nnz[i] = (dnnz[i]+onnz[i])*static_cast< int >(blockSize);

  createSeqBAIJ(blockSize,m,n,nz,nnz,name);

  delete[] nnz;
}
#endif

//////////////////////////////////////////////////////////////////////////////

void SAMGLSSMatrix::printToScreen() const
{
  CFout << "SAMGLSS: SAMGLSSMatrix\n";

  // cycle all rows
  std::vector< int > idx;
  for (int r=0;r<m_nnu;++r) {

    // push column indices and sort descending
    idx.reserve(m_ia[r+1]-m_ia[r]);
    for (int i=m_ia[r]-1;i<m_ia[r+1]-1;++i)
      idx.push_back(m_ja[i]-1);
    sort(idx.begin(),idx.end(),std::greater< int >());
    idx.erase(unique(idx.begin(),idx.end()),idx.end());

    // pop column indices in ascending order
    int c=0;
    while (!idx.empty()) {
      const int c_next = idx.back();
      idx.pop_back();
      if (c_next<0) continue;
      for (;c<c_next;++c)
        CFout << "~ ";
      // seek and output value at (r,c_next)
      for (int i=m_ia[r]-1;i<m_ia[r+1]-1;++i)
        if (c_next==m_ja[i]-1) {
          CFout << m_a[i] << " ";
          ++c;
          break;
        }
    }
    for (;c<m_nnu_nrhalo;++c)
      CFout << "~ ";
    CFout << "\n";
  } // cycle all rows

  CFout << "\n";
}

//////////////////////////////////////////////////////////////////////////////

void SAMGLSSMatrix::printToFile(const char* fileName) const
{
  std::string file(fileName);
  std::ofstream ofile(file.c_str());
  ofile.precision(16);

  // cycle all rows
  std::vector< int > idx;
  for (int r=0;r<m_nnu;++r) {

    // push column indices and sort descending, keeping unique indices
    idx.reserve(m_ia[r+1]-m_ia[r]);
    for (int i=m_ia[r]-1;i<m_ia[r+1]-1;++i)
      idx.push_back(m_ja[i]-1);
    sort(idx.begin(),idx.end(),std::greater< int >());
    idx.erase(unique(idx.begin(),idx.end()),idx.end());

    // pop 0-based column indices in ascending order
    int c=0;
    while (!idx.empty()) {
      const int c_next = idx.back();
      idx.pop_back();
      if (c_next<0)  continue;
      for (;c<c_next;++c)
        ofile << "~ ";
      // locate and output value at (r,c_next)
      for (int i=m_ia[r]-1;i<m_ia[r+1]-1;++i)
        if (c_next==m_ja[i]-1) {
          ofile << m_a[i] << " ";
          ++c;
          break;
        }
    }
    for (;c<m_nnu_nrhalo;++c)
      ofile << "~ ";
    ofile << std::endl;
  } // cycle all rows

  ofile << std::endl;
  ofile.close();

  // .amg file (CSR structure)
  file = std::string(fileName) + ".amg";
  ofile.open(file.c_str());
  ofile.precision(16);
  ofile << "ia" << std::endl;
  CFint* pia = m_ia;
  for (CFint i=0;i<=m_nnu;++i)
    ofile << *pia++ << std::endl;
  ofile << "ja" << std::endl;
  CFint* pja = m_ja;
  for (CFint i=0;i<m_nna;++i)
    ofile << *pja++  << std::endl;
  ofile << "a" << std::endl;
  CFreal* pa = m_a;
  for (CFint i=0;i<m_nna;++i)
    ofile << *pa++ << std::endl;
  ofile.close();

  // .frm file with fixed parameters:
  // f  data is formatted
  // 4  format specification is version 4 (...)
  // 0  (or 1) if optional ip file exists
  file = std::string(fileName) + ".frm";
  ofile.open(file.c_str());
  ofile << "f  4" << std::endl
    << "  " << m_nna << "  " << m_nnu << "  " << m_matrix << "  " << m_nsys
    << "  0" << std::endl;
  ofile.close();
}

//////////////////////////////////////////////////////////////////////////////

void SAMGLSSMatrix::setValue(
  const CFint im, const CFint in, const CFreal value )
{
  if (im>=m_nnu || in>=m_nnu_nrhalo) return;
  if (im!=in) {
    // off-diagonal element (check allocation)
    const int ki = m_ia[im];
    const int kf = m_ia[im+1]-1;
    if (kf-ki>0)
      for (int k=ki;k<kf;++k)
        if (in==m_ja[k]-1) {
          m_a[k] = value;
          return;
        }
  } else {
    // diagonal element
    m_a[ m_ia[im]-1 ] = value;
    return;
  }
  notFoundIndex(im,in);
  setValue(im,in,value);
}

//////////////////////////////////////////////////////////////////////////////

void SAMGLSSMatrix::setValues(
  const CFuint m, const CFint* im, const CFuint n, const CFint* in,
  const CFreal* values )
{
  const int mn = (int) (m*n);
  int* ik = new int[mn];
  getCSRIndices(m,im,n,in, ik);
  for (int i=0;i<mn;++i)
    if (ik[i]>=0)
      m_a[ik[i]] = values[i];
  delete[] ik;
}

//////////////////////////////////////////////////////////////////////////////

void SAMGLSSMatrix::addValue(const CFint im,const CFint in,const CFreal value)
{
  if (im>=m_nnu||in>=m_nnu_nrhalo) return;
  if (im!=in) {
    // off-diagonal element
    for (int k=m_ia[im];k<m_ia[im+1]-1;++k) {
      if (in==m_ja[k]-1) {
        m_a[k] += value;
        return;
      }
    }
  } else {
    // diagonal element
    m_a[ m_ia[im]-1 ] += value;
    return;
  }
  notFoundIndex(im,in);
  addValue(im,in,value);
}

//////////////////////////////////////////////////////////////////////////////

void SAMGLSSMatrix::addValues(
  const CFuint m, const CFint* im, const CFuint n, const CFint* in,
  const CFreal* values )
{
  const int mn = (int) (m*n);
  int* ik = new int[mn];
  getCSRIndices(m,im,n,in, ik);
  for (int i=0;i<mn;++i)
    if (ik[i]>=0)
      m_a[ik[i]] += values[i];
  delete[] ik;
}

//////////////////////////////////////////////////////////////////////////////

void SAMGLSSMatrix::getValue(
  const CFint im, const CFint in, CFreal& value )
{
  value = getValue(im,in);
}

//////////////////////////////////////////////////////////////////////////////

void SAMGLSSMatrix::getValues(
  const CFuint m, const CFint* im, const CFuint n, const CFint* in,
  CFreal* values )
{
  const int mn = (int) (m*n);
  int* ik = new int[mn];
  getCSRIndices(m,im,n,in, ik);
  for (int i=0;i<mn;++i)
    values[i] = (ik[i]<0? 0.:m_a[ik[i]]);
  delete[] ik;
}

//////////////////////////////////////////////////////////////////////////////

void SAMGLSSMatrix::setDiagonal(Framework::LSSVector& diag)
{
  CFreal value=0.;
  for (int i=0;i<m_nnu;++i) {
    diag.getValue(i,value);
    m_a[m_ia[i]-1] = value;
  }
}

//////////////////////////////////////////////////////////////////////////////

void SAMGLSSMatrix::addToDiagonal(Framework::LSSVector& diag)
{
  CFreal value=0.;
  for (int i=0;i<m_nnu;++i) {
    diag.getValue(i,value);
    m_a[m_ia[i]-1] += value;
  }
}

//////////////////////////////////////////////////////////////////////////////

void SAMGLSSMatrix::resetToZeroEntries()
{
  for (int i=0;i<m_nna;++i)
    m_a[i] = 0.;
}

//////////////////////////////////////////////////////////////////////////////

void SAMGLSSMatrix::setValues(const Framework::BlockAccumulator& acc)
{
  const int k = acc.getM()*acc.getNB()*acc.getN()*acc.getNB();
  int* ik = new int[k];
  getBlockAccumulatorCSRIndices(acc, ik);
  for (int i=0;i<k;++i)
    m_a[ik[i]] = acc.values[i];
  delete[] ik;
}

//////////////////////////////////////////////////////////////////////////////

void SAMGLSSMatrix::addValues(const Framework::BlockAccumulator& acc)
{
  const int k = acc.getM()*acc.getNB()*acc.getN()*acc.getNB();
  int* ik = new int[k];
  getBlockAccumulatorCSRIndices(acc, ik);
  for (int i=0;i<k;++i)
    m_a[ik[i]] += acc.values[i];
  delete[] ik;
}

//////////////////////////////////////////////////////////////////////////////

void SAMGLSSMatrix::freezeNonZeroStructure()
{
}

//////////////////////////////////////////////////////////////////////////////

// Non-abstract functions

//////////////////////////////////////////////////////////////////////////////

void SAMGLSSMatrix::getCSRIndices(
  const CFuint m, const CFint* im, const CFuint n, const CFint* in, int*& ik )
{
  // this routine assumes all rows of the matrix exist
  int iki=0;

  // cycle given rows
  for (CFuint imi=0;imi<m;++imi) {
    const int r = im[imi];

    // cycle given columns
    for (CFuint ini=0;ini<n;++ini) {
      const int c = in[ini];
      if (r<0||c<0||r>=m_nnu||c>=m_nnu_nrhalo) {
        ik[iki++] = -1;
        continue;
      }

      // seek column indices to see if c is present in the row ranging from
      // m_ia[r]-1 and m_ia[r+1]-1 (non-inclusively), r==c breaks immediately
      bool foundIndex = false;
      for (int seek=m_ia[r]-1;seek<m_ia[r+1]-1;++seek) {
        if (!m_ja[seek])  break;
        if (m_ja[seek]-1!=c)  continue;
        ik[iki++] = seek;
        foundIndex = true;
        break;
      }
      if (!foundIndex) {
        /* either break execution or recursively repeat until all indices are
           allocated, given enough free space */
        notFoundIndex(r,c);
        getCSRIndices(m,im,n,in, ik);
        return;
      }

    } // cycle given columns
  } // cycle given rows
}

//////////////////////////////////////////////////////////////////////////////

void SAMGLSSMatrix::getBlockAccumulatorCSRIndices(
  const Framework::BlockAccumulator& acc, int*& ik )
{
  // this routine assumes all rows of the matrix exist
  const CFuint nbEqs = acc.getNB();
  int iki=0;

  // cycle given rows
  for (CFuint imi=0; imi<acc.getM(); ++imi) { for (CFuint er=0; er<nbEqs; ++er) {
    const int r = (int) acc.getIM()[imi]*nbEqs+er;

    // cycle given columns
    for (CFuint ini=0;ini<acc.getN();++ini) { for (CFuint ec=0;ec<nbEqs;++ec) {
      const int c = (int) acc.getIN()[ini]*nbEqs+ec;

      if (r<0||c<0||r>=m_nnu||c>=m_nnu_nrhalo) {
        ik[iki++] = -1;
        continue;
      }

      // seek column indices to see if c is present in the row ranging from
      // m_ia[r]-1 and m_ia[r+1]-1 (non-inclusively), r==c breaks immediately
      bool foundIndex = false;
      for (int seek=m_ia[r]-1;seek<m_ia[r+1]-1;++seek) {
        if (!m_ja[seek])  break;
        if (m_ja[seek]-1!=c)  continue;
        ik[iki++] = seek;
        foundIndex = true;
        break;
      }
// std::string x = "looking for r,c = " +StringOps::to_str(r)+" , "+StringOps::to_str(c) + (!foundIndex?" not found!\n":"\n"); CFLog(INFO,x);
      if (!foundIndex && m_preconstruct)
        ik[iki++] = -1;
      else if (!foundIndex) {
        /* either break execution or recursively repeat until all indices are
           allocated, given enough free space */
        notFoundIndex(r,c);
        getBlockAccumulatorCSRIndices(acc,ik);
        return;
      }

    } }  // cycle given columns
  } }  // cycle given rows
}


//////////////////////////////////////////////////////////////////////////////
/*
void SAMGLSSMatrix::getBlockAccumulatorCSRIndices(
  const Framework::BlockAccumulator& acc,int*& ik )
{
  const CFint nb = static_cast< CFint >(acc.getNB());
  const CFuint m = acc.getM()*acc.getNB();
  CFint *im = new CFint[m];
  for (CFint i=0;i<acc.getM();++i) {
    const CFint acc_im_i = acc.getIM()[i];
    for (CFint e=0;e<nb;++e)
      im[i*nb+e] = (acc_im_i<0? -1:acc_im_i*nb+e);
  }

  const CFuint n = acc.getN()*acc.getNB();
  CFint *in = new CFint[n];
  for (CFint i=0;i<acc.getN();++i)
    for (CFint e=0;e<nb;++e)
      in[i*nb+e] = acc.getIN()[i]*nb+e;

  getCSRIndices(m,im,n,in, ik);
}
*/
//////////////////////////////////////////////////////////////////////////////

double SAMGLSSMatrix::getValue(
  const CFint im, const CFint in ) const
{
  cf_assert(std::max<CFint>(im,in)<m_nnu);
  if (im!=in) {
    // off-diagonal element
    for (int k=m_ia[im];k<m_ia[im+1]-1;++k)
      if (in==m_ja[k]-1)
        return m_a[k];
    return 0.;
  } else
    return getDiagonalValue(im);
}

//////////////////////////////////////////////////////////////////////////////

void SAMGLSSMatrix::notFoundIndex(const int r, const int c)
{
  CFAUTOTRACE;
  if (!m_preconstruct) {
    // allocate off-diagonal entry (r,c) where space is found in m_ja
    for (int j=m_ia[r];j<m_ia[r+1]-1;++j)
      if (!m_ja[j]) {
        m_ja[j]= c+1;
        return;
      }
  }
  CFLog(INFO,"SAMGLSS: matrix index not allocated ("<< r <<","<< c <<")\n");
  throw Common::BadValueException (FromHere(),"SAMGLSS: matrix index not allocated\n");
}

//////////////////////////////////////////////////////////////////////////////
/*
int SAMGLSSMatrix::findFreeSpaceOnRow(const int nr) const
{
  // assuming diagonal is allocated, starting from last element of nr row
  int free = 0;
  for (int j=m_ia[nr+1]-2;j>=m_ia[nr]+m_nsys;--j)
    if (!m_ja[j])
      ++free;
    else
      break;
  return free;
}
*/
//////////////////////////////////////////////////////////////////////////////
/*
void SAMGLSSMatrix::moveRows(const int r1, const int r2, const int off)
{
  // move rows [r1,r2] inclusive by offset, positive or negative, and new
  // positions are zeroed (this is a destructive operation and invalidates
  // current CSR indexes, use with care)
  const int r3 = r2+1;
  if (off>0) {
    for (int j=m_ia[r3]-2-off;j>=m_ia[r1]-1;--j) {
      m_ja[j+off] = m_ja[j];
      m_a[j+off]  = m_a[j];
    }
    for (int j=m_ia[r1]-1;j<m_ia[r1]-1+off;++j) {
      m_ja[j] = 0;
      m_a[j]  = 0.;
    }
    for (int r=r2;r>=r1;--r)
      m_ia[r]+=off;
  } else {
    for (int j=m_ia[r1]-1;j<=m_ia[r3]-2-off;++j) {
      m_ja[j+off] = m_ja[j];
      m_a[j+off]  = m_a[j];
    }
    for (int j=m_ia[r3]-2;j>m_ia[r3]-2+off;--j) {
      m_ja[j] = 0;
      m_a[j]  = 0.;
    }
    for (int r=r2;r>=r1;--r)
      m_ia[r]+=off;
  }
}
*/
//////////////////////////////////////////////////////////////////////////////

void SAMGLSSMatrix::compressStructure()
{
  //FIXME implement (because of n-coupled systems with <n bc's!)
}

//////////////////////////////////////////////////////////////////////////////

void SAMGLSSMatrix::getMatrix( int& nnu, int& nna, int& nsys, int& matrix,
  int*& ia, int*& ja, double*& a ) const
{
  nnu    = m_nnu;
  nna    = m_nna;
  nsys   = m_nsys;
  matrix = m_matrix;
  ia     = m_ia;
  ja     = m_ja;
  a      = m_a;
}

//////////////////////////////////////////////////////////////////////////////

void SAMGLSSMatrix::setMatrix( int nnu, int nna, int nsys, int matrix,
  int* ia, int* ja, double* a )
{
  if (m_nnu) {
    delete[] m_ia;
    delete[] m_ja;
    delete[] m_a;
  }
  m_nnu    = nnu;
  m_nna    = nna;
  m_nsys   = nsys;
  m_matrix = matrix;

  m_ia = new int[m_nnu+1];
  m_ja = new int[m_nna];
  m_a  = new double[m_nna];

  for (int i=0;i<m_nnu+1;++i)
    m_ia[i] = ia[i];
  for (int i=0;i<m_nna;++i)
    m_ja[i] = ja[i];
  for (int i=0;i<m_nna;++i)
    m_a[i] = a[i];
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SAMGLSS
}  // namespace COOLFluiD

