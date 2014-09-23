// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#include "Common/StringOps.hh"
#include "Common/CFLog.hh"
#include "Common/BadValueException.hh"
#include "Framework/BlockAccumulator.hh"
#include "Framework/LSSVector.hh"
#include "Pardiso/PardisoMatrix.hh"

namespace COOLFluiD {
  namespace Pardiso {

//////////////////////////////////////////////////////////////////////////////

void PardisoMatrix::createSeqAIJ(
  const CFint m, const CFint n, const CFint nz, const CFint* nnz,
  const char* name)
{
  createSeqBAIJ(1,m,n,nz,nnz,name);
}

//////////////////////////////////////////////////////////////////////////////

void PardisoMatrix::createSeqBAIJ(
  const CFuint blockSize, const CFint m, const CFint n, const CFint nz,
  const CFint* nnz, const char* name)
{
  CFAUTOTRACE;

  cf_assert(m==n);
  m_name = std::string(name);  // why not?
  m_nnu = m;
  m_nnu_nrhalo = n;
  m_nsys = static_cast< int >(blockSize);

  // matrix row indices construction
  m_ia = new int[m_nnu+1];
  m_ia[0] = 1;
  for (int r=0; r<m; ++r)
    m_ia[r+1] = m_ia[r]+(nz? nz:nnz[r]);

  // number of non-zero elements
  m_nna = m_ia[m_nnu]-1;

  // matrix column indices allocation
  m_ja = new int[m_nna];
  for (int i=0; i<m_nna; ++i)
    m_ja[i] = 0;

  // allocate memory for matrix values and reset
  m_a = new double[m_nna];
  resetToZeroEntries();
}

//////////////////////////////////////////////////////////////////////////////

#ifdef CF_HAVE_MPI
void PardisoMatrix::createParAIJ(
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
void PardisoMatrix::createParBAIJ(
  MPI_Comm comm, const CFuint blockSize,
  const CFint m, const CFint n, const CFint M, const CFint N,
  const CFint dnz, const int* dnnz, const CFint onz, const int* onnz,
  const char* name)
{
  CFAUTOTRACE;
  const int nz = (dnz+onz)*static_cast< int >(blockSize);
  int *nnz = new int[m];
  for (int i=0; i<m; ++i)
    nnz[i] = (dnnz[i]+onnz[i])*static_cast< int >(blockSize);

  createSeqBAIJ(blockSize,m,n,nz,nnz,name);

  delete[] nnz;
}
#endif

//////////////////////////////////////////////////////////////////////////////

void PardisoMatrix::printToScreen() const
{
  // cycle all rows
  std::vector< int > idx;
  for (int r=0; r<m_nnu; ++r) {

    // push column indices and sort descending
    idx.reserve(m_ia[r+1]-m_ia[r]);
    for (int i=m_ia[r]-1; i<m_ia[r+1]-1; ++i)
      idx.push_back(m_ja[i]-1);
    sort(idx.begin(),idx.end(),std::greater< int >());
    idx.erase(unique(idx.begin(),idx.end()),idx.end());

    // pop column indices in ascending order
    int c=0;
    while (!idx.empty()) {
      const int c_next = idx.back();
      idx.pop_back();
      if (c_next<0) continue;
      for (; c<c_next; ++c)
        CFout << "~ ";
      // seek and output value at (r,c_next)
      for (int i=m_ia[r]-1; i<m_ia[r+1]-1; ++i)
        if (c_next==m_ja[i]-1) {
          CFout << m_a[i] << " ";
          ++c;
          break;
        }
    }
    for (; c<m_nnu_nrhalo; ++c)
      CFout << "~ ";
    CFout << "\n";
  } // cycle all rows

  CFout << "\n";
}

//////////////////////////////////////////////////////////////////////////////

void PardisoMatrix::printToFile(const char* fileName) const
{
  std::ofstream ofile;
  ofile.precision(16);
  std::string file;

#if 0
  // dense (as opposed to sparse) output
  file = std::string(fileName);
  ofile.open(file.c_str());

  // cycle all rows
  std::vector< int > idx;
  for (int r=0; r<m_nnu; ++r) {

    // push column indices and sort descending, keeping unique indices
    idx.reserve(m_ia[r+1]-m_ia[r]);
    for (int i=m_ia[r]-1; i<m_ia[r+1]-1; ++i)
      idx.push_back(m_ja[i]-1);
    sort(idx.begin(),idx.end(),std::greater< int >());
    idx.erase(unique(idx.begin(),idx.end()),idx.end());

    // pop 0-based column indices in ascending order
    int c=0;
    while (!idx.empty()) {
      const int c_next = idx.back();
      idx.pop_back();
      if (c_next<0)  continue;
      for (; c<c_next; ++c)
        ofile << "~ ";
      // locate and output value at (r,c_next)
      for (int i=m_ia[r]-1; i<m_ia[r+1]-1; ++i)
        if (c_next==m_ja[i]-1) {
          ofile << m_a[i] << " ";
          ++c;
          break;
        }
    }
    for (; c<m_nnu_nrhalo; ++c)
      ofile << "~ ";
    ofile << std::endl;
  }

  ofile << std::endl;
  ofile.close();
#endif

  // csr structure output
  file = std::string(fileName) + ".csr";
  ofile.open(file.c_str());
  ofile.precision(16);
  ofile << "ia" << std::endl;
  int *pia = m_ia;
  for (CFint i=0; i<=m_nnu; ++i, ++pia)
    ofile << *pia << " ";
  ofile << std::endl;
  ofile << "ja" << std::endl;
  int *pja = m_ja;
  for (CFint r=0; r<m_nnu; ++r) {
    for (CFint j=m_ia[r]; j<m_ia[r+1]; ++j, ++pja)
      ofile << *pja << " ";
    ofile << std::endl;
  }
  ofile << "a" << std::endl;
  double* pa = m_a;
  for (CFint r=0; r<m_nnu; ++r) {
    for (CFint j=m_ia[r]; j<m_ia[r+1]; ++j, ++pa)
      ofile << *pa << " ";
    ofile << std::endl;
  }
  ofile.close();
}

//////////////////////////////////////////////////////////////////////////////

void PardisoMatrix::setValue(
  const CFint im, const CFint in, const CFreal value )
{
  cf_assert(im<m_nnu && in<m_nnu_nrhalo);
  for (int k=m_ia[im]-1; k<m_ia[im+1]-1; ++k)
    if (in==m_ja[k]-1) {
      m_a[k] = value;
      return;
    }
  notFoundIndex(im,in);
}

//////////////////////////////////////////////////////////////////////////////

void PardisoMatrix::setValues(
  const CFuint m, const CFint* im, const CFuint n, const CFint* in,
  const CFreal* values )
{
  int* ik = new int[m*n];
  getCSRIndices(m,im,n,in, ik);
  for (CFuint i=0; i<m*n; ++i)
    if (ik[i]>=0)
      m_a[ik[i]] = values[i];
  delete[] ik;
}

//////////////////////////////////////////////////////////////////////////////

void PardisoMatrix::addValue(const CFint im,const CFint in,const CFreal value)
{
  cf_assert(im<m_nnu && in<m_nnu_nrhalo);
  for (int k=m_ia[im]-1; k<m_ia[im+1]-1; ++k)
    if (in==m_ja[k]-1) {
      m_a[k] += value;
      return;
    }
  notFoundIndex(im,in);
}

//////////////////////////////////////////////////////////////////////////////

void PardisoMatrix::addValues(
  const CFuint m, const CFint* im, const CFuint n, const CFint* in,
  const CFreal* values )
{
  int* ik = new int[m*n];
  getCSRIndices(m,im,n,in, ik);
  for (CFuint i=0; i<m*n; ++i)
    if (ik[i]>=0)
      m_a[ik[i]] += values[i];
  delete[] ik;
}

//////////////////////////////////////////////////////////////////////////////

void PardisoMatrix::getValue(
  const CFint im, const CFint in, CFreal& value )
{
  value = getValue(im,in);
}

//////////////////////////////////////////////////////////////////////////////

void PardisoMatrix::getValues(
  const CFuint m, const CFint* im, const CFuint n, const CFint* in,
  CFreal* values )
{
  int* ik = new int[m*n];

  getCSRIndices(m,im,n,in, ik);
  for (CFuint i=0; i<m*n; ++i)
    values[i] = (ik[i]<0? 0.:m_a[ik[i]]);
  delete[] ik;
}

//////////////////////////////////////////////////////////////////////////////

void PardisoMatrix::setDiagonal(Framework::LSSVector& diag)
{
  CFreal value=0.;
  for (int i=0; i<m_nnu; ++i) {
    diag.getValue(i,value);
    m_a[m_ia[i]-1] = value;
  }
}

//////////////////////////////////////////////////////////////////////////////

void PardisoMatrix::addToDiagonal(Framework::LSSVector& diag)
{
  CFreal value=0.;
  for (int i=0; i<m_nnu; ++i) {
    diag.getValue(i,value);
    m_a[m_ia[i]-1] += value;
  }
}

//////////////////////////////////////////////////////////////////////////////

void PardisoMatrix::resetToZeroEntries()
{
  for (int i=0; i<m_nna; ++i)
    m_a[i] = 0.;
}

//////////////////////////////////////////////////////////////////////////////

void PardisoMatrix::setValues(const Framework::BlockAccumulator& acc)
{
  const int k = acc.getM()*acc.getNB()*acc.getN()*acc.getNB();
  int* ik = new int[k];
  getBlockAccumulatorCSRIndices(acc, ik);
  for (int i=0; i<k; ++i)
    m_a[ik[i]] = acc.values[i];
  delete[] ik;
}

//////////////////////////////////////////////////////////////////////////////

void PardisoMatrix::addValues(const Framework::BlockAccumulator& acc)
{
  const int k = acc.getM()*acc.getNB()*acc.getN()*acc.getNB();
  int* ik = new int[k];
  getBlockAccumulatorCSRIndices(acc, ik);
  for (int i=0; i<k; ++i)
    m_a[ik[i]] += acc.values[i];
  delete[] ik;
}

//////////////////////////////////////////////////////////////////////////////

void PardisoMatrix::freezeNonZeroStructure()
{
}

//////////////////////////////////////////////////////////////////////////////

// Non-abstract functions

//////////////////////////////////////////////////////////////////////////////

void PardisoMatrix::createJAStructure(
  std::vector< std::vector< CFuint > >& nz )
{
  CFAUTOTRACE;

  cf_assert( m_nnu && m_nna && m_nsys && (m_ia!=CFNULL) && (m_ja!=CFNULL) );

  // construct ja entries (r/c are 0-based row/column indices)
  for (int r=0; r<m_nnu; ++r) {  // row index (skipping diagonal entry)
    int j = m_ia[r]-1;
    for (CFuint n=0; n<nz[r].size(); ++n) {        // for each neighbor
      const int c = static_cast< int >(nz[r][n]);  // column index, 0-based
      m_ja[j++] = c+1;
    }  // for each neigbor (index n)
  }  // for each row (r)
}

//////////////////////////////////////////////////////////////////////////////

void PardisoMatrix::getCSRIndices(
  const CFuint m, const CFint* im, const CFuint n, const CFint* in, int*& ik )
{
  // this routine assumes all rows of the matrix exist
  int iki=0;

  // cycle given rows
  for (CFuint imi=0; imi<m; ++imi) {
    const int r = im[imi];

    // cycle given columns
    for (CFuint ini=0; ini<n; ++ini) {
      const int c = in[ini];
      if (r<0||c<0||r>=m_nnu||c>=m_nnu_nrhalo) {
        ik[iki++] = -1;
        continue;
      }

      // seek column indices to see if c is present in the row r
      bool foundIndex = false;
      for (int k=m_ia[r]-1; k<m_ia[r+1]-1; ++k)
        if (m_ja[k]-1==c) {
          ik[iki++] = k;
          foundIndex = true;
          break;
        }
      if (!foundIndex)
        notFoundIndex(r,c);

    } // cycle given columns
  } // cycle given rows
}

//////////////////////////////////////////////////////////////////////////////

void PardisoMatrix::getBlockAccumulatorCSRIndices(
  const Framework::BlockAccumulator& acc, int*& ik )
{
  // this routine assumes all rows of the matrix exist
  const int nbEqs = (int) acc.getNB();
  int iki=0;

  // cycle given rows
  for (CFuint imi=0; imi<acc.getM(); ++imi) {
    for (int er=0; er<nbEqs; ++er) {
      const int r = acc.getIM()[imi]*nbEqs+er;

      // cycle given columns
      for (CFuint ini=0; ini<acc.getN(); ++ini) {
        for (int ec=0; ec<nbEqs; ++ec) {
          const int c = acc.getIN()[ini]*nbEqs+ec;

          if (r<0||c<0||r>=m_nnu||c>=m_nnu_nrhalo) {
            ik[iki++] = -1;
            continue;
          }

          // seek column indices to see if c is present in the row r
          bool foundIndex = false;
          for (int k=m_ia[r]-1; k<m_ia[r+1]-1; ++k)
            if (m_ja[k]-1==c) {
              ik[iki++] = k;
              foundIndex = true;
              break;
            }
          if (!foundIndex)
            notFoundIndex(r,c);

        }
      }  // cycle given columns

    }
  }  // cycle given rows
}


//////////////////////////////////////////////////////////////////////////////

double PardisoMatrix::getValue(
  const CFint im, const CFint in ) const
{
  cf_assert(std::max(im,in)<m_nnu);
  for (int k=m_ia[im]; k<m_ia[im+1]-1; ++k)
    if (in==m_ja[k]-1)
      return m_a[k];
  return 0.;
}

//////////////////////////////////////////////////////////////////////////////

void PardisoMatrix::notFoundIndex(const int r, const int c)
{
  std::string e("PardisoMatrix index (" + Common::StringOps::to_str(r) + "," +
    Common::StringOps::to_str(c) + ") not allocated");
  CFLog(INFO,e << "\n");
  throw Common::BadValueException (FromHere(),e +"\n");
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Pardiso
} // namespace COOLFluiD

