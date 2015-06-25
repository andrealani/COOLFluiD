// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.



#include "Common/PE.hh"
#include "TrilinosMatrix.hh"
#include "Common/NotImplementedException.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace Trilinos {

//////////////////////////////////////////////////////////////////////////////

TrilinosMatrix::TrilinosMatrix() :
  LSSMatrix(),
  _mat(NULL),
  _map(NULL),
  _toBeDestroyed(false),
  _isAssembled(false),
  _name("")
{
}

//////////////////////////////////////////////////////////////////////////////

TrilinosMatrix::TrilinosMatrix(const TrilinosMatrix& other) :
  LSSMatrix()
{
  cf_assert(other._toBeDestroyed);
  _map = other._map;
  _mat = new Epetra_CrsMatrix(*(other._mat));
  _toBeDestroyed = true;
}
      
//////////////////////////////////////////////////////////////////////////////

const TrilinosMatrix& TrilinosMatrix::operator= (const TrilinosMatrix& other)
{
  cf_assert(other._toBeDestroyed);
  destroy();
  _map = other._map;
  _mat = new Epetra_CrsMatrix(*(other._mat));
  _toBeDestroyed = true;
  return *this;
}

//////////////////////////////////////////////////////////////////////////////

TrilinosMatrix::~TrilinosMatrix()
{
  destroy();
}

//////////////////////////////////////////////////////////////////////////////

void TrilinosMatrix::destroy()
{
  if (_toBeDestroyed) {
    delete _mat;
  }
  _mat = NULL;
  _map = NULL;
  _toBeDestroyed = false;
  _isAssembled = false;
}

//////////////////////////////////////////////////////////////////////////////

void TrilinosMatrix::createSeqAIJ(const CFint m,
                               const CFint n,
                               const CFint nz,
                               const int* nnz,
                               const char* name)
{
  // checks
  cf_assert(!_toBeDestroyed);                 // otherwise, the user is doing something wrong
  cf_assert(_map!=NULL);
  cf_assert(_map->ElementSize() == 1);
  cf_assert(_map->NumGlobalElements() == m);  // size of the matrix and the map must be compatible
  cf_assert(m == n);

  _name = name;

  if (nnz == NULL) {
    _mat = new Epetra_CrsMatrix(Copy, *_map, (int)nz);
  }
  else {
    int *copy_nnz = new int[m];
    for (int i=0; i<m; i++) copy_nnz[i] = nnz[i];
    _mat = new Epetra_CrsMatrix(Copy, *_map, copy_nnz);
    delete[] copy_nnz;
    //_mat = new Epetra_CrsMatrix(Copy, *_map, const_cast<int*>(nnz));  // --> should be more efficient
    // REMARK: using staticProfile==true is dangerous: if the user underestimates nnz,
    //         Epetra can abort ...
  }

  _toBeDestroyed = true;
}

//////////////////////////////////////////////////////////////////////////////

void TrilinosMatrix::createSeqBAIJ(const CFuint blockSize,
                                const CFint m,
                                const CFint n,
                                const CFint nz,
                                const int* nnz,
                                const char* name)
{
  int nbRows = m;
  int nbBlockRows = static_cast<int>(m) / static_cast<int>(blockSize);  // integer division
  int nbCols = n;
  int new_nz = nz*blockSize;
  int *new_nnz;

  _name = name;

  if (nnz==NULL) new_nnz=NULL;
  else {
    new_nnz = new int[nbRows];
    for (int R=0; R<nbBlockRows; R++) {
      for (CFuint r=0; r<blockSize; r++) {
        new_nnz[R*blockSize+r] = nnz[R]*blockSize;
      }
    }
  }

  createSeqAIJ(nbRows, nbCols, new_nz, new_nnz, name);

  if (new_nnz!=NULL) delete[] new_nnz;
}

//////////////////////////////////////////////////////////////////////////////

#ifdef CF_HAVE_MPI
void TrilinosMatrix::createParAIJ(MPI_Comm comm,
                               const CFint m,
                               const CFint n,
                               const CFint M,
                               const CFint N,
                               const CFint dnz,
                               const int* dnnz,
                               const CFint onz,
                               const int* onnz,
                               const char* name)
{
  // checks
  cf_assert(!_toBeDestroyed);
  cf_assert(_map!=NULL);
  cf_assert(_map->ElementSize() == 1);
  cf_assert(_map->NumGlobalElements() == M);       // sizes of matrix and map must be compatible
  cf_assert(M == N);
  cf_assert(_map->NumMyElements() == m);           // sizes of matrix and map must be compatible
  cf_assert(m == n);

  _name = name;

  if (dnnz == NULL) {
    _mat = new Epetra_CrsMatrix(Copy, *_map, dnz + onz);
  }
  else {

    // translation of petsc style dnnz and onnz to epetra style nnz
    int *nnz = new int[m];
    for (int r=0; r<m; r++) nnz[r] = dnnz[r] + onnz[r];

    // creation of the matrix
    _mat = new Epetra_CrsMatrix(Copy, *_map, nnz);
    // REMARK: using staticProfile==true is dangerous: if the user does something wrong
    //         Epetra aborts ...

    // deallocation
    delete[] nnz;
  }

  _toBeDestroyed = true;
}
#endif

//////////////////////////////////////////////////////////////////////////////

#ifdef CF_HAVE_MPI
void TrilinosMatrix::createParBAIJ(MPI_Comm comm,
        const CFuint blockSize,
        const CFint m,
        const CFint n,
        const CFint M,
        const CFint N,
        const CFint dnz,
        const int* dnnz,
        const CFint onz,
        const int* onnz,
        const char* name)
{
  // checks
  cf_assert(!_toBeDestroyed);
  cf_assert(_map!=NULL);
  cf_assert(_map->ElementSize() == 1);
  cf_assert(_map->NumGlobalElements() == M);       // sizes of matrix and map must be compatible
  cf_assert(M == N);
  cf_assert(_map->NumMyElements() == m);           // sizes of matrix and map must be compatible
  cf_assert(n == m);
  cf_assert(M % blockSize == 0);
  cf_assert(m % blockSize == 0);

  _name = name;

  if (dnnz == NULL) {
    _mat = new Epetra_CrsMatrix(Copy, *_map, (dnz + onz) * blockSize);  // dnz, onz are block oriented --> * blockSize
  }
  else {

    // translation of petsc block style dnnz and onnz to epetra individual style nnz
    int *nnz = new int[m*blockSize];
    int counter = 0;
    for (int R=0; R<m; R++) {
      for (CFuint r=0; r<blockSize; r++) {
        nnz[counter] = (dnnz[R] + onnz[R]) * blockSize;
        counter++;                                      // TODO: move ++ into previous line
      }
    }

    // creation of the matrix
    _mat = new Epetra_CrsMatrix(Copy, *_map, nnz);
    // REMARK: using staticProfile==true is dangerous: if the user does something wrong
    //         Epetra aborts ...

    // deallocation
    delete[] nnz;
  }

  _toBeDestroyed = true;
}
#endif

//////////////////////////////////////////////////////////////////////////////

void TrilinosMatrix::printToFile(const char* fileName) const
{
  Common::PE::GetPE().GetCommunicator();

  throw Common::NotImplementedException (FromHere(),"TrilinosMatrix::printToFile()");
}

//////////////////////////////////////////////////////////////////////////////

void TrilinosMatrix::set1DPoisson(const Epetra_Map *map) {

  // checks
  cf_assert(!_toBeDestroyed);

#ifdef CF_HAVE_MPI

  _map = map;

  const std::string nsp = MeshDataStack::getActive()->getPrimaryNamespace();
  MPI_Comm comm = PE::GetPE().GetCommunicator(nsp);
  
  // communication
  int nbProcessors, myRank;
  MPI_Comm_size(comm, &nbProcessors);
  MPI_Comm_rank(comm, &myRank);
  
  int sendbuf[2] = {_map->GID(0), _map->GID(_map->NumMyElements()-1)};
  int *recvbuf = new int[2*nbProcessors];
  MPI_Allgather(sendbuf, 2, MPI_INT, recvbuf, 2, MPI_INT, comm);
  
  _mat = new Epetra_CrsMatrix(Copy, *map, 3);

  int     gcid[3], grid;
  double  val[3] = {-1.0, 2.0, -1.0};

  if (myRank==0) {
    grid    = _map->GID(0);
    gcid[0] = grid;
    gcid[1] = _map->GID(1);
    _mat->InsertGlobalValues(grid, 2, &val[1], gcid);
  }
  else {
    grid    = _map->GID(0);
    gcid[0] = recvbuf[2*myRank-1];
    gcid[1] = grid;
    gcid[2] = _map->GID(1);
    _mat->InsertGlobalValues(grid, 3, val, gcid);
  }

  for (int lrid=1; lrid<_map->NumMyElements()-1; lrid++) {
    grid    = _map->GID(lrid);
    gcid[0] = _map->GID(lrid-1);
    gcid[1] = grid;
    gcid[2] = _map->GID(lrid+1);
    _mat->InsertGlobalValues(grid, 3, val, gcid);
  }

  if (myRank==nbProcessors-1) {
    grid    = _map->GID(_map->NumMyElements()-1);
    gcid[0] = _map->GID(_map->NumMyElements()-2);
    gcid[1] = grid;
    _mat->InsertGlobalValues(grid, 2, val, gcid);
  }
  else {
    grid    = _map->GID(_map->NumMyElements()-1);
    gcid[0] = _map->GID(_map->NumMyElements()-2);
    gcid[1] = grid;
    gcid[2] = recvbuf[2*myRank+2];
    _mat->InsertGlobalValues(grid, 3, &val[0], gcid);
  }

  delete[] recvbuf;

#else // CF_HAVE_MPI

  int gcid[3], grid;
  double val[3] = {-1.0, 2.0, -1.0};

  // first row
  grid    = _map->GID(0);
  gcid[0] = grid;
  gcid[1] = _map->GID(1);
  _mat->InsertGlobalValues(grid, 2, &val[1], gcid);

  // internal rows
  for (int i=1; i<_map->NumMyElements(); i++) {
    grid    = _map->GID(i);
    gcid[0] = _map->GID(i - 1);
    gcid[1] = grid;
    gcid[2] = _map->GID(i + 1);
    _mat->InsertGlobalValues(grid, 3, val, gcid);
  }

  // last row
  grid    = _map->GID(_map->NumMyElements()-1);
  gcid[0] = _map->GID(_map->NumMyElements()-2);
  gcid[1] = grid;
  _mat->InsertGlobalValues(grid, 2, val, gcid);

#endif // CF_HAVE_MPI

  _mat->FillComplete();
  _mat->OptimizeStorage();
  _isAssembled = true;
  _toBeDestroyed = true;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace Trilinos

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
