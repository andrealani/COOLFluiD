// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Common/PE.hh"
#include "Framework/BlockAccumulator.hh"
#include "Framework/DataHandle.hh"
#include "Paralution/ParalutionMatrix.hh"

#include "Framework/DataSocketSink.hh"

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace Paralution {

//////////////////////////////////////////////////////////////////////////////

ParalutionMatrix::ParalutionMatrix() :
  Framework::LSSMatrix(),
  socket_states("states"),
  socket_nodes("nodes")
{

 // Framework::DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
 // Framework::DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
 // bool useNodeBased = getMethodData().useNodeBased();
 // CFuint nbStates = (!useNodeBased) ? states.size() : nodes.size();
 firstIter = true;
 CFLog(VERBOSE,"ParalutionMatrix created \n");

}
      
//////////////////////////////////////////////////////////////////////////////

ParalutionMatrix::~ParalutionMatrix()
{
}

//////////////////////////////////////////////////////////////////////////////

void ParalutionMatrix::createSeqAIJ(const CFint m,
                               const CFint n,
                               const CFint nz,
                               const CFint* nnz,
                               const char* name)
{
  // // creation of the matrix
  // CF_CHKERRCONTINUE(MatCreateSeqAIJ(PETSC_COMM_SELF, m, n, nz, nnz, &m_mat));

  // if (name != CFNULL) {
  //   CF_CHKERRCONTINUE(ParalutionObjectSetName((ParalutionObject) m_mat, name));
  // }
}

//////////////////////////////////////////////////////////////////////////////

void ParalutionMatrix::createSeqBAIJ(const CFuint blockSize,
                                const CFint m,
                                const CFint n,
                                const CFint nz,
                                const CFint* nnz,
                                const char* name)
{

}

//////////////////////////////////////////////////////////////////////////////

#ifdef CF_HAVE_MPI
void ParalutionMatrix::createParAIJ(MPI_Comm comm,
				    const CFint m,
				    const CFint n,
				    const CFint M,
				    const CFint N,
				    const CFint dnz,
				    const CFint* dnnz,
				    const CFint onz,
				    const CFint* onnz,
				    const char* name)
{
}
#endif

//////////////////////////////////////////////////////////////////////////////

#ifdef CF_HAVE_MPI
void ParalutionMatrix::createParBAIJ(MPI_Comm comm,
				     const CFuint blockSize,
				     const CFint m,
				     const CFint n,
				     const CFint M,
				     const CFint N,
				     const CFint dnz,
				     const CFint* dnnz,
				     const CFint onz,
				     const CFint* onnz,
				     const char* name)
{
}
#endif

//////////////////////////////////////////////////////////////////////////////

void ParalutionMatrix::createAIJ(const CFuint BlockSize,
                                 const CFuint nbRows,
                                 const CFuint globalSize)
{
   CFuint ValPerBlock = BlockSize*BlockSize;
   CFuint ValPerRow = ValPerBlock*5;   //4 neighbour plus the cell itself!! REALLY HARDCODED AND AWFUL!!
   CFuint TotalSize = ValPerRow*globalSize; //Maximun number of nonzerovalues of the matrix

  // _val.resize(TotalSize);   
  // _col.resize(TotalSize);   
  // _row.resize(TotalSize);   
  abort();
  //_size = TotalSize;
  //std::cout << "Total size (nnz): " << TotalSize << "\n";
  //_nnz = 0;
}


//////////////////////////////////////////////////////////////////////////////


void ParalutionMatrix::createCSR(std::valarray<CFint> allNonZero, CFuint nbEqs)
{
   
   CFuint nbStates = allNonZero.size();
   _rowlength = nbStates+1;

   _rowoff = new CFint[_rowlength];  //.resize(_rowlength);
   _size = 0;
   _rowoff[0] = 0;
/*
  The allNonZero can store {0,3,5,...} indicating the number of neigbours,
  or {0, 54, 126} indicating the exact index of the offset. In the first one it is needed to multiply by nbEqs

   if (allNonZero[1]<=nbEqs){
     for (CFuint i=0; i<nbStates; i++){ 
        _rowoff[i+1] = _size+allNonZero[i]*nbEqs;
        _size += allNonZero[i]*nbEqs;
       // std::cout << _rowoff[i] << "\t";
     }
   }else{
*/
     for (CFuint i=0; i<nbStates; i++){ 
        _rowoff[i+1] = _size+allNonZero[i];
        _size += allNonZero[i];
      //  std::cout << _rowoff[i] << "\t";
     }
//   }
 //  _row.resize(_size);
   _val = new CFreal[_size]; //.resize(_size);   
   _col = new CFint[_size];   

   CFLog(NOTICE, "Allocated: " << _size << "\t nbStates: " << _rowoff[nbStates] << "\n");
 // _nnz = 0;
}


//////////////////////////////////////////////////////////////////////////////

void ParalutionMatrix::setValues(const Framework::BlockAccumulator& acc)
{
  // CFLog(DEBUG_MIN, "ParalutionMatrix::setValues()\n");
  // CF_CHKERRCONTINUE( MatSetValuesBlocked(m_mat,acc.getM(),&acc.getIM()[0],acc.getN(),&acc.getIN()[0],
  //   const_cast<Framework::BlockAccumulator&>(acc).getPtr(), INSERT_VALUES) );
}
      
//////////////////////////////////////////////////////////////////////////////

void ParalutionMatrix::addValues(const Framework::BlockAccumulator& acc)
{
  using namespace std;

  CFuint m = acc.getM();    //Number of neighbors
  CFuint n = acc.getN();    //Central cell
  CFuint nb = acc.getNB();  //nb eqs
  const std::vector<int>& nbID = acc.getIM(); //Array storing indexes of neigbours
  const std::vector<int>& IDs = acc.getIN(); //Array storing neigbours
  CFuint IndexCSR;
  CFreal val;
  CFuint RowPosition;
  CFuint RowPositionPlusOne;
  CFuint mm;
  CFuint index = 0; //Also number of nonzerovalues
  CFreal* values = const_cast<Framework::BlockAccumulator&>(acc).getPtr();
      CFuint counter = 0;                //Number of non-connected neigbours
      for (CFint mi=0; mi<m; mi++){      //Loop over the neigbours
         if (nbID[mi] == -1) {
           counter++;
         }else{
           for (CFint nbj=0;nbj<nb; nbj++){
              RowPosition = _rowoff[nbID[mi]*nb + nbj];  //Neighbour row position
              RowPositionPlusOne = _rowoff[nbID[mi]*nb + nbj + 1]; //index of next row
              mm = (RowPositionPlusOne-RowPosition)/nb;  //Number of nonzero values of the neighbour
              IndexCSR = -1;
              for (CFint mii=0; mii<mm; mii++){  //Look for the correct index looping over the neighbors of the neigbour
                 if (_col[RowPosition+mii*nb] == IDs[0]*nb || _col[RowPosition+mii*nb] == -1){
                    IndexCSR = RowPosition+mii*nb;
                 }
                
              }

              for (CFint nbi=0; nbi<nb; nbi++){
		
                val = values[mi*nb*nb + nbj*nb + nbi]; 
                _col[IndexCSR+nbi] = IDs[0]*nb + nbi;
                _val[IndexCSR+nbi] += val;
              }
           }
         }
      }
}
      


//////////////////////////////////////////////////////////////////////////////

void ParalutionMatrix::printToScreen() const
{
  // CF_CHKERRCONTINUE(MatAssemblyBegin(m_mat,MAT_FINAL_ASSEMBLY));
  // CF_CHKERRCONTINUE(MatAssemblyEnd(m_mat,MAT_FINAL_ASSEMBLY));
  // CF_CHKERRCONTINUE(MatView(m_mat, PETSC_VIEWER_STDOUT_WORLD));
}

//////////////////////////////////////////////////////////////////////////////

void ParalutionMatrix::printToFile(const char* fileName) const
{

  m_mat.WriteFileMTX(fileName);
  // ParalutionViewer viewer;
  // CF_CHKERRCONTINUE(MatAssemblyBegin(m_mat,MAT_FINAL_ASSEMBLY));
  // CF_CHKERRCONTINUE(MatAssemblyEnd(m_mat,MAT_FINAL_ASSEMBLY));
  
  // const std::string nsp = Framework::MeshDataStack::getActive()->getPrimaryNamespace();
  // CF_CHKERRCONTINUE(ParalutionViewerASCIIOpen(Common::PE::GetPE().GetCommunicator(nsp),fileName,&viewer));
  // CF_CHKERRCONTINUE(ParalutionViewerSetType(viewer,PETSCVIEWERASCII));
  // CF_CHKERRCONTINUE(ParalutionViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATRIXMARKET));
  // // PETSC_VIEWER_ASCII_COMMON, PETSC_VIEWER_ASCII_INDEX, PETSC_VIEWER_ASCII_MATLAB
  // CF_CHKERRCONTINUE(MatView(m_mat, viewer));
  // CF_CHKERRCONTINUE(ParalutionViewerDestroy(&viewer));
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace Paralution

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
