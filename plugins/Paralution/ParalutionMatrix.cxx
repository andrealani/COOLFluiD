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

   _val = new CFreal[TotalSize];   
   _col = new CFint[TotalSize];   
   _row = new CFint[TotalSize];   

  _size = TotalSize;
  std::cout << "Total size (nnz): " << TotalSize << "\n";
  _nnz = 0;
}


//////////////////////////////////////////////////////////////////////////////


void ParalutionMatrix::createCSR(std::valarray<CFint> allNonZero, CFuint nbEqs)
{
   
   CFuint nbStates = allNonZero.size();

   _rowoff = new CFint[nbStates+1];
   _size = 0;
   _rowoff[0] = 0;
   for (CFuint i=0; i<nbStates; i++){ 
      _rowoff[i+1] = _size+allNonZero[i];
      _size += allNonZero[i];
      //std::cout << _rowoff[i] << "\t";
   }
  
   _row = new CFint[_size];
   _val = new CFreal[_size];   
   _col = new CFint[_size];   

   std::cout << "Total _size: " << _size << "\n";
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
 // CFLog(NOTICE, "ParalutionMatrix::addValues() \n");
  using namespace std;
  RealVector rows;
  RealVector cols;
  RealVector values;
  CFuint size = acc.size(); //number of entries
  rows.resize(size);
  cols.resize(size);
  values.resize(size);  

  CFuint n = acc.getM();    //Number of neighbors
  CFuint m = acc.getN();    //Central cell
  CFuint nb = acc.getNB();  //nb eqs

  const std::vector<int> nbID = acc.getIM(); //Array storing indexes of neigbours
  const std::vector<int> IDs = acc.getIN(); //Array storing neigbours
//cout << IDs.size() << nbID.size() << endl;
  CFuint IndexCOO;
  CFuint index = 0; //Also number of nonzerovalues


//  const std::vector<CFint>& workspace = acc.getWorkspace();
  
//  cout << workspace.size() << "\n";

  for (CFuint mi = 0; mi<m; mi++){
          CFint counter=0;
    for (CFuint ni = 0; ni<n; ni++){

          if(nbID[ni] == -1){
             counter++;
          }
       for (CFuint nbi=0; nbi<nb; nbi++){
         for(CFuint nbj=0; nbj<nb; nbj++){
          if (nbID[ni] != -1){
           CFreal val = acc.getValue(ni, mi, nbi, nbj);

           IndexCOO = _rowoff[IDs[mi]*nb+nbi] + (ni-counter)*nb + nbj;  
           
           _row[IndexCOO] = IDs[mi]*nb+nbi;
           _col[IndexCOO] = nbID[ni]*nb+nbj;
       //if (_val[IndexCOO] == 0){
           _val[IndexCOO] = val;
//}

          //     cout << IndexCOO << " " << IDs[mi]*nb+nbi << " " << nbID[ni]*nb+nbj <<  " " <<  val << "\n";
   //        if (IDs[mi] == 0) {
   //             cout << mi << " " << IDs[mi] << " " << ni << " " << counter << " " << nbID[ni] << " " << nbi << " " << nbj <<  " " <<  IndexCOO << " " << _col[IndexCOO]<<"\n";
    //       }
           
  //         cout << "(" << _row[IndexCOO] << "," << _col[IndexCOO] << ",): " << _val[IndexCOO] << " \t" << IndexCOO << "\n";
         //  if (val != 0.0){
         //    rows[index] = IDs[mi]*nb+nbi;       //Maybe nbi and nbj the other way around?
         //    cols[index] = nbID[ni]*nb+nbj; 
         //    values[index] = val; 
         //    index++; 
         //  }
          }
         }
       }
    }
  }
//abort();
//  cout << "Non-zero values: " << index << ", total size: " << size << " \n";
//  for (CFuint i=0; i<index; i++){
//     cout << "(" << rows[i] << "," << cols[i] << ",): " << values[i] << " \t";
//  }
  //addtoCOO
  
//[artesla2:00438] Signal: Segmentation fault (11)
//[artesla2:00438] Signal code: Address not mapped (1)
//[artesla2:00438] Failing at address: (nil)

//    std::cout << _nnz << "\n"; //to see the out of range??? 0!!
/*
    for (CFuint i=0; i<index; i++){
    bool found = false;
      for (CFuint _i=0; _i<_nnz; _i++){
         if (rows[i]==_row[_i] && cols[i]==_col[_i]){
            _val[_i] += values[i]; found=true; 
         }
      }
      if (found==false){
        _row[_nnz] = rows[i];
        _col[_nnz] = cols[i];
        _val[_nnz] = values[i];
        _nnz++;
      }
   }
*/
  // CFLog(DEBUG_MIN, "ParalutionMatrix::addValues()\n");
  // CF_CHKERRCONTINUE( MatSetValuesBlocked(m_mat,acc.getM(),&acc.getIM()[0],acc.getN(),&acc.getIN()[0],
  // 					 const_cast<Framework::BlockAccumulator&>(acc).getPtr(), ADD_VALUES) );
//  cout << _nnz << endl;
//  for (CFuint i=0; i<_nnz; i++){
//     cout << "(" << _row[i] << "," << _col[i] << "): " << _val[i] << " \t";
//  }

}
      
//////////////////////////////////////////////////////////////////////////////

void ParalutionMatrix::addtoCOO(CFuint* col, CFuint* row, CFreal* val, CFuint nnz)
{

   for (CFuint i=0; i<nnz; i++){
      for (CFuint _i=0; _i<_nnz; _i++){
         if (row[i]==_row[_i] && col[i]==_col[_i]){
            _val[_i] += val[i];
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
