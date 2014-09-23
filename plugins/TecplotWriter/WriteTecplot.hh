// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_TecplotWriter_WriteTecplot_hh
#define COOLFluiD_TecplotWriter_WriteTecplot_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataProcessingData.hh"
#include "Framework/DataSocketSink.hh"
#include "Common/NonCopyable.hh"
#include "MathTools/RealVector.hh"
#include "Framework/Framework.hh"
#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/FileHandlerOutput.hh"
#include <iostream>

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace TecplotWriter {

//////////////////////////////////////////////////////////////////////////////

 /*   This class computes the Vector Potential
  *   usage:
  *   WriteTecplot::getInstance().[...]
  *
  *   @author Emanuele Sartori
  */
class WriteTecplot : public Common::NonCopyable<WriteTecplot> {

//////////////////////////////////////////////////////////////////////////////
//                      Class Write Tecplot File                            //
//                                                                          //
//  USAGE:                                                                  //
//  - execute .setup() method                                               //
//  - set volumes and nodes using .setDataSockets(.,.) method               //
//  - set title using .setTitle("title") method                             //
//  - set extrapolation strategy (from cell center values to node center    //
//    velues): 1=1/4, 2=volume mean.                                        //
//                                                                          //
//  To add a variable to the data set:                                      //
//  - use .addOutput("variableName",vector[0],scaling,alreasyInNode) to     //
//    add variables to the tecplot file: vector[0] is the first element     //
//    of the vector you wanto to plot, scaling is a scaling factor you      //
//    may need to use, alreadyInNode is a boolean value (true if the        //
//    vector contains data in nodes and false if it contains values in      //
//    cell centers.                                                         //
//                                                                          //
//  At this point you can:                                                  //
//  - run .writeFileStructure("FileName") to prepare the tecplot file       //
//    and then .writeOutput("FileName") to append data to the file;         //
//  - run .writeTecplot("FileName") to write structure and data at the      //
//    same time.                                                            //
//                                                                          //
//  Remember that if you are passing a local vector as a parameter of       //
//  method .addOutput, you have to execute .writeOutput before the          //
//  local vector is deallocated.                                            //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

public:
  /// Default destructor
  ~WriteTecplot()
  {
  }

  /// Returns the instance of this VectorPotential
  /// This is the access point to the Singleton
  /// @return the instance of the singleton
  static WriteTecplot& getInstance()
  {
    static WriteTecplot aWriteTecplot;
    return aWriteTecplot;
  } 


  /*
   *  Add an output variable
   */
  void addOutput(std::string Name, CFreal* aVectorPointer, CFreal Scale, bool inNodes = false)
  {
    //std::cout << " adding variable...\n" ;

    l_variableNames.push_back(Name);
    l_variableVectorPointers.push_back(aVectorPointer);
    l_variableScale.push_back(Scale);
    l_inNodes.push_back(inNodes);
  }


  /*
   *  Reset
   */
  void resetOutput()
  {
    //l_variableNames.clear;
    //l_variableVectorPointers.clear;
    l_variableNames.resize(0);
    l_variableVectorPointers.resize(0);
    l_variableScale.resize(0);
    l_inNodes.resize(0);
    l_title = "";
  }


  /*
   *  Set filename
   */
  void setTitle(std::string Title)
  {
    //std::cout << " setting title...\n" ;
    l_title = Title;
  }


  /*
   *  Set strategy to compute node values
   */
  void setNodeExtrapolation(CFuint strategy)
  {
    //std::cout << " setting node extrapolation...\n" ;

    // volumes or 1/4 or distance?
    if (strategy==0) //default
          strategy=1;
    l_strategy = strategy;
  }


  /*
   *  Sets the DataSockets
   */
  void setDataSockets
  (COOLFluiD::Framework::DataSocketSink<CFreal> volumesSocket,
  COOLFluiD::Framework::DataSocketSink<COOLFluiD::Framework::Node*, COOLFluiD::Framework::GLOBAL> nodesSocket);


  /*
   *  setup
   */
  void setup();


  /*
   *  Call this method when you want to write a tecplot file.
   *
   *  It is equivalent to call first writeFileStructure() and
   *  then writeOutput().
   *  Please note that if you use this function, the output queue will be cleared.
   */
  void writeFile(std::string FileName);


  /*
   *  Call this method when you want to prepare a tecplot file.
   *
   *  It writes the tecplot file title, structure, cells connectivity...
   */
  void writeFileStructure(std::string FileName);


  /*
   *  Call this method when you want to add new variables to a tecplot file.
   *
   *  It writes variables output when the file structure is already there.
   *  Please note that if you use this function, the output queue will be cleared.
   *  
   */
  void writeOutput(std::string FileName);

  void ReIm_TO_ModPhase(const CFreal Re, const CFreal Im, CFreal &Modulo, CFreal &Phase);

  CFreal normalizeAngle(CFreal DegreeAngle);

private:
  /**
   * Constructor
   */
  WriteTecplot();

  /**
   * from cell center values to node values
   */
  void prepareNodeValues(CFuint strategy, std::vector<std::vector<CFreal> > &outputVector);


private: //data

  /// socket for volumes
  COOLFluiD::Framework::DataSocketSink<CFreal> socket_volumes;

  /// socket for nodes
  COOLFluiD::Framework::DataSocketSink < COOLFluiD::Framework::Node* , 
         COOLFluiD::Framework::GLOBAL > socket_nodes;

  // variable names
  std::vector<std::string> l_variableNames;

  // vector pointers
  std::vector<CFreal*> l_variableVectorPointers;

  // vector: scale factors
  std::vector<CFreal> l_variableScale;

  // vector: values in cells or in nodes?!?
  std::vector<bool> l_inNodes;

  // output file name
  std::string l_title;

  // strategy to be used to compute values in nodes
  CFuint l_strategy;

  /// builder of geometric entities
  COOLFluiD::Framework::GeometricEntityPool<COOLFluiD::Framework::StdTrsGeoBuilder> m_geoBuilder;

  ///flag to known if the socket have already been set
  bool m_isSocketsSet;

}; // end class WriteTecplot

//////////////////////////////////////////////////////////////////////////////

    } // namespace TecplotWriter
  
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_TecplotWriter_WriteTecplot_hh
