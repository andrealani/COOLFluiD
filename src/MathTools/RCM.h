#ifndef COOLFluiD_MathTools_RCM_h
#define COOLFluiD_MathTools_RCM_h

#include <iostream>
#include <fstream>
#include <string>

#include "Common/ConnectivityTable.hh"
#include "MathTools/MathTools.hh"

namespace COOLFluiD {

class MathTools_API RCM
{
 public:
  /// Transforms a connectivity from Cell to Node to Node to Node
  static void transformCellNode2NodeNode (const Common::ConnectivityTable<CFuint>& cellstate, 
					  Common::ConnectivityTable<CFuint>& nodenode);
  
  /// Transforms a connectivity from Cell to Node to Node to Node
  static void transformCellNode2NodeNodeMedianDual (const Common::ConnectivityTable<CFuint>& cellstate, 
						    const Common::ConnectivityTable<CFuint>& cellnode,
						    Common::ConnectivityTable<CFuint>& nodenode);
  
  /// Applies the Reverse Cuthill-McKee algorithm to the graph of the connectivity passed
  /// @param new_id is a vector with the new id numbers
  static void renumber (Common::ConnectivityTable<CFuint>& cellstate, 
			Common::ConnectivityTable<CFuint>& cellnode, 
			std::valarray<CFuint>& new_id,
			const bool useMedianDual);
  
  /// reads the a cell to node connectivity from the file
  static int read_input (const std::string& filename, 
			 Common::ConnectivityTable<CFuint>& cellnode);
  
  /// prints the connectivity to a file
  static void print_table (const std::string& filename, 
			   const Common::ConnectivityTable<CFuint>& cellstate, 
			   const Common::ConnectivityTable<CFuint>& cellnode, 
			   const bool useMedianDual);
};

} // namespace COOLFluiD

#endif // COOLFluiD_MathTools_RCM_h
