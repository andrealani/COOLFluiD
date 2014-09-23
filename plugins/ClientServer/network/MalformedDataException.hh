#ifndef COOLFluiD_network_MalformedDataException_h
#define COOLFluiD_network_MalformedDataException_h

/////////////////////////////////////////////////////////////////////////////

#include <QString>

#include "Common/Exception.hh"

#include "ClientServer/network/NetworkAPI.hh"

namespace COOLFluiD 
{
  namespace network 
  {
    
    /////////////////////////////////////////////////////////////////////////////
    
    /// @brief Exception thrown when data do not have excpected format.
    
    /// @author Quentin Gasper.
    
    class Network_API MalformedDataException : public COOLFluiD::Common::Exception
    {
      
    public:
      
      /// Constructor
      MalformedDataException(const COOLFluiD::Common::CodeLocation& where, 
                             const std::string& what);
      
      /// Copy constructor
      MalformedDataException(const MalformedDataException& e) throw ();
      
    };
    
    /////////////////////////////////////////////////////////////////////////////
    
  } //namespace Network 
} // namespace COOLFLuiD

/////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_network_MalformedDataException_h
