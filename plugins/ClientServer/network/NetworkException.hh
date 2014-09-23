#ifndef COOLFluiD_network_NetworkException_h
#define COOLFluiD_network_NetworkException_h

/////////////////////////////////////////////////////////////////////////////

#include <QString>

#include "Common/Exception.hh"
#include "Common/CodeLocation.hh"

#include "ClientServer/network/NetworkAPI.hh"

namespace COOLFluiD 
{
  namespace network 
  {
    
    namespace Common
    {
      class CodeLocation;
      class Exception;
    }
    
    /////////////////////////////////////////////////////////////////////////////
    
    /// @brief Exception thrown when the server can not open its m_socket.
    
    /// @author Quentin Gasper.
    
    class Network_API NetworkException : public COOLFluiD::Common::Exception
    {
      
    public:
      
      /// Constructor
      NetworkException(const COOLFluiD::Common::CodeLocation& where, 
                       const std::string& what);
      
      /// Copy constructor
      NetworkException(const NetworkException& e) throw ();
      
    };
    
    /////////////////////////////////////////////////////////////////////////////
    
  } //namespace Network 
} // namespace COOLFLuiD

/////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_network_NetworkException_h
