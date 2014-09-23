#ifndef COOLFluiD_network_IOException_h
#define COOLFluiD_network_IOException_h

/////////////////////////////////////////////////////////////////////////////

#include <QString>

#include "Common/Exception.hh"

#include "ClientServer/network/NetworkAPI.hh"

namespace COOLFluiD 
{
  namespace network 
  {
    
    /////////////////////////////////////////////////////////////////////////////
    
    /// @brief Exception thrown when an In/Out error occurs.
    
    /// @author Quentin Gasper.
    
    class Network_API IOException : public COOLFluiD::Common::Exception
    {
      
    public:
      
      /// Constructor
      IOException(const COOLFluiD::Common::CodeLocation& where, 
                  const std::string& what);
      
      /// Copy constructor
      IOException(const IOException& e) throw ();
      
    };
    
    /////////////////////////////////////////////////////////////////////////////
    
  } //namespace Network 
} // namespace COOLFLuiD

/////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_network_IOException_h
