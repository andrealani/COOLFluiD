#ifndef COOLFluiD_server_UnknownClientIdException_h
#define COOLFluiD_server_UnknownClientIdException_h

/////////////////////////////////////////////////////////////////////////////

#include "Common/Exception.hh"

namespace COOLFluiD
{
  
  namespace server
  {
    
    /////////////////////////////////////////////////////////////////////////////
    /// @brief Exception thrown when a given client id could not be associated to
    /// any existing client.
    
    /// @author Quentin Gasper.
    class UnknownClientIdException : public COOLFluiD::Common::Exception
    {
    public:
      
      /// Constructor
      UnknownClientIdException(const COOLFluiD::Common::CodeLocation& where, 
                               const std::string& what);
      
      /// Copy constructor
      UnknownClientIdException(const UnknownClientIdException& e) throw ();
      
    }; // class UnknownClientIdException
    
    /////////////////////////////////////////////////////////////////////////////
    
  } // namespace server
} // namespace COOLFluiD

/////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_server_UnknownClientIdException_h
