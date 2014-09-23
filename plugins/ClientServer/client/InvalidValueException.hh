#ifndef COOLFluiD_client_InvalidValueException_h
#define COOLFluiD_client_InvalidValueException_h

/////////////////////////////////////////////////////////////////////////////

#include "Common/Exception.hh"

namespace COOLFluiD
{
  
  namespace client
  {
    
    /////////////////////////////////////////////////////////////////////////////
    /// @brief Exception thrown by @c GraphicalOption class when a QVariant 
    /// object can not be convert to the desired type of data.
    
    /// @author Quentin Gasper.
    class InvalidValueException : public COOLFluiD::Common::Exception
    {
    public:
      
      /// Constructor
      InvalidValueException(const COOLFluiD::Common::CodeLocation& where, 
                            const std::string& what);
      
      /// Copy constructor
      InvalidValueException(const InvalidValueException& e) throw ();
      
    }; // class InvalidValueException
    
    /////////////////////////////////////////////////////////////////////////////
    
  } // namespace client
} // namespace COOLFluiD

/////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_client_InvalidValueException_h
