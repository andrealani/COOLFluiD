#ifndef COOLFluiD_server_TypesNotFoundException_h
#define COOLFluiD_server_TypesNotFoundException_h

/////////////////////////////////////////////////////////////////////////////

#include "Common/Exception.hh"

namespace COOLFluiD
{
  
  namespace server
  {
    
    /////////////////////////////////////////////////////////////////////////////
    
    /// @brief Exception thrown when types could not be loaded.
    
    /// @author Quentin Gasper.
    
    class TypesNotFoundException : public COOLFluiD::Common::Exception
    {
    public:
      
      /// Constructor
      TypesNotFoundException(const COOLFluiD::Common::CodeLocation& where, 
                             const std::string& what);
      
      /// Copy constructor
      TypesNotFoundException(const TypesNotFoundException& e) throw ();
      
    }; // class TypesNotFoundException
    
    /////////////////////////////////////////////////////////////////////////////
    
  } // namespace server
} // namespace COOLFluiD

/////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_server_TypesNotFoundException_h
