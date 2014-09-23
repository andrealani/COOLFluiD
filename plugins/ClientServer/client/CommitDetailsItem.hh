#ifndef COOLFluiD_client_CommiDetailsItem_h
#define COOLFluiD_client_CommiDetailsItem_h

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD 
{
  namespace client 
  {
    
    /////////////////////////////////////////////////////////////////////////////
    
    class CommitDetailsItem 
    {
      
    public:
      
      CommitDetailsItem(const QString & optionName, const QString & oldValue, 
                        const QString & currentValue);
      
      CommitDetailsItem(const QString & oldValue, const QString & currentValue);
      
      bool isNewOption() const;
      
      QString getCurrentValue() const;
      
      QString getOldValue() const;
      
      QString getOptionName() const;
      
    private:
      
      QString m_optionName;
      
      QString m_oldValue;
      
      QString m_currentValue;
      
      bool m_newOption;
      
      
    }; // class CommitDetailsItem
    
    ////////////////////////////////////////////////////////////////////////////
    
  } // namespace client
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif