#ifndef COOLFluiD_client_HostListPanel_h
#define COOLFluiD_client_HostListPanel_h

//////////////////////////////////////////////////////////////////////////////

#include <QListView>
#include <QHash>

class QString;
class QStandardItemModel;
class QStandardItem;
class QDomDocument;

namespace COOLFluiD
{
  namespace client
  {
    
    //////////////////////////////////////////////////////////////////////////////
    
    /// @brief List m_view that allows user to select hosts.
    
    /// The data for the internal model are provided through a XML tree, with the
    /// following format:
    /// @code
    /// <hosts>
    ///  <item name="hostname" selected="false">name to display</item>
    ///  <!-- (...) -->
    /// </host>
    /// @endcode
    /// Where: @n
    /// @li @c name gives the hostname
    /// @li @c selected gives the the checkstate of the associated checkbox. 
    /// @c true for a checked state, @c false for an unchecked state.
    /// @li tag data are the string to display next to the checkbox
    
    /// @author Quentin Gasper.
    
    class HostListPanel : public QListView
    {
    public:
      
      /// @brief Constructor
      HostListPanel();
      
      /// @brief Destructor.
      
      /// Frees all allocated memory.
      ~HostListPanel();
      
      /// @brief Sets model m_items.
      
      /// Model is previously cleared.
      void setItems(const QDomDocument & doc);
      
      /// @brief Gives a modified deep copy of internal of the data tree. 
      
      /// This tree represents modifcations that were done since the last call
      /// of @c #setItems.
      /// @return Returns a modified deep copy of internal of the data tree. 
      QDomDocument getDocument() const;
      
      /// @brief Indicates whether checkboxe states where modified since the
      /// last call to @c #setItems.
      
      /// @return Returns @c true if at least one checkbox state has changed;
      /// otherwise returns @c false.
      bool isModified() const;
      
    private:
      
      /// @brief Model.
      QStandardItemModel * m_viewModel;
      
      /// @brief Model m_items.
      QList<QStandardItem *> m_viewItems;
      
      /// @brief Host list.
      QDomDocument m_hosts;
      
      /// @brief Clears model and deletes m_items.
      void clear();
      
      
      
    }; // class HostListPanel
    
    //////////////////////////////////////////////////////////////////////////////
    
  } // namespace client
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_client_HostListPanel_h
