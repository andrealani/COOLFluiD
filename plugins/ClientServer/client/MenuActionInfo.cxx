#include <QtCore>
#include <QtGui>

#include "ClientServer/client/MenuActionInfo.hh"

using namespace COOLFluiD::client;

MenuActionInfo::MenuActionInfo()
{
  this->initDefaults();
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void MenuActionInfo::initDefaults()
{
  m_text = "";
  m_enabled = true;
  m_checkable = false;
  m_icon = QIcon();
  m_shortcut = QKeySequence();
  
  m_slot = NULL;
  m_menu = NULL;
  m_toolbar = NULL;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

QAction * MenuActionInfo::buildAction(QObject * parent)
{
  QAction * action = NULL;
  
  if(parent != NULL)
  {
    action = new QAction(m_text, parent);
    action->setEnabled(m_enabled);
    
    if(!m_shortcut.isEmpty())
      action->setShortcut(m_shortcut);
    
    if(m_slot != NULL)
      connect(action, SIGNAL(triggered()), parent, m_slot);
    
    if(m_menu != NULL)
      m_menu->addAction(action);
    
    if(m_toolbar != NULL)
      m_toolbar->addAction(action);
    
    action->setIcon(m_icon);
    action->setCheckable(m_checkable);
    action->setIconVisibleInMenu(true);
  }
  
  return action;
}
