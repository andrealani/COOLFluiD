#include <iostream>

#include <QtCore>
#include <QApplication>

#include "ClientServer/client/MainWindow.hh"

using namespace std;

using namespace COOLFluiD::client;

int main(int argc, char *argv[])
{
  QApplication app(argc, argv);
  
  MainWindow window;
  window.showMaximized();
  
  return app.exec();
}
