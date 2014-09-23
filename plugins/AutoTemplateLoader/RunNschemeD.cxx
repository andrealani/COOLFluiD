#include "AutoTemplateLoader/NschemeD.hh"
#include "AutoTemplateLoader/AutoTemplateLoader.hh"
#include "AutoTemplateLoader/RunNscheme.hh"
#include "Environment/ObjectProvider.hh"
       
COOLFluiD::Environment::ObjectProvider< RunNscheme< COOLFluiD::Benchmark::NschemeD<5,4> > ,
                                        BaseRunLib,
                                        COOLFluiD::AutoTemplateLoader::AutoTemplateLoaderModule >
aRunNschemeD_5_4_provider("NschemeD_5_4");

COOLFluiD::Environment::ObjectProvider< RunNscheme< COOLFluiD::Benchmark::NschemeD<4,3> > ,
                                        BaseRunLib,
                                        COOLFluiD::AutoTemplateLoader::AutoTemplateLoaderModule >
aRunNschemeD_4_3_provider("NschemeD_4_3");

COOLFluiD::Environment::ObjectProvider< RunNscheme< COOLFluiD::Benchmark::NschemeD<3,3> > ,
                                        BaseRunLib,
                                        COOLFluiD::AutoTemplateLoader::AutoTemplateLoaderModule >
aRunNschemeD_3_3_provider("NschemeD_3_3");



