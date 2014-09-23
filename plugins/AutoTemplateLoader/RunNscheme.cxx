#include "AutoTemplateLoader/Nscheme.hh"
#include "AutoTemplateLoader/AutoTemplateLoader.hh"
#include "AutoTemplateLoader/RunNscheme.hh"
#include "Environment/ObjectProvider.hh"

COOLFluiD::Environment::ObjectProvider< RunNscheme< COOLFluiD::Benchmark::Nscheme<5,4> > ,
                                        BaseRunLib,
                                        COOLFluiD::AutoTemplateLoader::AutoTemplateLoaderModule >
aRunNscheme_5_4_provider("Nscheme_5_4");

COOLFluiD::Environment::ObjectProvider< RunNscheme< COOLFluiD::Benchmark::Nscheme<4,3> > ,
                                        BaseRunLib,
                                        COOLFluiD::AutoTemplateLoader::AutoTemplateLoaderModule >
aRunNscheme_4_3_provider("Nscheme_4_3");

COOLFluiD::Environment::ObjectProvider< RunNscheme< COOLFluiD::Benchmark::Nscheme<3,3> > ,
                                        BaseRunLib,
                                        COOLFluiD::AutoTemplateLoader::AutoTemplateLoaderModule >
aRunNscheme_3_3_provider("Nscheme_3_3");



