#include "AutoTemplateLoader/NschemeTVMET.hh"
#include "AutoTemplateLoader/AutoTemplateLoader.hh"
#include "AutoTemplateLoader/RunNscheme.hh"
#include "Environment/ObjectProvider.hh"

COOLFluiD::Environment::ObjectProvider< RunNscheme< COOLFluiD::Benchmark::NschemeTVMET<4,3> > ,
                                        BaseRunLib,
                                        COOLFluiD::AutoTemplateLoader::AutoTemplateLoaderModule >
aRunNschemeTVMET_4_3_provider("NschemeTVMET_4_3");

COOLFluiD::Environment::ObjectProvider< RunNscheme< COOLFluiD::Benchmark::NschemeTVMET<3,3> > ,
                                        BaseRunLib,
                                        COOLFluiD::AutoTemplateLoader::AutoTemplateLoaderModule >
aRunNschemeTVMET_3_3_provider("NschemeTVMET_3_3");


