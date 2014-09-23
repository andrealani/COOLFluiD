#include <cassert>

#include "AutoTemplateLoader/AutoTemplateLoader.hh"
#include "AutoTemplateLoader/RunRoeFluxFlens.hh"
#include "Environment/ObjectProvider.hh"
       
COOLFluiD::Environment::ObjectProvider< RunRoeFluxFlens<double,4> ,
                                        BaseRunLib,
                                        COOLFluiD::AutoTemplateLoader::AutoTemplateLoaderModule >
aRunRoeFluxFlens_provider("RunRoeFluxFlens");
