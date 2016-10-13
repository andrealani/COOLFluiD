cf_add_case( MPI default PCASE WireSteady/wireSteadyImpl.CFcase )
cf_add_case( MPI 1 2     PCASE wireSteadyImplProjBigCircle100X30.CFcase )
cf_add_case( MPI	1	CASEDIR WireSteady		PCASE wireSteadyImplRegressionCUDA.CFcase  		CASEFILES circle80x110.neu ) 
