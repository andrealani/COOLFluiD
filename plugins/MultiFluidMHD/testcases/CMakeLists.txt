cf_add_case( MPI	1	CASEDIR 1_HartmannFlow		PCASE hartmannRegressionTest.CFcase  		CASEFILES hartmann.inter halfTube40X250.neu ) 
cf_add_case( MPI	8	CASEDIR 2_MagneticReconnection 	PCASE reconnectionPsi05RegressionTest.CFcase 	CASEFILES chromosphericRec2D.inter Mesh_36x6_400x250.neu )
cf_add_case( MPI        8       CASEDIR 3_Chromostrat		PCASE UNSTEADYCASE_KhMagneticFieldRegressionTest.CFcase   CASEFILES Mesh_2x2_100x1500.neu MFinteractive_unsteady.inter )
cf_add_case( MPI        8       CASEDIR 4_ChromostratExtended   PCASE UNSTEADYCASE_KhMagneticFieldRegressionTest.CFcase   CASEFILES Mesh_2x2_300x300_extended.neu MFinteractive_unsteady.inter )
