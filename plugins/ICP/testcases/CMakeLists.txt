cf_add_case( MPI 8  CASEDIR ICP_EQProbe_Step1_Torch_NoElectricField PCASE ICP_EQProbe_Step1_Torch_NoElectricField.CFcase CASEFILES spherical_torch.CFmesh spherical_torch.neu torch.inter )
cf_add_case( MPI 8  CASEDIR ICP_EQProbe_Step2_Torch_ElectricField   PCASE ICP_EQProbe_Step2_Torch_ElectricField.CFcase   CASEFILES torch_1st_out.CFmesh )
cf_add_case( MPI 8  CASEDIR ICP_EQProbe_Step4_FullDomain_1stOrder   PCASE ICP_EQProbe_Step4_FullDomain_1stOrder.CFcase   CASEFILES afterextrap.CFmesh )
cf_add_case( MPI 8  CASEDIR ICP_EQProbe_Step5_FullDomain_2ndOrder   PCASE ICP_EQProbe_Step5_FullDomain_2ndOrder.CFcase   CASEFILES restart.CFmesh )
