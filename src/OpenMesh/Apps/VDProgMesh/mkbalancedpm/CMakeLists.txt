include (VCICommon)

vci_add_executable (mkbalancedpm mkbalancedpm.cc)

find_package(Eigen3)

target_link_libraries (mkbalancedpm
  OpenMeshCore
  OpenMeshTools
  Eigen3::Eigen
)

