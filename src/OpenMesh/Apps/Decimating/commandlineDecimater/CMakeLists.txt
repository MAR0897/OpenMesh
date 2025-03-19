include (VCICommon)

vci_add_executable (commandlineDecimater ../decimater.cc)

find_package(Eigen3)

target_link_libraries (commandlineDecimater
  OpenMeshCore
  OpenMeshTools
  Eigen3::Eigen
)

