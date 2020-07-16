#===========================================================
ADD_SUBDIRECTORY("${CMAKE_SOURCE_DIR}/Python/esys/lsm/util")
ADD_SUBDIRECTORY("${CMAKE_SOURCE_DIR}/Python/esys/lsm/geometry")

ADD_LIBRARY(libLsmPy SHARED LsmMpiPy.cpp RunnablePy.cpp CheckPointParamsPy.cpp InteractionParamsPy.cpp InteractionGroupPy.cpp BondInteractionGroupPy.cpp BondedTriMeshPrmsPy.cpp BondedMesh2DPrmsPy.cpp MeshBuildParamsPy.cpp WallPrmsPy.cpp ParticlePy.cpp FieldSaverPrmsPy.cpp ParticleFieldSaverPrmsPy.cpp InteractionFieldSaverPrmsPy.cpp ElasticTriMeshPrmsPy.cpp WallFieldSaverPrmsPy.cpp TriangleFieldSaverPrmsPy.cpp RotParticlePy.cpp LmParticleAdder.cpp TriggerPrmsPy.cpp ElasticMesh2DPrmsPy.cpp ParticleIdPairSetPy.cpp ParticleIdPairPy.cpp ParticleIdPairVectorPy.cpp RotParticleViPy.cpp RotThermalParticlePy.cpp)
ADD_LIBRARY(LsmPy SHARED LsmMpiPy.cpp RunnablePy.cpp CheckPointParamsPy.cpp InteractionParamsPy.cpp InteractionGroupPy.cpp BondInteractionGroupPy.cpp BondedTriMeshPrmsPy.cpp BondedMesh2DPrmsPy.cpp MeshBuildParamsPy.cpp WallPrmsPy.cpp ParticlePy.cpp FieldSaverPrmsPy.cpp ParticleFieldSaverPrmsPy.cpp InteractionFieldSaverPrmsPy.cpp ElasticTriMeshPrmsPy.cpp WallFieldSaverPrmsPy.cpp TriangleFieldSaverPrmsPy.cpp RotParticlePy.cpp LmParticleAdder.cpp TriggerPrmsPy.cpp ElasticMesh2DPrmsPy.cpp ParticleIdPairSetPy.cpp ParticleIdPairPy.cpp ParticleIdPairVectorPy.cpp RotParticleViPy.cpp RotThermalParticlePy.cpp)

TARGET_LINK_LIBRARIES(libLsmPy Parallel Foundation Geometry Model libFoundationPy BoostPythonUtil ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})
TARGET_LINK_LIBRARIES(LsmPy Parallel Foundation Geometry Model libFoundationPy BoostPythonUtil ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})
INSTALL(TARGETS libLsmPy DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
INSTALL(TARGETS LsmPy DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
