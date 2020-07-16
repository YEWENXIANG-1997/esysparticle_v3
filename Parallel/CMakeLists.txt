#===========================================================
ADD_LIBRARY(Parallel SHARED ASubLattice.cpp LatticeMaster.cpp mpia2abuf.cpp mpibarrier.cpp mpibuf.cpp mpicmdbuf.cpp mpisgbuf.cpp mpisgvbuf.cpp mpivbuf.cpp SubLatticeControler.cpp SubLattice.cpp CheckPointer.cpp CheckPointController.cpp CheckPointParams.cpp MpiInfo.cpp GeometryReader.cpp CheckPointInfo.cpp CheckPointLoader.cpp LatticeParam.cpp BroadCast_cmd.cpp BMesh2D_cmd.cpp GetRef_cmd.cpp Mesh2DReader.cpp)
TARGET_LINK_LIBRARIES(Parallel Foundation Geometry TmlMessage TmlType TmlComm Fields ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})
INSTALL(TARGETS Parallel DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
