#############################################################
##                                                         ##
## Copyright (c) 2003-2017 by The University of Queensland ##
## Centre for Geoscience Computing                         ##
## http://earth.uq.edu.au/centre-geoscience-computing      ##
##                                                         ##
## Primary Business: Brisbane, Queensland, Australia       ##
## Licensed under the Open Software License version 3.0    ##
## http://www.apache.org/licenses/LICENSE-2.0              ##
##                                                         ##
#############################################################

##

include $(top_srcdir)/Config/Cpp.am

lib_LTLIBRARIES=libModel.la

AM_CPPFLAGS += \
      $(SILO_CPPFLAGS) \
    -I$(top_srcdir)/Foundation \
    -I$(top_srcdir)/Parallel   \
    -I$(top_srcdir)/Fields     \
    -I$(top_srcdir)/Parser

libModel_la_SOURCES = \
	ABCDampingIGP.cpp \
	ABCDampingIGP.h \
	AdhesiveFriction.cpp \
	AdhesiveFriction.h \
	BasicParticle.cpp \
	BasicParticle.h \
	BodyForceGroup.cpp \
	BodyForceGroup.h \
	BodyForceGroup.hpp \
	BondedInteraction.cpp \
	BondedInteraction.h \
	BTriangleInteraction.cpp \
	BTriangleInteraction.h \
	BEdge2DInteraction.cpp \
	BEdge2DInteraction.h \
	BCorner2DInteraction.cpp \
	BCorner2DInteraction.h \
	BTriMeshInteraction.h \
	BTriMeshInteractionCpData.h \
	BTriMeshInteractionCpData.cpp \
	BTriMeshIP.h \
	BMesh2DIP.h \
	BWallInteraction.cpp \
	BWallInteraction.h \
	BWallInteraction.hpp \
	BWallInteractionGroup.cpp \
	BWallInteractionGroup.h \
	BWallInteractionGroup.hpp \
	CappedBondedInteraction.cpp \
	CappedBondedInteraction.h \
	DampingIGP.cpp \
	DampingIGP.h \
	Damping.hpp \
	Damping.h \
	ECornerInteraction.cpp \
	ECornerInteraction.h \
	ECorner2DInteraction.cpp \
	ECorner2DInteraction.h \
	EEdgeInteraction.cpp \
	EEdgeInteraction.h \
	EEdge2DInteraction.cpp \
	EEdge2DInteraction.h \
	EMesh2DInteraction.h \
	ElasticInteraction.cpp \
	ElasticInteraction.h \
	ElasticInteractionGroup.cpp \
	ElasticInteractionGroup.h \
	ElasticInteractionGroup.hpp \
	ESphereBodyInteraction.h \
	ESphereBodyInteraction.hpp \
	ESphereBodyInteraction.cpp \
	ESphereBodyInteractionGroup.hpp \
	ESphereBodyInteractionGroup.cpp \
	ETriMeshIP.h \
	ETriMeshInteraction.h \
	ETriangleInteraction.cpp \
	ETriangleInteraction.h \
	EWallInteractionGroup.cpp \
	EWallInteractionGroup.h \
	EWallInteractionGroup.hpp \
	EWallInteraction.h \
	EWallInteraction.hpp \
	FractalFriction.cpp \
	FractalFriction.h \
	FrictionInteraction.cpp \
	FrictionInteraction.h \
	HertzianElasticInteraction.cpp \
	HertzianElasticInteraction.h \
	HertzianViscoElasticFrictionInteraction.cpp \
	HertzianViscoElasticFrictionInteraction.h \
	HertzianViscoElasticInteraction.cpp \
	HertzianViscoElasticInteraction.h \
	HertzMindlinInteraction.cpp \
	HertzMindlinInteraction.h \
	HertzMindlinViscoInteraction.cpp \
	HertzMindlinViscoInteraction.h \
	IGParam.cpp \
	IGParam.h \
	Interaction.cpp \
	Interaction.h \
	InteractionGroup.h \
	InteractionGroup.hpp \
	InteractionParam.cpp \
	InteractionParam.h \
	LinearDashpotInteraction.cpp \
	LinearDashpotInteraction.h \
	LocalDampingIGP.h \
	LocalDampingIGP.cpp \
	LocalDamping.h \
	LocalDamping.hpp \
	MeshData.cpp \
	MeshData.h \
	MeshData2D.cpp \
	MeshData2D.h \
	Mesh2D.cpp \
	Mesh2D.h \
	Mesh2D.hpp \
	Particle.cpp \
	Particle.h \
	RotBondedInteraction.cpp \
	RotBondedInteraction.h \
	RotDamping.hpp \
	RotDamping.h \
	RotElasticInteraction.cpp \
	RotElasticInteraction.h \
	RotElasticInteractionGroup.cpp \
	RotElasticInteractionGroup.h \
	RotElasticInteractionGroup.hpp \
	RotFricInteraction.cpp \
	RotFricInteraction.h \
	RotLocalDamping.h \
	RotLocalDamping.hpp \
	RotPairInteraction.cpp \
	RotPairInteraction.h \
	RotParticle.cpp \
	RotParticle.h \
	ShortBondedInteraction.h \
	ShortBondedInteraction.cpp \
	ShortBondedInteractionCpData.h \
	ShortBondedInteractionCpData.cpp \
	SoftBWallInteractionGroup.cpp \
	SoftBWallInteractionGroup.h \
	SoftBWallInteractionGroup.hpp \
	SoftBWallInteraction.h \
	SoftBWallInteraction.hpp \
	SphereBody.h \
	SphereBody.cpp \
	SphereBodyInteraction.h \
	SphereBodyInteraction.hpp \
	SphereBodyIG.h \
	TempPartStore.h \
	TempPartStore.hpp \
	TriMesh.cpp \
	TriMesh.h \
	TriMesh.hpp \
	ViscWallIG.cpp \
	ViscWallIG.h \
	ViscWallIG.hpp \
	ViscWallInteraction.h \
	ViscWallInteraction.hpp \
	VWFrictionInteraction.cpp \
	VWFrictionInteraction.h \
	Wall.cpp \
	Wall.h \
	WallIG.h \
	WallInteraction.h \
	WallInteraction.hpp \
	BondedInteractionCpData.h \
	BondedInteractionCpData.cpp \
	BMesh2DInteractionCpData.cpp \
	BMesh2DInteractionCpData.h \
	BMesh2DInteraction.h \
	BTriMeshIP.cpp \
	RotParticleVi.cpp \
	RotParticleVi.h \
	RotThermBondedInteraction.cpp \
	RotThermBondedInteraction.h \
	RotThermElasticInteraction.cpp \
	RotThermElasticInteractionGroup.cpp \
	RotThermElasticInteractionGroup.h \
	RotThermElasticInteractionGroup.hpp \
	RotThermElasticInteraction.h \
	RotThermFricInteraction.cpp \
	RotThermFricInteraction.h \
	RotThermPairInteraction.cpp \
	RotThermPairInteraction.h \
	RotThermParticle.cpp \
	RotThermParticle.h \
	TaggedEWallInteractionGroup.h \
	TaggedEWallInteractionGroup.hpp \
	ThermParticle.cpp \
	ThermParticle.h \
	$(top_srcdir)/Fields/field_const.h \
	$(top_srcdir)/Fields/FieldMaster.cpp \
	$(top_srcdir)/Fields/FieldMaster.h \
	$(top_srcdir)/Fields/FieldSlave.cpp \
	$(top_srcdir)/Fields/FieldSlave.h \
	$(top_srcdir)/Fields/VectorWallFieldSlave.h \
	$(top_srcdir)/Fields/VectorWallFieldSlave.hpp \
	$(top_srcdir)/Fields/WallFieldMaster.h \
	$(top_srcdir)/Fields/WallFieldMaster.cpp \
	$(top_srcdir)/Fields/WallFieldSlave.h \
	$(top_srcdir)/Fields/WallFieldSlave.cpp \
	$(top_srcdir)/Parallel/CheckPointable.h \
	$(top_srcdir)/Parallel/CheckPointable.cpp \
	$(top_srcdir)/Parallel/IterativeReader.h \
	$(top_srcdir)/Parallel/IterativeReader.hpp \
	$(top_srcdir)/Parallel/Mesh2DReader.h \
	$(top_srcdir)/Parallel/Mesh2DReader.cpp \
	$(top_srcdir)/Parallel/Mesh2DReaderImpl.cpp \
	$(top_srcdir)/Parallel/MeshReader.cpp \
	$(top_srcdir)/Parallel/MeshReader.h \
	$(top_srcdir)/Parallel/MeshReaderImpl.cpp \
	$(top_srcdir)/Parallel/vec3_mpi.cpp \
	$(top_srcdir)/Parallel/vec3_mpi.h

libModel_la_LIBADD =                                   \
  $(top_builddir)/Foundation/libFoundation.la          \
  $(top_builddir)/Geometry/libGgGeometry.la            \
  $(top_builddir)/tml/message/libTmlMessage.la         \
  $(top_builddir)/tml/type/libTmlType.la               \
  $(top_builddir)/tml/comm/libTmlComm.la               
