/////////////////////////////////////////////////////////////
//                                                         //
// Copyright (c) 2003-2017 by The University of Queensland //
// Centre for Geoscience Computing                         //
// http://earth.uq.edu.au/centre-geoscience-computing      //
//                                                         //
// Primary Business: Brisbane, Queensland, Australia       //
// Licensed under the Open Software License version 3.0    //
// http://www.apache.org/licenses/LICENSE-2.0              //
//                                                         //
/////////////////////////////////////////////////////////////

#ifndef __SUBLATTICE_H
#define __SUBLATTICE_H

// -- project includes --

#include "Parallel/mpibuf.h"
#include "Parallel/mpivbuf.h"
#include "Parallel/ASubLattice.h"
#include "Parallel/LatticeParam.h"

#include "Foundation/Timer.h"

#include "Model/TriMesh.h"
#include "Model/Mesh2D.h"

// -- interactions --
#include "Model/Interaction.h"
#include "Model/EWallInteractionGroup.h"
#include "Model/BWallInteractionGroup.h"
#include "Model/ViscWallIG.h"
#include "Model/SoftBWallInteractionGroup.h"

#include "Model/RotBondedInteraction.h"
#include "Model/RotElasticInteraction.h"
#include "Model/RotFricInteraction.h"

#include "Model/SphereBodyIG.h"
#include "Model/ESphereBodyInteractionGroup.h"

// -- fields --
#include "Fields/FieldSlave.h"
#include "Fields/ScalarParticleFieldSlave.h"
#include "Fields/VectorParticleFieldSlave.h"
#include "Fields/ScalarInteractionFieldSlave.h"

//--- TML includes ---
#include "tml/comm/comm.h"
#include "tml/comm/comm_world.h"

//--- MPI includes ---
#include <mpi.h>

//--- STL includes ---
#include <vector>
#include <map>
#include <string>

using std::vector;
using std::map;
using std::string;

template <class T> class ParallelParticleArray;
class AParallelInteractionStorage;

/*!
  \brief Class for a Sublattice.
*/
template <typename T>
class TSubLattice : public ASubLattice
{
 protected:
  typedef ParallelParticleArray<T> ParticleArray;
  // -- particles storage--
  ParticleArray *m_ppa;

  // -- interaction storage --
  map<string,AParallelInteractionStorage*> m_dpis; //! dynamic interactions
  map<string,AParallelInteractionStorage*> m_bpis; //! persistent interactions
  
  typedef map<string,AInteractionGroup<T>*> NameIGroupMap;
  NameIGroupMap m_singleParticleInteractions; //! Body/gravitational force interactions

  // -- special interaction groups --
  map<string,AParallelInteractionStorage*> m_damping;
  map<string,AWallInteractionGroup<T>*> m_WIG;
  map<string,ASphereBodyInteractionGroup<T>*> m_SIG;
  map<string,CWall*> m_walls;
  map<string,CSphereBody*> m_spheres;
  map<string,TriMesh*> m_mesh; 
  map<string,Mesh2D*> m_mesh2d;

  // -- parameters --
  double m_dt; // time step for integration  // sawano:this is not a time step. This is the time which has unit!
  double m_t; // time step // sawano
  double m_nrange; // search range for neighbor table
  //  double m_rmax; // max. particle radius
  double m_alpha; // search range padding

  int m_last_ns; // timestamp of last neighborsearch

  // sawano
  bool m_bPacking; // flag for packing 
  double m_totalVolume; // total volume of all particles
  double m_iniFactor; // initial factor for radius expantion
  double m_cuml_Factor; // cumlation of factor for radius expantion
  double m_beta; // params for radius expansion
  double m_gamma; // params for radius expansion

  // Temporary storage of connections
  map<int, vector<int> > m_temp_conn;

  // -- MPI --
  int m_rank; //!< rank in m_comm
  MPI_Comm m_comm;
  TML_Comm m_tml_comm;
  MPI_Comm m_worker_comm; //!< MPI communicator between workers (excl. master)
  TML_Comm m_tml_worker_comm; //!< TML version of the communicator between workers (excl. master)
  //  int m_num_workers;
  /**
   * Process partitioning, non-zero entries are taken literally,
   * zero entries are computed via MPI_Dims_create.
   */
  esys::lsm::CLatticeParam::ProcessDims m_dims;

  // -- times --
  double packtime; 
  double unpacktime;
  double commtime;
  double forcetime;

  // -- fields --

  // -- "new" fields --
  map<int,AFieldSlave*> m_field_slaves;
  
  MpiWTimers *m_pTimers;

  // -- private functions --
  void calcForces();
  void zeroForces();
  void integrate(double);
  void integrateTherm(double dt);
  void zeroHeat();
  void thermExpansion();
  void calcHeatFrict();
  void calcHeatTrans();
  const double calcFactor(); // sawano
  const double calcTotalVolume(); // sawano

  // functions doing the actual work adding interaction groups
  virtual bool doAddPIG(const string&,const string&,CVarMPIBuffer&,bool tagged=false);
  virtual bool doAddDamping(const string&,CVarMPIBuffer&);

 public:
  TSubLattice(const esys::lsm::CLatticeParam &prm, int rank, MPI_Comm comm, MPI_Comm worker_comm);
  virtual ~TSubLattice();

  virtual void setTimeStepSize(double dt);

  virtual vector<int> getCommCoords() const {return m_ppa->getCommCoords();};
  virtual vector<int> getCommDims() const {return m_ppa->getCommDims();};

  virtual void receiveParticles();
  virtual void receiveConnections();

  virtual const MPI_Comm &getWorkerComm() const {return m_worker_comm;}
  
  virtual void initNeighborTable(const Vec3&,const Vec3&);
  virtual void initNeighborTable(const Vec3&,const Vec3&,const vector<bool>&);
  virtual void addWall();
  virtual void addElasticWIG();
  virtual void addBondedWIG();
  virtual void addDirBondedWIG();
  virtual void addViscWIG();
  virtual void addTaggedElasticWIG();
  virtual void addTriMesh();
  virtual void addTriMeshIG();
  virtual void addBondedTriMeshIG();
  virtual void addMesh2D();
  virtual void addMesh2DIG();
  virtual void addBondedMesh2DIG();
  virtual void addPairIG();
  virtual void addTaggedPairIG();
  virtual void addSingleIG();
  virtual void addBondedIG();
  virtual void addCappedBondedIG();
  virtual void addShortBondedIG();
  virtual void addRotBondedIG();
  virtual void addRotThermBondedIG();
  virtual void addDamping();
  virtual void setExIG();
  virtual void removeIG();
  virtual void setBondBrokenSwitch(); // sawano
  virtual void getWallPos();
  virtual void getWallForce();
  virtual void addSphereBody();
  virtual void addESphereBodyIG();
  virtual void getSphereBodyPos();
  virtual void getSphereBodyForce();

  virtual void oneStep();
  virtual void oneStepTherm();
  virtual void exchangePos();
  void setTimeStep(double dt){m_dt=dt;};

  virtual void rebuildParticleArray();
  virtual void rebuildInteractions();
  virtual void searchNeighbors();
  virtual void updateInteractions();
  virtual void checkNeighbors();
  void resetDisplacements();
  
  virtual void do2dCalculations(bool do2d);
  virtual int getNumParticles();

  std::pair<double, int> findParticleNearestTo(const Vec3 &pt);

  std::pair<int, Vec3> getParticlePosn(int particleId);

  AParallelInteractionStorage *findPis(const std::string &groupName);
  
//  virtual ParticleIdPairVector getBondGroupIdPairs(const std::string &groupName);

  template <typename TmplPis>
  ParticleIdPairVector getBondGroupIdPairs(TmplPis *pBondPis);

  virtual void moveParticleTo();
  virtual void moveTaggedParticlesBy();
  virtual void moveSingleParticleTo(int particleId, const Vec3 &posn);
  virtual void moveWallBy();
  virtual void moveSphereBodyBy();
  virtual void setWallNormal();
  virtual void applyForceToWall();
  virtual void setParticleVelocity();
  virtual void setParticleFluidForce(); // sawano
  virtual void setParticleRadiusFactor(); // sawano
  virtual void setRadiusExpansionParams(); // sawano
  virtual void setParticleRadiusInitFactor(); // sawano
  virtual void setParticleRadiusFactor_inner(const double fac); // sawano
  virtual void setFlagforRadiusExpansion(); // sawano
  virtual void setParticleTag(); // sawano
  virtual void getDistance(); // sawano
  virtual void getTotalVolume(); // sawano
  virtual void setVelocityOfWall();
  virtual void setParticleNonDynamic();
  virtual void setParticleNonRot();
  virtual void setParticleNonTrans();
  virtual void setParticleDensity();
  virtual void setTaggedParticleVel();
  virtual void tagParticleNearestTo();
  virtual void moveSingleNode();
  virtual void moveTaggedNodes();
  virtual void translateMeshBy(const std::string &meshName, const Vec3 &translation);

  virtual void sendDataToMaster();
  virtual void countParticles();

  virtual void printStruct();
  virtual void printData();
  virtual void printTimes();

  //-- "new" field functions --
  virtual void addScalarParticleField();
  virtual void addVectorParticleField();
  virtual void addScalarInteractionField();
  virtual void addScalarHistoryInteractionField();
  virtual void addVectorInteractionField();
  virtual void addScalarTriangleField();
  virtual void addVectorTriangleField();
  virtual void sendFieldData();
  virtual void addVectorWallField();

  // -- checkpointing --
  virtual void saveSnapShotData(std::ostream &oStream);
  virtual void saveCheckPointData(std::ostream &oStream);
  virtual void loadCheckPointData(std::istream &iStream);
  
  // -- timer --
  virtual void setTimer(MpiWTimers &timers) {m_pTimers = &timers;}

  // -- mesh data exchange --
  virtual void getMeshNodeRef();
  virtual void getMeshFaceRef();
  virtual void getMesh2DStress();
  virtual void getTriMeshForce();
  virtual void getParticleData(const IdVector &particleIdVector);
};

#include "Parallel/SubLattice.hpp"

#endif
