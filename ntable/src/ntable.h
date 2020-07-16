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

#ifndef __NTABLE_H
#define __NTABLE_H

#include "Foundation/console.h"

//--- project includes ---
#include "Foundation/vec3.h"
#include "ntable/src/dslice.h"
#include "ntable/src/handle.h"
#include "Geometry/Triangle.h"
#include "Geometry/AEdge.h"



//--- STL includes ---
#include <utility>
#include <valarray>
#include <vector>
#include <list>

#if HAVE_CONFIG_H
#include <config.h>
#endif

/* #if HAVE_HASH_MAP */
/* #include <hash_map> */
/* #else */
#include <map>
/* #endif */

using std::valarray;
using std::vector;
using std::list;
using std::pair;
using std::make_pair;

//--- IO includes ---
#include <iostream>

//--- forward declarations ---
template<typename T>
class NTSlab;

template<typename T>
class NTBlock;

template<typename T>
class NTSlab_iter;

/*!
  \class NeighborTable
  \brief class for neighbor search 
*/
template <typename T>
class NeighborTable
{
 public: // types
  typedef valarray<vector<typename list<T>::iterator> > arraytype;
  typedef vector<typename list<T>::iterator>            pointtype;
  typedef pair<int,int>                                 indextype;
  typedef list<pair<T*,T*> >                            pairlist; 
  typedef list<T*>                                      particlelist;
/* #if HAVE_HASH_MAP */
/*   typedef std::hash_map<int, T*>                        IdParticleMap; */
/* #else */
  typedef std::map<int, T*>                             IdParticleMap;
/* #endif */
  typedef typename list<T>::iterator                    iterator;

 private:
  // member variables
  list<T>       m_list;          //!< list of particles
  arraytype     m_array;         //!< search array
  IdParticleMap m_idParticleMap; //!< mapping between particle-id and particle-pointer
  Vec3          m_p0_global;     //!< minimum corner of global search space
  double        m_dim;           //!< grid spacing of search array
  double        m_alpha;         //!< padding factor (dim=2*rmax+alpha)
  int           m_global_idx;    //!< minimum corner index (x component)
  int           m_global_idy;    //!< minimum corner index (y component)
  int           m_global_idz;    //!< minimum corner index (z component)
  int           m_xsize;         //!< number of grid point, x direction
  int           m_ysize;         //!< number of grid point, y direction
  int           m_zsize;         //!< number of grid point, z direction
  bool          m_valid; 
  Vec3          m_min_corner;    //!< minimum corner of the search array
  Vec3          m_max_corner;    //!< maximum corner of the search array
  

  // private member functions
  void clear_search_array();
  int index(const Vec3&); // index from position
 
  void addPairsToList(T_Handle<pairlist>,int,int);
  void addPairsToListLocal(T_Handle<pairlist>,int);
  void addPairsToListFlagged(T_Handle<pairlist>,int,int);
  void addPairsToListLocalFlagged(T_Handle<pairlist>,int);

  arraytype *array(){return &m_array;};

 public:
  //! Constructors
  NeighborTable();
  NeighborTable(int,int,int,double,double,const Vec3&,int,int,int);
  ~NeighborTable();

  void insert(const T&); //!< particle insertion
  void insert(iterator i,const T& data){insert(data);}; //!< STL compat. insert
  void build(); //!< build search array
  inline int index(int,int,int) const; // total index from x,y,z indices
  
  //!< iterators
  iterator begin(){return m_list.begin();};
  iterator end(){return m_list.end();};
 
  //!< number of particles at a given gridpoint 
  unsigned int nparts_at_gridpoint(unsigned int idx) const
  {
    return (m_array[idx]).size();
  }

  // check in posn. is in inner part
  bool isInInner(const Vec3&);

  //!< dimensions
  int xsize(){return m_xsize;};
  int ysize(){return m_ysize;};
  int zsize(){return m_zsize;};
  int size(){return m_list.size();};
  Vec3 base_point() const {return m_p0_global;};
  int base_idx_x() const {return m_global_idx;};
  int base_idx_y() const {return m_global_idy;};
  int base_idx_z() const {return m_global_idz;};

  //!< partial access functions
  NTSlab<T> xy_slab(int);
  NTSlab<T> xz_slab(int);
  NTSlab<T> yz_slab(int);
  NTBlock<T> block(int,int,int,int,int,int);
  NTBlock<T> block(const Vec3&,const Vec3&);
  NTBlock<T> inner();
  
  // access ops
  T* ptr(NeighborTable<T>::indextype);
  T& ref(NeighborTable<T>::indextype);
  T* ptr_by_id(int);
  T *getNearestPtr(const Vec3&); 

  // erase particle (not meant to be used directly)
  void erase(NeighborTable<T>::indextype);

  // neigbor list creation
  T_Handle<pairlist> getFullList();
  T_Handle<pairlist> getNewList();
  T_Handle<particlelist> getParticlesAtPlane(const Vec3&,const Vec3&); // distance as parameter ??
  T_Handle<particlelist> getParticlesNearSphere(const Vec3&,const double&); 
  T_Handle<particlelist> getParticlesNearTriangle(const Triangle&);
  T_Handle<particlelist> getParticlesNearEdge(const AEdge*);
  T_Handle<particlelist> getParticlesNearPoint(const Vec3&);
  T_Handle<particlelist> getAllParticles();

  //!output
  template <typename TT>
  friend std::ostream& operator<<(std::ostream &, const NeighborTable<TT> &);
};

#include "ntable/src/ntable.hpp"

#endif // __NTABLE_H
