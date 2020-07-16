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

#include "nt_block.h"

/*!
  construct an empty, uninitialized NeighborTable -> not usable
*/
template<typename T>
NeighborTable<T>::NeighborTable()
  : m_list(),
    m_array(),
    m_idParticleMap(),
    m_p0_global(),
    m_dim(0),
    m_alpha(0),
    m_global_idx(0),
    m_global_idy(0),
    m_global_idz(0),
    m_xsize(0),
    m_ysize(0),
    m_zsize(0),
    m_valid(false)
{
}

/*!
  construct neigbortable with known dimensions

  \param x nr. of grid points in x-direction
  \param y nr. of grid points in y-direction
  \param z nr. of grid points in z-direction
  \param range grid spacing
  \param alpha pair search cutoff
  \param p0_global minimal corner (origin) of the global search space
  \param ix x-index of the local origin 
  \param iy y-index of the local origin
  \param iz z-index of the local origin
*/
template<typename T>
NeighborTable<T>::NeighborTable(
  int x,
  int y,
  int z,
  double range,
  double alpha,
  const Vec3& p0_global,
  int ix,
  int iy,
  int iz
)
  : m_p0_global(p0_global),
    m_dim(range),
    m_alpha(alpha),
    m_xsize(x),
    m_ysize(y),
    m_zsize(z),
    m_valid(true)
{
  m_global_idx=ix;
  m_global_idy=iy;
  m_global_idz=iz;	
  m_min_corner=Vec3(m_p0_global.X()+m_dim*double(m_global_idx),
		    m_p0_global.Y()+m_dim*double(m_global_idy),
		    m_p0_global.Z()+m_dim*double(m_global_idz));
  m_max_corner=Vec3(m_p0_global.X()+m_dim*double(m_global_idx+m_xsize),
		    m_p0_global.Y()+m_dim*double(m_global_idy+m_ysize),
		    m_p0_global.Z()+m_dim*double(m_global_idz+m_zsize));
  console.Debug() << "Ntable corners : " << m_min_corner << " - " << m_max_corner << "\n";
  m_array.resize(m_xsize*m_ysize*m_zsize);
}

/*!
  destruct NeighborTable
*/
template<typename T>
NeighborTable<T>::~NeighborTable()
{
}

/*!
  clean up the search array
*/
template<typename T>
void NeighborTable<T>::clear_search_array()
{
  for(unsigned int iter=0;iter<m_array.size();iter++){
    m_array[iter].erase(m_array[iter].begin(),m_array[iter].end());
  }
}

/*!
  helper function which returns the index in the search array
  into which a given position is mapped. Returns -1 if pos is 
  outside the search space

  \param the position
*/
template<typename T>
int NeighborTable<T>::index(const Vec3& pos)
{
  int idx=-1;
  int ix=int(floor((pos.X()-m_p0_global.X())/m_dim))-m_global_idx;
  int iy=int(floor((pos.Y()-m_p0_global.Y())/m_dim))-m_global_idy;
  int iz=int(floor((pos.Z()-m_p0_global.Z())/m_dim))-m_global_idz;
  if((ix>=0)&&(ix<m_xsize)&&
     (iy>=0)&&(iy<m_ysize)&&
     (iz>=0)&&(iz<m_zsize)){
    idx=m_ysize*m_zsize*ix+m_zsize*iy+iz;
  }
  return idx;
}

/*!
  helper function which returns the index in the search array
  from given x-,y- and z-indices.

  \param x the x-index
  \param y the y-index
  \param z the z-index
  \warning no checks
*/
template<typename T>
int NeighborTable<T>::index(int x,int y,int z) const
{
  return m_ysize*m_zsize*x+m_zsize*y+z;
}

/*!
  check if a position is in the inner part 

  \param pos the position
*/
template<typename T>
bool NeighborTable<T>::isInInner(const Vec3& pos)
{
  bool res;

  int ix=int(floor((pos.X()-m_p0_global.X())/m_dim))-m_global_idx;
  int iy=int(floor((pos.Y()-m_p0_global.Y())/m_dim))-m_global_idy;
  int iz=int(floor((pos.Z()-m_p0_global.Z())/m_dim))-m_global_idz;

  if(m_zsize<3) { // 2D-> ignore Z
    res=((ix>0) && (ix<m_xsize-1) && (iy>0) && (iy<m_ysize-1));
  } else {
    res=((ix>0) && (ix<m_xsize-1) && (iy>0) && (iy<m_ysize-1)&& (iz>0) && (iz<m_zsize-1));
  }

  return res;  
}


/*!
  insert a particle into the NeighborTable

  \param t the particle
*/
template<typename T>
void NeighborTable<T>::insert(const T& t)
{
  // put into list
  typename list<T>::iterator iter=m_list.insert(m_list.end(),t);
  // put into array
  int idx=index(iter->getPos());
  if(idx!=-1){
//     m_idParticleMap.insert(IdParticleMap::value_type(iter->getID(),&(*iter)));
    m_idParticleMap.insert(make_pair(iter->getID(),&(*iter)));
    m_array[idx].push_back(iter);
  } else { // outside -> delete
    typename list<T>::iterator h=iter;
    iter--;
    m_list.erase(h);
  }
}

/*!
  Build or rebuild the search array, calls clean_search_array. 
  Particles outside the search space are removed from the list
*/ 
template<typename T>
void NeighborTable<T>::build()
{
  clear_search_array();
  // clean out id map
  m_idParticleMap.clear();
  for(typename list<T>::iterator iter=m_list.begin();
      iter!=m_list.end();
      iter++){
    int idx=index(iter->getPos());
    int id=iter->getID();
    if(idx!=-1){
      m_array[idx].push_back(iter);
      m_idParticleMap[id]=&(*iter);
    } else { // outside -> delete
      typename list<T>::iterator h=iter;
      iter--;
      // delete from list
      m_list.erase(h);
    }
  }
}

/*!
   Return pointer to particle at index. 
*/ 
template<typename T>
T* NeighborTable<T>::ptr(NeighborTable<T>::indextype idx)
{
  return &(*(m_array[idx.first][idx.second]));
}

/*!
  Return reference to particle at index.
*/
template<typename T>
T& NeighborTable<T>::ref(NeighborTable<T>::indextype idx)
{
  return *(m_array[idx.first][idx.second]);
}

/*!
   Return pointer to particle with given id. Return NULL if the
   table doesn't contain the particle. 

   \param id the id of the particle
*/ 
template<typename T>
T* NeighborTable<T>::ptr_by_id(int id)
{
  T* pParticle = NULL;
  typename IdParticleMap::iterator it = m_idParticleMap.find(id);

  if (it != m_idParticleMap.end()) {
    pParticle = it->second;
  }
  // debug output
  if(pParticle!=NULL){
    if(id!=pParticle->getID()){
      console.Debug() << "inconsistent idParticleMap: " << id << " vs. " << pParticle->getID() << "\n";
    }
  }
  return pParticle;
}

/*!
  Return pointer to particle nearest to given position. Returns
  NULL if position is outside search area.

  \param pos position
*/
template<typename T>
T *NeighborTable<T>::getNearestPtr(const Vec3& pos) 
{
  T* res=NULL;

  // get grid index for pos
  int idx=index(pos);

  double dist=3.0*(m_dim*m_dim); // squared max. dist (grid diagonal) 
  if(idx!=-1){
    for(typename pointtype::iterator iter=m_array[idx].begin();
	iter!=m_array[idx].end();
	iter++){
      double ndist=(pos-(*iter)->getPos()).norm2();
      res=(ndist<dist) ? &(**iter) : res;
      dist=(ndist<dist) ? ndist : dist;
    }
  }

  return res;
}

/*!
  erase particle (not meant to be used directly)
*/
template<typename T>
void NeighborTable<T>::erase(NeighborTable<T>::indextype idx)
{
  m_idParticleMap.erase(m_array[idx.first][idx.second]->getID());
  m_list.erase(m_array[idx.first][idx.second]);
  m_array[idx.first].erase(m_array[idx.first].begin()+idx.second);
}

/*!
  Return representation for a slab of the search array in the xy-plane.

  \param z position of the slab in z-direction
*/
template<typename T>
NTSlab<T> NeighborTable<T>::xy_slab(int z)
{
  return NTSlab<T>(this,DSlice(z,m_ysize,m_zsize,m_xsize,m_ysize*m_zsize));
}

/*!
  Return representation for a slab of the search array in the xz-plane.

  \param y position of the slab in y-direction
*/
template<typename T>
NTSlab<T> NeighborTable<T>::xz_slab(int y)
{
  return NTSlab<T>(this,DSlice(y*m_zsize,m_zsize,1,m_xsize,m_ysize*m_zsize));
}

/*!
  Return representation for a slab of the search array in the yz-plane.

  \param x position of the slab in x-direction
*/
template<typename T>
NTSlab<T> NeighborTable<T>::yz_slab(int x)
{
  return NTSlab<T>(this,DSlice(x*m_ysize*m_zsize,m_zsize,1,m_ysize,m_zsize));
}

/*
  Return block of the search array

  \param xmin minimum index in x-dimension
  \param xmax maximum index in x-dimension
  \param ymin minimum index in y-dimension
  \param ymax maximum index in y-dimension
  \param zmin minimum index in z-dimension
  \param zmax maximum index in z-dimension
*/
template<typename T>
NTBlock<T> NeighborTable<T>::block(int xmin,int xmax,int ymin,int ymax,int zmin,int zmax)
{
  return NTBlock<T>(this,xmin,xmax,ymin,ymax,zmin,zmax);
}

/*
  Return block of the search array.

  \param vmin minimum position
  \param vmax maximum position
*/
template<typename T>
NTBlock<T> NeighborTable<T>::block(const Vec3& vmin,const Vec3& vmax)
{
  // minimum corner index
  int xmin=int(floor((vmin.X()-m_p0_global.X())/m_dim))-m_global_idx;
  int ymin=int(floor((vmin.Y()-m_p0_global.Y())/m_dim))-m_global_idy;
  int zmin=int(floor((vmin.Z()-m_p0_global.Z())/m_dim))-m_global_idz;
  // check for minimum outside 
  xmin=(xmin < 0) ? 0 : xmin;
  ymin=(ymin < 0) ? 0 : ymin;
  zmin=(zmin < 0) ? 0 : zmin;
  
  // maximum corner index
  int xmax=int(floor((vmax.X()-m_p0_global.X())/m_dim))-m_global_idx;
  int ymax=int(floor((vmax.Y()-m_p0_global.Y())/m_dim))-m_global_idy;
  int zmax=int(floor((vmax.Z()-m_p0_global.Z())/m_dim))-m_global_idz;
  // check for maximum outside 
  xmax=(xmax < m_xsize) ? xmax : m_xsize-1;
  ymax=(ymax < m_ysize) ? ymax : m_ysize-1;
  zmax=(zmax < m_zsize) ? zmax : m_zsize-1;
 
  return NTBlock<T>(this,xmin,xmax,ymin,ymax,zmin,zmax);
}



/*!
  Return the inner region, i.e. excluding the boundary, as a NTBlock
*/ 
template<typename T>
NTBlock<T> NeighborTable<T>::inner()
{
  return NTBlock<T>(this,1,m_xsize-2,1,m_ysize-2,1,m_zsize-2);
}



/*!
  Add pairs containing the particles at two given gridpoints to list.

  \param list the pairlist
  \param idx1 the index of the first gridpoint
  \param idx2 the index of the second gridpoint
*/
template<typename T>
void NeighborTable<T>::addPairsToList(T_Handle<pairlist> list,int idx1,int idx2)
{  
  
  for(typename pointtype::iterator iter=m_array[idx1].begin();
      iter!=m_array[idx1].end();
      iter++){
    for(typename pointtype::iterator iter2=m_array[idx2].begin();
	iter2!=m_array[idx2].end();
	iter2++){
      double dist2=((*iter)->getPos()-(*iter2)->getPos()).norm2();
      double dmax=(*iter)->getRad()+(*iter2)->getRad()+m_alpha;
      if(dist2<=(dmax*dmax)){
	if((*iter)->getID()<(*iter2)->getID()){
	  list->push_back(make_pair(&(**iter),&(**iter2)));
	} else {
	  list->push_back(make_pair(&(**iter2),&(**iter)));
	} 
      }
    }
  }
}

/*!
  Add pairs containing the particles at one given gridpoint to list.

  \param list the pairlist
  \param idx the index of the  gridpoint
*/
template<typename T>
void NeighborTable<T>::addPairsToListLocal(T_Handle<pairlist> list,int idx)
{
  if(m_array[idx].size()>=2){ // at least 2 particles here, or no pairs
    for(typename pointtype::iterator iter=m_array[idx].begin();
	iter!=m_array[idx].end()-1;
	iter++){
      for(typename pointtype::iterator iter2=iter+1;
	  iter2!=m_array[idx].end();
	  iter2++){
	double dist2=((*iter)->getPos()-(*iter2)->getPos()).norm2();
	double dmax=(*iter)->getRad()+(*iter2)->getRad()+m_alpha;
	if(dist2<=(dmax*dmax)){
	  if((*iter)->getID()<(*iter2)->getID()){
	    list->push_back(make_pair(&(**iter),&(**iter2)));
	  } else {
	    list->push_back(make_pair(&(**iter2),&(**iter)));
	  } 
	}
      }
    }
  }
}

/*!
  Add pairs containing at least one flagged particle at two 
  given gridpoints to list. 

  \param list the pairlist
  \param idx1 the index of the first gridpoint
  \param idx2 the index of the second gridpoint
*/
template<typename T>
void NeighborTable<T>::addPairsToListFlagged(T_Handle<pairlist> list,int idx1,int idx2)
{  
  
  for(typename pointtype::iterator iter=m_array[idx1].begin();
      iter!=m_array[idx1].end();
      iter++){
    if((*iter)->isFlagged()){
      for(typename pointtype::iterator iter2=m_array[idx2].begin();
	  iter2!=m_array[idx2].end();
	  iter2++){
	if((*iter2)->isFlagged()){
	  if((*iter)->getID()<(*iter2)->getID()){
	    list->push_back(make_pair(&(**iter),&(**iter2)));
	  } else {
	    list->push_back(make_pair(&(**iter2),&(**iter)));
	  } 
	}
      }
    }
  }
}

/*!
  Add pairs containing at least one flagged particle at one 
  given gridpoint to list.

  \param list the pairlist
  \param idx the index of the  gridpoint
*/
template<typename T>
void NeighborTable<T>::addPairsToListLocalFlagged(T_Handle<pairlist> list,int idx)
{
  if(m_array[idx].size()>=2){ // at least 2 particles here, or no pairs
    for(typename pointtype::iterator iter=m_array[idx].begin();
	iter!=m_array[idx].end()-1;
	iter++){
      if((*iter)->isFlagged()){
	for(typename pointtype::iterator iter2=iter+1;
	    iter2!=m_array[idx].end();
	    iter2++){
	  if((*iter2)->isFlagged()){
	    if((*iter)->getID()<(*iter2)->getID()){
	      list->push_back(make_pair(&(**iter),&(**iter2)));
	    } else {
	      list->push_back(make_pair(&(**iter2),&(**iter)));
	    } 
	  }
	}
      }
    }
  }
}


/*!
  Create a full list of pair of neighboring particles and return handle to it.
*/
template<typename T>
T_Handle<typename NeighborTable<T>::pairlist> NeighborTable<T>::getFullList()
{
  T_Handle<typename NeighborTable<T>::pairlist> list=new typename NeighborTable<T>::pairlist;

  for(int ix=0;ix<m_xsize;ix++){
    for(int iy=0;iy<m_ysize;iy++){
      for(int iz=0;iz<m_zsize;iz++){
	// add pairs within the gridpoint
	addPairsToListLocal(list,index(ix,iy,iz));
	// get search range, considering boundaries
	int xmax=(ix<m_xsize-1) ? ix+1 : ix;
	int ymax=(iy<m_ysize-1) ? iy+1 : iy;
	int zmax=(iz<m_zsize-1) ? iz+1 : iz;
	int xmin=(ix>0) ? ix-1 : ix;
	int ymin=(iy>0) ? iy-1 : iy;
	int zmin=(iz>0) ? iz-1 : iz;
	for(int i=xmin;i<=xmax;i++){
	  for(int j=ymin;j<=ymax;j++){
	    for(int k=zmin;k<=zmax;k++){
	      int idx1=index(ix,iy,iz);
	      int idx2=index(i,j,k);
	      if(idx2>idx1){
		addPairsToList(list,idx1,idx2);
	      }
	    }
	  }
	}
      }
    }
  }
  return list;
}

/*!
  Create a list of all pairs of neighboring particles involving a "flagged" particles 
  and return handle to it.
*/
template<typename T>
T_Handle<typename NeighborTable<T>::pairlist> NeighborTable<T>::getNewList()
{
  T_Handle<typename NeighborTable<T>::pairlist> nlist=new typename NeighborTable<T>::pairlist; 

  for(int ix=0;ix<m_xsize;ix++){
    for(int iy=0;iy<m_ysize;iy++){
      for(int iz=0;iz<m_zsize;iz++){
	// add pairs within the gridpoint
	addPairsToListLocalFlagged(nlist,index(ix,iy,iz));
	// get search range, considering boundaries
	int xmax=(ix<m_xsize-1) ? ix+1 : ix;
	int ymax=(iy<m_ysize-1) ? iy+1 : iy;
	int zmax=(iz<m_zsize-1) ? iz+1 : iz;
	int xmin=(ix>0) ? ix-1 : ix;
	int ymin=(iy>0) ? iy-1 : iy;
	int zmin=(iz>0) ? iz-1 : iz;
	for(int i=xmin;i<=xmax;i++){
	  for(int j=ymin;j<=ymax;j++){
	    for(int k=zmin;k<=zmax;k++){
	      int idx1=index(ix,iy,iz);
	      int idx2=index(i,j,k);
	      if(idx2>idx1){
		addPairsToListFlagged(nlist,idx1,idx2);
	      }
	    }
	  }
	}
      }
    }
  }
  for(typename list<T>::iterator iter=m_list.begin();
      iter!=m_list.end();
      iter++){
    iter->setFlag(false);
  }
  return nlist;
}

/*!
  Get list of all particles along a given plane. Naive implementation, i.e. check all particles 
  for distance to plane. The plane is given by one point an the normal. 

  \param orig The origin of the plane
  \param normal The normal of the plane
*/
template<typename T>
T_Handle<typename NeighborTable<T>::particlelist> NeighborTable<T>::getParticlesAtPlane(const Vec3& orig, const Vec3& normal)
{
  T_Handle<typename NeighborTable<T>::particlelist> nlist=new typename NeighborTable<T>::particlelist;

  console.Debug() << "NeighborTable<T>::getParticlesAtPlane: m_dim = " << m_dim << "\n";
  for(typename list<T>::iterator iter=m_list.begin();
      iter!=m_list.end();
      iter++){
    double dist=(iter->getPos()-orig)*normal;
    if(fabs(dist)<m_dim){
      nlist->push_back(&(*iter));
    }
  }
  console.Debug() << "NeighborTable<T>::getParticlesAtPlane: found = " << nlist->size() << "particles\n";

  return nlist;
}

/*!
  Get list of all particles near a given sphere. 
  Naive implementation, i.e. check all particles for distance to sphere.

  \param centre The centre of the sphere
  \param radius The radius of the sphere
*/
template<typename T>
T_Handle<typename NeighborTable<T>::particlelist> NeighborTable<T>::getParticlesNearSphere(const Vec3& centre, const double& radius)
{
  T_Handle<typename NeighborTable<T>::particlelist> nlist=new typename NeighborTable<T>::particlelist;

  console.Debug() << "NeighborTable<T>::getParticlesNearSphere: m_dim = " << m_dim << "\n";
  for(typename list<T>::iterator iter=m_list.begin();
      iter!=m_list.end();
      iter++){
    double dist=(iter->getPos()-centre).norm() - radius;
    if(fabs(dist)<m_dim){
      nlist->push_back(&(*iter));
    }
  }
  console.Debug() << "NeighborTable<T>::getParticlesNearSphere: found = " << nlist->size() << "particles\n";

  return nlist;
}

/*!
  Get list of all particles near a given triangle. Checks all particles at grid within range
  of the bounding box of the triangle; (could be improved)

  \param T the triangle
*/
template<typename T>
T_Handle<typename NeighborTable<T>::particlelist> NeighborTable<T>::getParticlesNearTriangle(const Triangle& Tr)
{
  T_Handle<typename NeighborTable<T>::particlelist> nlist=new typename NeighborTable<T>::particlelist;
  
  // work out search area from bounding box of the triangle
  Vec3 v_min=Tr.getBoundingBoxMin()-Vec3(m_dim,m_dim,m_dim);
  Vec3 v_max=Tr.getBoundingBoxMax()+Vec3(m_dim,m_dim,m_dim);
  // check if completely outside
  if((v_min.X()<m_max_corner.X())&&(v_min.Y()<m_max_corner.Y())&&(v_min.Z()<m_max_corner.Z())&&
     (v_max.X()>m_min_corner.X())&&(v_max.Y()>m_min_corner.Y())&&(v_max.Z()>m_min_corner.Z())){
  
    NTBlock<T> TriangleBlock=block(v_min,v_max);
    for(typename NTBlock<T>::iterator iter=TriangleBlock.begin();
	iter!=TriangleBlock.end();
	iter++){
      if(Tr.sep(iter->getPos()) < iter->getRad()+m_alpha){
	nlist->push_back(&(*iter));
      }
    }
  }
  
  return nlist;
}

/*!
  Get list of all particles near a given edge. Checks all particles at grid within range
  of the bounding box of the triangle; (could be improved - search area via Bresenham)

  \param  E the Edge
*/  
template<typename T>
T_Handle<typename NeighborTable<T>::particlelist> NeighborTable<T>::getParticlesNearEdge(const AEdge* E)
{
  T_Handle<typename NeighborTable<T>::particlelist> nlist=new typename NeighborTable<T>::particlelist;

  // get search area via bounding box
  Vec3 v_min=E->getBoundingBoxMin()-Vec3(m_dim,m_dim,m_dim);
  Vec3 v_max=E->getBoundingBoxMax()+Vec3(m_dim,m_dim,m_dim);
  // check if completely outside
  if((v_min.X()<m_max_corner.X())&&(v_min.Y()<m_max_corner.Y())&&(v_min.Z()<m_max_corner.Z())&&
     (v_max.X()>m_min_corner.X())&&(v_max.Y()>m_min_corner.Y())&&(v_max.Z()>m_min_corner.Z())){
    
    NTBlock<T> SearchBlock=block(v_min,v_max);
       
    for(typename NTBlock<T>::iterator iter=SearchBlock.begin();
	iter!=SearchBlock.end();
	iter++){
      if(E->sep(iter->getPos()) < iter->getRad()+m_alpha){
	nlist->push_back(&(*iter));
      }
    }
  }
  
  return nlist;
}

/*!
  Get list of all particles near a given point. 

  \param  p the point
*/  
template<typename T>
T_Handle<typename NeighborTable<T>::particlelist> NeighborTable<T>::getParticlesNearPoint(const Vec3& p)
{
  T_Handle<typename NeighborTable<T>::particlelist> nlist=new typename NeighborTable<T>::particlelist;

  // get search area via bounding box
  Vec3 v_min=p-Vec3(m_dim,m_dim,m_dim);
  Vec3 v_max=p+Vec3(m_dim,m_dim,m_dim);
  // check if completely outside
  if((v_min.X()<m_max_corner.X())&&(v_min.Y()<m_max_corner.Y())&&(v_min.Z()<m_max_corner.Z())&&
     (v_max.X()>m_min_corner.X())&&(v_max.Y()>m_min_corner.Y())&&(v_max.Z()>m_min_corner.Z())){
    
    NTBlock<T> SearchBlock=block(v_min,v_max);
       
    for(typename NTBlock<T>::iterator iter=SearchBlock.begin();
	iter!=SearchBlock.end();
	iter++){
      if((p-iter->getPos()).norm() < iter->getRad()+m_alpha){
	nlist->push_back(&(*iter));
      }
    }
  }

  return nlist;
}

/*!
  Get list of all particles
*/
template<typename T>
T_Handle<typename NeighborTable<T>::particlelist> NeighborTable<T>::getAllParticles()
{
  T_Handle<typename NeighborTable<T>::particlelist> nlist=new typename NeighborTable<T>::particlelist;
  
  for(typename list<T>::iterator iter=m_list.begin();
      iter!=m_list.end();
      iter++){
    nlist->push_back(&(*iter));
  }

  return nlist;
}

/*!
  output operator for NeighborTable

  \param ost the output stream
  \param NT the NeighborTable
*/
template<typename T>
std::ostream& operator<<(std::ostream& ost, const NeighborTable<T>& NT)
{
  ost << "---NeighborTable---" << endl;
  ost << "3d array dimensions (x,y,z),size : (" << NT.m_xsize << "," << NT.m_ysize<< "," << NT.m_zsize << "), " << NT.m_array.size() << endl;
  ost << "search range : " << NT.m_dim << endl;
  ost << "--list--" << endl;
  ost << (NT.m_list).size() << " elements" << endl;
  for(typename list<T>::const_iterator iter=(NT.m_list).begin();
      iter!=(NT.m_list).end();
      iter++){
    ost << iter->getID() << " " << iter->getPos() << endl;
  }
  ost << "---search array---" << endl;
  for(int ix=0;ix<NT.m_xsize;ix++){
    for(int iy=0;iy<NT.m_ysize;iy++){
      for(int iz=0;iz<NT.m_zsize;iz++){
	unsigned int np=NT.nparts_at_gridpoint(NT.index(ix,iy,iz));
	int idx=NT.index(ix,iy,iz);
	//ost << "---" << endl;
	ost << "(" << ix << "," << iy << "," << iz << ") , [" << idx << "], " << np << " : ";
	for(unsigned int i=0;i<np;i++){
	  ost << (NT.m_array[idx])[i]->getID() << " ";
	}
	ost << endl;
      }
    }
  }
  return ost;
}
