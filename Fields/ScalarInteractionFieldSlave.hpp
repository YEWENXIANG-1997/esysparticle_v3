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

//-- STL includes --
#include <vector>
#include <utility>

using std::vector;
using std::pair;

// -- IO includes --
#include <iostream>

using std::cout;
using std::endl;

#include "Foundation/quadtuple.h"
#include "pis/pi_storage.h"

/*!
  constructor 

  \param comm the TML communicator used for sending the data back to the master
  \param pis a pointer to the interaction storage
  \param rdf the particle member function to access the data
*/
template <typename T>
ScalarInteractionFieldSlave<T>::ScalarInteractionFieldSlave(TML_Comm* comm,TParallelInteractionStorage<T>* pis,typename T::ScalarFieldFunction rdf):InteractionFieldSlave<T>(comm,pis)
{
  m_rdf=rdf;
} 

/*!
  send full field data and position of the interaction 
*/
template <typename T>
void ScalarInteractionFieldSlave<T>::SendDataFull()
{
  vector<pair<Vec3,double> > data;

  data=this->m_pis->forAllInnerInteractionsGetWithPos(m_rdf);

  // send data to master
  this->m_comm->send_gather(data,0);
}

/*!
  send full field data and position of the interaction 
*/
template <typename T>
void ScalarInteractionFieldSlave<T>::SendDataFull2()
{
  vector<pair<esys::lsm::quintuple<Vec3,double,Vec3,double,Vec3>,double> > data;

  data=this->m_pis->forAllInnerInteractionsGetRaw2(m_rdf);

  // send data to master
  this->m_comm->send_gather(data,0);
}

/*!
  send full field data, position of the interaction and pos. and Id of the particles. Data format is pid1, pid2, pos1, pos2, ipos, val  
*/
template <typename T>
void ScalarInteractionFieldSlave<T>::SendDataWithPosID()
{
  vector<pair<esys::lsm::quintuple<int,int,Vec3,Vec3,Vec3>, double> > data;

  // debug output 
  console.XDebug() << "ScalarInteractionFieldSlave<T>::SendDataWithPosID()\n";

  data=this->m_pis->forAllInnerInteractionsGetDataWithPosID(m_rdf);

  // debug output 
  console.XDebug() << "sending " << data.size() << " data\n"; 

  // send data to master
  this->m_comm->send_gather(data,0);
}

/*!
  send sum only
*/
template <typename T>
void ScalarInteractionFieldSlave<T>::SendDataSum()
{
  vector<double> data_vec;

  // get data from interactions
  this->m_pis->forAllInnerInteractionsGet(data_vec,m_rdf);

  // sum data
  double sum=0.0;
  for(vector<double>::iterator iter=data_vec.begin();
      iter!=data_vec.end();
      iter++){
    sum+=*iter;
  }

  vector<double> sum_vec;
  sum_vec.push_back(sum);
  this->m_comm->send_gather(sum_vec,0);
}

/*!
  send full field data and position of the interaction 
*/
template <typename T>
void ScalarInteractionFieldSlave<T>::SendDataWithID()
{
  vector<pair<esys::lsm::triplet<int,int,Vec3>, double> > data;

  // debug output 
  console.XDebug() << "ScalarInteractionFieldSlave<T>::SendDataWithID()\n";

  data=this->m_pis->forAllInnerInteractionsGetDataWithID(m_rdf);

  // debug output 
  console.XDebug() << "sending " << data.size() << " data\n"; 

  // send data to master
  this->m_comm->send_gather(data,0);
}


/*!
  send maximum only
*/
template <typename T>
void ScalarInteractionFieldSlave<T>::SendDataMax()
{
  vector<double> data_vec;

  // get data from interactions
  this->m_pis->forAllInnerInteractionsGet(data_vec,m_rdf);

  // sum data
  double max=*(data_vec.begin());
  for(vector<double>::iterator iter=data_vec.begin();
      iter!=data_vec.end();
      iter++){
    max=(*iter > max) ? *iter : max;
  }

  vector<double> max_vec;
  max_vec.push_back(max);
  this->m_comm->send_gather(max_vec,0);
}
