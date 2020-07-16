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

#ifndef __ROTBONDEDINTERACTION_H
#define __ROTBONDEDINTERACTION_H

// -- project includes --
#include "Model/RotPairInteraction.h"
#include "Model/RotParticle.h"
#include "Model/BondedInteractionCpData.h"
#include "Model/IGParam.h"
#include "Foundation/vec3.h"

// -- I/O includes --
#include <iostream>
using std::ostream;

/*!
  \struct CRotBondedIGP
  \brief Interaction parameters for bonded interaction between rotational particles
  \author Shane Latham, Steffen Abe
  $Revision$
  $Date$
*/

 double calc_angle( double , double ) ;

class CRotBondedIGP : public AIGParam
{
public:
  CRotBondedIGP();
  CRotBondedIGP(
    const  std::string &name,
    double kr,
    double ks,
    double kt,
    double kb,
    double max_nForce,
    double max_shForce,
    double max_tMoment,
    double max_bMoment,
    int    tag,
    bool   scaling,
    bool   AmeanR_scaling,
    double truncated,
    double beta1,
    double beta2
  );

  CRotBondedIGP(
    const  std::string &name,
    double youngsModulus,
    double poissonsRatio,
    double cohesion,
    double tanAngle,
    int    tag,
    bool   AmeanR_scaling,
    double truncated,
    double beta1,
    double beta2
  );

  virtual std::string getTypeString() const
  {
    return "RotBonded";
  }
  
  double kr,ks,kt,kb ;
  double max_nForce, max_shForce,max_tMoment, max_bMoment;
  int    tag;
  bool   scaling;
  bool   meanR_scaling;
  double truncated;
  double beta1;
  double beta2;
};

/*!
   \class CRotBondedInteraction
   \brief Elastic interaction between bonded particles between rotational particles
   \author Shane Latham, Steffen Abe
   $Revision$
   $Date$
*/
class CRotBondedInteraction : public ARotPairInteraction
{
 public: // types
  typedef CRotBondedIGP ParameterType;

  /**
   * Used by PIS to save/load check-point data for objects of this type.
   */
  typedef BondedInteractionCpData CheckPointable;

  typedef double (CRotBondedInteraction::* ScalarFieldFunction)() const; 
  typedef pair<bool,double> (CRotBondedInteraction::* CheckedScalarFieldFunction)() const;
  typedef Vec3 (CRotBondedInteraction::* VectorFieldFunction)() const; 

 private:

  //   protected:
  double m_dist;  //!< current distance, cached from last calcForces()
  double m_r0;  //!< equilibrium separation

  double m_kr ;     //!< spring constant
  double m_ks ;
  double m_kb ;
  double m_kt ;

  double m_max_nForce;   // always >0
  double m_max_shForce ;
  double m_max_tMoment ;
  double m_max_bMoment ;

  double m_nForce;  //  >0, pulling; <0 , compressing
  double m_shForce ; // always >0
  double m_tMoment ;
  double m_bMoment ;

  Vec3 m_force;   //!< current force, cached for E_pot calculation
  Vec3 m_moment ;

  Vec3 m_cpos; // ?
  Vec3 m_D; //!< initial positions of the particles
  int  m_tag;
  bool m_scaling;
  bool m_meanR_scaling;
  double m_truncated;
  double m_beta1;
  double m_beta2;

 public:

  CRotBondedInteraction();
  CRotBondedInteraction(CRotParticle*,CRotParticle*,const CRotBondedIGP&);
  virtual ~CRotBondedInteraction();
                                                                                
  static ScalarFieldFunction getScalarFieldFunction(const string&);
  static CheckedScalarFieldFunction getCheckedScalarFieldFunction(const string&);
  static VectorFieldFunction getVectorFieldFunction(const string&);

  static string getType(){return "RotBonded";};

  int getTag() const;
  void setTag(int tag);

  void calcForces();
  //void setBreak(double);
  bool broken();

  double getPotentialEnergy() const;
  double getNormalPotentialEnergy() const;
  double getShearPotentialEnergy() const;
  double getTwistPotentialEnergy() const;
  double getBendPotentialEnergy() const;
  double getCriterion() const;
  Vec3   getForce() const;
  Vec3   getNormalForce() const;
  Vec3   getTangentialForce() const;
  virtual Vec3 getPos() const {return m_cpos;};

  Vec3 getCentrePtDiff() const;
  Vec3 getInitialCentrePtDiff() const;
  Vec3 getInitialMidPoint() const;

  Vec3 getP2ShearForcePt() const;
  Vec3 getP1ShearForcePt() const;
  Vec3 getContactPoint() const;

  Vec3 getShearDiff() const;

  friend ostream& operator<<(ostream&,const CRotBondedInteraction&);
  friend class TML_PackedMessageInterface;
  
  virtual void saveCheckPointData(std::ostream &oStream);

  virtual void loadCheckPointData(std::istream &iStream);

  // save/load of restart parameters
  virtual void saveRestartData(std::ostream &oStream);
  virtual void loadRestartData(std::istream &iStream);
};

#endif //__BONDEDINTERACTION_H
