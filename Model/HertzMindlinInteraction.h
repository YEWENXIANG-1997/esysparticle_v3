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

#ifndef __HERTZMINDLININTERACTION_H
#define __HERTZMINDLININTERACTION_H

#include "Model/IGParam.h"
#include "Model/Interaction.h"
#include "Model/Particle.h"
#include "Foundation/vec3.h"

#include <iostream>
#include <utility>

/*!
  \brief Interaction group parameters for Hertz-Mindlin interactions
*/
class CHertzMindlinIGP : public AIGParam
{
public:
  CHertzMindlinIGP();

  CHertzMindlinIGP(
    const std::string &name,
    double E,
    double nu,
    double fricCoef,
    double dT
  );

  virtual std::string getTypeString() const
  {
    return "HertzMindlin";
  }

  void setTimeStepSize(double dt);

  double m_E;  // Young's modulus
  double m_nu; // Poisson ratio
  double mu;   // Friction coefficient
  double dt;
};

/*!
   \class CHertzMindlinInteraction
   \brief Frictional+Elastic interaction between particles
   \author Qi Shao & Dion Weatherley
   $Date: 2016-12-14 10:00:00 (wed, 14 dec 2016) $
*/
class CHertzMindlinInteraction : public APairInteraction
{
 public: // types
  typedef CHertzMindlinIGP ParameterType;

  typedef
    double (CHertzMindlinInteraction::* ScalarFieldFunction)(
    ) const;
  typedef
    std::pair<bool,double>
    (CHertzMindlinInteraction::* CheckedScalarFieldFunction)(
    ) const;
  typedef
    Vec3 (CHertzMindlinInteraction::* VectorFieldFunction)(
    ) const;

  static ScalarFieldFunction getScalarFieldFunction(const string&);
  static CheckedScalarFieldFunction getCheckedScalarFieldFunction(
    const string&
  );
  static VectorFieldFunction getVectorFieldFunction(const string&);

 protected:
  double m_E;     //!< Young's modulus
  double m_nu;    //!< Poisson ratio
  double m_r0;    //!< equilibrium distance
  double m_mu;    //!< coefficient of friction
  double m_dt;    //!< time step
  Vec3 m_Ffric;   //!< current frictional force
  Vec3 m_force_deficit; //!< difference between fric. force & force necessary for slip
  Vec3 m_cpos; //!< contact position
  Vec3 m_normal_force; //!< current normal force
  double m_dn; //!< caching displacement for E_pot
  bool m_is_slipping; //!< static/dynamic status of the interaction
  bool m_is_touching; //!< contact status of the interaction
  double m_E_diss; //!< dissipated energy

 public:
  CHertzMindlinInteraction();
  CHertzMindlinInteraction(CParticle*,CParticle*);
  CHertzMindlinInteraction(
    CParticle*,
    CParticle*,
    const CHertzMindlinIGP&
  );
  virtual ~CHertzMindlinInteraction();

  static string getType() {return "CHertzMindlinInteraction";};

  virtual void calcForces();
  virtual bool isPersistent();

  void setTimeStepSize(double dt);

  std::pair<bool,double> getAbsFrictionalForce() const;
  std::pair<bool,double> getAbsFN() const;
  std::pair<bool,double> getAbsMuFN() const;
  std::pair<bool,double> getSlipVelocity() const;
  std::pair<bool,double> getNormalStress() const;
  std::pair<bool,double> getMaxFricStress() const;
  std::pair<bool,double> getAbsFrictionalStress() const;

  double getAbsForceDeficit() const;
  double getPotentialEnergy() const;
  double getSlipping()const;
  double getSticking()const;
  double getDissipatedEnergy() const;
  double Count() const;
  Vec3 getForce() const;
  Vec3 getNormalForce() const;
  virtual Vec3 getPos() const {return m_cpos;};

  std::pair<bool,double> getMuEff(const Vec3&,const Vec3&) const;
  std::pair<bool,double> getMuEffXY() const
  {
    return getMuEff(Vec3(1.0,0.0,0.0),Vec3(0.0,1.0,0.0));
  };
  std::pair<bool,double> getMuEffXZ() const
  {
    return getMuEff(Vec3(1.0,0.0,0.0),Vec3(0.0,0.0,1.0));
  };

  friend std::ostream& operator<<(
    std::ostream&,
    const CHertzMindlinInteraction&
  );
  friend class TML_PackedMessageInterface;
};

#endif //__HERTZMINDLININTERACTION_H
