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


#include "Tools/StressCalculator/Contact.h"

namespace esys
{
  namespace lsm
  {
    bool ParticleData::s_is3d = true;
    const double ParticleData::FOUR_THIRDS_PI = ((4.0/3.0)*M_PI); // sawano
    //const double ParticleData::FOUR_THIRDS_PI = ((4.0/3.0)*M_PIl);
  }
}

