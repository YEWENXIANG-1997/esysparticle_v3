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


#include "tml/message/packed_message_interface.h"

//----
// partial specialisations of the pack/unpack functions
// for simple data types.

//--- int ---
template<>
void TML_PackedMessageInterface::pack<int>(const int& i)
{
  append(i);
}

template<>
void TML_PackedMessageInterface::unpack<int>(int& i)
{
  i=pop_int();
}

//--- double ----
template<>
void TML_PackedMessageInterface::pack<double>(const double& d)
{
  append(d);
}

template<>
void TML_PackedMessageInterface::unpack<double>(double& d)
{
  d=pop_double();
}


/*!
  Pack a Vec3 into a TML packed message

  \param p the Vec3
*/
template<>
void TML_PackedMessageInterface::pack<Vec3>(const Vec3& v)
{
  append(v.X());
  append(v.Y());
  append(v.Z());
}

/*!
  Unpack a Vec3 from a TML packed message

  \param p the Vec3
*/
template<>
void TML_PackedMessageInterface::unpack<Vec3>(Vec3& v)
{
  double db[3];

  pop_doubles(db,3);
  v=Vec3(db[0],db[1],db[2]);
}
  
