//
// Copyright (c) 2014 CNRS
// Authors: Florent Lamiraux
//
//
// This file is part of hpp-constraints.
// hpp-constraints is free software: you can redistribute it
// and/or modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation, either version
// 3 of the License, or (at your option) any later version.
//
// hpp-constraints is distributed in the hope that it will be
// useful, but WITHOUT ANY WARRANTY; without even the implied warranty
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// General Lesser Public License for more details. You should have
// received a copy of the GNU Lesser General Public License along with
// hpp-constraints. If not, see
// <http://www.gnu.org/licenses/>.

#include <hpp/util/debug.hh>
#include <hpp/constraints/relative-com.hh>
#include <hpp/model/device.hh>
#include <hpp/model/joint.hh>

namespace hpp {
  namespace eigen {
    void convert (const constraints::vector3_t& v, model::vectorOut_t res)
    {
      res [0] = v[0]; res [1] = v[1]; res [2] = v[2];
    }

    void convert (const constraints::matrix3_t& m, matrix3_t& res)
    {
      res (0,0) = m (0,0); res (0,1) = m (0,1); res (0,2) = m (0,2);
      res (1,0) = m (1,0); res (1,1) = m (1,1); res (1,2) = m (1,2);
      res (2,0) = m (2,0); res (2,1) = m (2,1); res (2,2) = m (2,2);
    }
  } // namespace eigen

  namespace constraints {
    
    static vector3_t zero3d (0, 0, 0);

    RelativeComPtr_t RelativeCom::create (const DevicePtr_t& robot,
					  const JointPtr_t& joint,
					  const vector3_t reference, 
                                          std::vector <bool> mask)
    {
      RelativeCom* ptr = new RelativeCom (robot, joint, reference, mask);
      RelativeComPtr_t shPtr (ptr);
      return shPtr;
    }

    RelativeCom::RelativeCom (const DevicePtr_t& robot, const JointPtr_t& joint,
			      const vector3_t reference, std::vector <bool> mask) :
      DifferentiableFunction (robot->configSize (), robot->numberDof (), 
                               3, "RelativeCom"),
      robot_ (robot), joint_ (joint), reference_ (reference), SBT_()
    {
      crossR_.setZero ();
      crossL_.setZero ();
      //-----------Added part -------------------------------mask [0] + mask [1] + mask [2]
      //size_type nbRows = 3; //mask [0] + mask [1] + mask [2];
      //matrix_t selection (nbRows, 3); selection.setZero ();
      //size_type row = 0;
      for (std::size_t i=0; i<3; ++i) 
      {
        for (std::size_t j=0; j<3; ++j) 
           SBT_ (i, j) = 0.0;

	if (mask [i]) 
        {
	  SBT_ (i, i) = 1;
	  //++row;
	}
        //++row;     //leave the matrix squared (just leave in zeros de non-fixed axis row)
      }
      //SBT_ = selection;
    }

    void RelativeCom::impl_compute (vectorOut_t result,
				    ConfigurationIn_t argument)
      const throw ()
    {
      robot_->currentConfiguration (argument);
      robot_->computeForwardKinematics ();
      const Transform3f& Ml = joint_->currentTransformation ();
      const vector3_t& x = robot_->positionCenterOfMass ();
      //fcl::Matrix3f RT = Ml.getRotation (); RT.transpose ();
      //const fcl::Vec3f& tl = Ml.getTranslation ();

      JointPtr_t ra = robot_->getJointByName("RLEG_JOINT5");
      const Transform3f& Mr = ra->currentTransformation ();
      const fcl::Vec3f& pl = Ml.getTranslation ();
      const fcl::Vec3f& pr = Mr.getTranslation ();
      const fcl::Vec3f& center = (pl + pr)*0.5; 
      //const fcl::Vec3f& center(c);
      //const fcl::Vec3f& x_t = SBT_ * (x - tl); 
      //const fcl::Vec3f& c_t = SBT_ * (center - tl);
      fcl::Matrix3f RoT (SBT_);
      eigen::convert (RoT * (x - center), result);  
      //eigen::convert (RT * (RoT * (x - center)), result);                       
      //Original:   
      //eigen::convert (RT * (x - tl) - reference_, result);
      //result = SBT_ * result;
    }

    void RelativeCom::impl_jacobian (matrixOut_t jacobian,
				     ConfigurationIn_t arg) const throw ()
    {
      robot_->currentConfiguration (arg);
      robot_->computeForwardKinematics ();
      const ComJacobian_t& Jcom = robot_->jacobianCenterOfMass ();
      //const JointJacobian_t& Jjoint (joint_->jacobian ());
      //const Transform3f& M = joint_->currentTransformation ();
      //fcl::Matrix3f RT (M.getRotation ()); RT.transpose ();
      const vector3_t& x = robot_->positionCenterOfMass ();
      const JointJacobian_t& Jl (joint_->jacobian ());
      JointPtr_t ra = robot_->getJointByName("RLEG_JOINT5");
      const JointJacobian_t& Jr (ra->jacobian ());

      //const Transform3f& Ml = joint_->currentTransformation ();
      //const Transform3f& Mr = ra->currentTransformation ();
      //const vector3_t& xr (Mr.getTranslation ());
      //const vector3_t& xl (Ml.getTranslation ());
      /*const vector3_t& t (M.getTranslation ());
      cross_ (0,1) = -x [2] + t [2]; cross_ (1,0) = x [2] - t [2];
      cross_ (0,2) = x [1] - t [1]; cross_ (2,0) = -x [1] + t [1];
      cross_ (1,2) = -x [0] + t [0]; cross_ (2,1) = x [0] - t [0];  //*/
      //crossR_ (0,1) = -xr [2];  crossR_ (1,0) = xr [2];
      //crossR_ (0,2) = xr [1];   crossR_ (2,0) = -xr [1];
      //crossR_ (1,2) = -xr [0];  crossR_ (2,1) = xr [0];
      //crossL_ (0,1) = -xl [2];  crossL_ (1,0) = xl [2];
      //crossL_ (0,2) = xl [1];   crossL_ (2,0) = -xl [1];
      //crossL_ (1,2) = -xl [0];  crossL_ (2,1) = xl [0];

      //eigen::matrix3_t eigenRT; eigen::convert (RT, eigenRT);
      eigen::matrix3_t eigenSBT; eigen::convert (SBT_, eigenSBT);
      //eigenRT = eigenSBT * eigenRT;
      jacobian.leftCols (Jl.cols ()) =
	eigenSBT * (Jcom - 0.5* ( Jl.topRows (3) + Jr.topRows (3) ) );
/*	eigenSBT * (Jcom - 0.5* ( (Jl.topRows (3) - crossL_ * Jl.bottomRows (3)) + 
                                 (Jr.topRows (3) - crossR_ * Jr.bottomRows (3)) ) ); //*/
      jacobian.rightCols (jacobian.cols () - Jl.cols ()).setZero ();
      hppDout (info, "Jcom = " << std::endl << Jcom);
      hppDout (info, "Jw = " << std::endl << Jjoint.bottomRows (3));
      hppDout (info, "Jv = " << std::endl << Jjoint.topRows (3));
    }

  } // namespace constraints
} // namespace hpp
