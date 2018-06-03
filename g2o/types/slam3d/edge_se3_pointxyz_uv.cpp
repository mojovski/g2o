// g2o - General Graph Optimization
// Copyright (C) 2011 R. Kuemmerle, G. Grisetti, W. Burgard
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// * Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright
//   notice, this list of conditions and the following disclaimer in the
//   documentation and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
// IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
// TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
// PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
// TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "edge_se3_pointxyz_uv.h"

namespace g2o {
  using namespace g2o;

  // point to camera projection, monocular
  EdgeSE3PointXYZUV::EdgeSE3PointXYZUV() : BaseBinaryEdge<2, Vector2D/*Measurement of the point IN the image (u,v,depth)*/, VertexPointXYZ, VertexSE3>() {
    resizeParameters(1);
    installParameter(params, 0);
    information().setIdentity();
    /*information()(2,2)=100; //Whats this?? The depth uncertainty is high??
    J.fill(0);
    J.block<3,3>(0,0) = -Matrix3D::Identity();
     */
  }

  bool EdgeSE3PointXYZUV::resolveCaches(){
    ParameterVector pv(1);
    pv[0]=params;
    resolveCache(cache, (OptimizableGraph::Vertex*)_vertices[1],"CACHE_CAMERA",pv);
    return cache != 0;
  }




  bool EdgeSE3PointXYZUV::read(std::istream& is) {
    int pid;
    is >> pid;
    setParameterId(0,pid);

    // measured keypoint
    Vector2D meas;
    for (int i=0; i<2; i++) is >> meas[i];
    setMeasurement(meas);
    // don't need this if we don't use it in error calculation (???)
    // information matrix is the identity for features, could be changed to allow arbitrary covariances    
    if (is.bad()) {
      return false;
    }
    for ( int i=0; i<information().rows() && is.good(); i++)
      for (int j=i; j<information().cols() && is.good(); j++){
  is >> information()(i,j);
  if (i!=j)
    information()(j,i)=information()(i,j);
      }
    if (is.bad()) {
      //  we overwrite the information matrix
      information().setIdentity();
      //information()(2,2)=10/_measurement(2); // scale the info by the inverse of the measured depth
    } 
    return true;
  }

  bool EdgeSE3PointXYZUV::write(std::ostream& os) const {
    os << params->id() << " ";
    for (int i=0; i<2; i++) os  << measurement()[i] << " ";
    for (int i=0; i<information().rows(); i++)
      for (int j=i; j<information().cols(); j++) {
        os <<  information()(i,j) << " ";
      }
    return os.good();
  }


  void EdgeSE3PointXYZUV::computeError() {
    // from cam to point (track)
    //VertexSE3 *cam = static_cast<VertexSE3*>(_vertices[0]);
    VertexPointXYZ *point = static_cast<VertexPointXYZ*>(_vertices[0]);

	g2o::Affine3D w2i = cache->w2i(); //world to image plane transform
    Vector3D p = w2i * point->estimate();
    Vector3D perr;
    perr = p/p(2);
    //perr(2) = 1; //p(2);

    // error, which is backwards from the normal observed - calculated
    // _measurement is the measured 3d point xyz
    _error = perr.head<2>() - _measurement;
    //    std::cout << _error << std::endl << std::endl;
  }

  void EdgeSE3PointXYZUV::linearizeOplus() {
      //ef: copied from types_six_dof_exmap.cpp,  EdgeProjectXYZ2UV::linearizeOplus()
      VertexSE3 * vj = static_cast<VertexSE3 *>(_vertices[1]);
      SE3Quat T(vj->estimate().rotation(), vj->estimate().translation());

      VertexPointXYZ* vi = static_cast<VertexPointXYZ*>(_vertices[0]);
      Vector3D xyz = vi->estimate();
      Vector3D xyz_trans = T.map(xyz);

      double x = xyz_trans[0];
      double y = xyz_trans[1];
      double z = xyz_trans[2];
      double z_2 = z*z;

      const CameraParameters * cam = static_cast<const CameraParameters *>(parameter(0));

      Eigen::Matrix<double,2,3,Eigen::ColMajor> tmp;
      tmp(0,0) = cam->focal_length;
      tmp(0,1) = 0;
      tmp(0,2) = -x/z*cam->focal_length;

      tmp(1,0) = 0;
      tmp(1,1) = cam->focal_length;
      tmp(1,2) = -y/z*cam->focal_length;

      Eigen::Matrix<double,2,3> vtmp=-1./z * tmp * T.rotation().toRotationMatrix();
      _jacobianOplusXi =  vtmp;

      _jacobianOplusXj(0,0) =  x*y/z_2 *cam->focal_length;
      _jacobianOplusXj(0,1) = -(1+(x*x/z_2)) *cam->focal_length;
      _jacobianOplusXj(0,2) = y/z *cam->focal_length;
      _jacobianOplusXj(0,3) = -1./z *cam->focal_length;
      _jacobianOplusXj(0,4) = 0;
      _jacobianOplusXj(0,5) = x/z_2 *cam->focal_length;

      _jacobianOplusXj(1,0) = (1+y*y/z_2) *cam->focal_length;
      _jacobianOplusXj(1,1) = -x*y/z_2 *cam->focal_length;
      _jacobianOplusXj(1,2) = -x/z *cam->focal_length;
      _jacobianOplusXj(1,3) = 0;
      _jacobianOplusXj(1,4) = -1./z *cam->focal_length;
      _jacobianOplusXj(1,5) = y/z_2 *cam->focal_length;

    /*//VertexSE3 *cam = static_cast<VertexSE3 *>(_vertices[0]);
    VertexPointXYZ *vp = static_cast<VertexPointXYZ *>(_vertices[1]);

    const Vector3D& pt;
    pt.head<2>()= vp->estimate();
    pt(2)=1;

    Vector3D Zcam = cache->w2l() * pt;

    //  J(0,3) = -0.0;
    J(0,4) = -2*Zcam(2);
    J(0,5) = 2*Zcam(1);

    J(1,3) = 2*Zcam(2);
    //  J(1,4) = -0.0;
    J(1,5) = -2*Zcam(0);

    J(2,3) = -2*Zcam(1);
    J(2,4) = 2*Zcam(0);
    //  J(2,5) = -0.0;

    J.block<3,3>(0,6) = cache->w2l().rotation();

    Eigen::Matrix<double,3,9,Eigen::ColMajor> Jprime = params->Kcam_inverseOffsetR()  * J;
    Vector3D Zprime = cache->w2i() * pt;

    Eigen::Matrix<double,3,9,Eigen::ColMajor> Jhom;
    Jhom.block<2,9>(0,0) = 1/(Zprime(2)*Zprime(2)) * (Jprime.block<2,9>(0,0)*Zprime(2) - Zprime.head<2>() * Jprime.block<1,9>(2,0));
    Jhom.block<1,9>(2,0) = Jprime.block<1,9>(2,0);

    _jacobianOplusXi = Jhom.block<3,6>(0,0);
    _jacobianOplusXj = Jhom.block<3,3>(0,6);
     */
  }


  bool EdgeSE3PointXYZUV::setMeasurementFromState(){
    //VertexSE3 *cam = static_cast<VertexSE3*>(_vertices[0]);
    VertexPointXYZ *point = static_cast<VertexPointXYZ*>(_vertices[0]);

    // calculate the projection
    const Vector3D& pt = point->estimate();

    Vector3D p = cache->w2i() * pt;
    Vector3D perr;
    perr = p/p(2);
    //perr(2) = p(2);
    _measurement = perr.head<2>();
    return true;
  }


  void EdgeSE3PointXYZUV::initialEstimate(const OptimizableGraph::VertexSet& from, OptimizableGraph::Vertex* /*to_*/)
  {
    (void) from;
    assert(from.size() == 1 && from.count(_vertices[0]) == 1 && "Can not initialize VertexDepthCam position by VertexTrackXYZ");

    VertexSE3 *cam = dynamic_cast<VertexSE3*>(_vertices[0]);
    VertexPointXYZ *point = dynamic_cast<VertexPointXYZ*>(_vertices[1]);
    const Eigen::Matrix<double, 3, 3, Eigen::ColMajor>& invKcam = params->invKcam();
    Vector3D p;
    p(2) = _measurement(2);
    p.head<2>() = _measurement.head<2>()*p(2);
    p=invKcam*p;
    point->setEstimate(cam->estimate() * (params->offset() * p));
  }

}