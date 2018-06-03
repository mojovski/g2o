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

#ifndef G2O_EDGE_PROJECT_UV_H_
#define G2O_EDGE_PROJECT_UV_H_

#include "g2o/core/base_binary_edge.h"

#include "vertex_se3.h"
#include "vertex_pointxyz.h"
#include "parameter_camera.h"
#include "g2o_types_slam3d_api.h"
#include <g2o/types/sba/types_six_dof_expmap.h>

namespace g2o {

  /*! \class EdgeSE3PointXYZUV
   * \brief g2o edge from a track to a mono camera. The measuremnt is passed
   * in image coordinates (u,v)
   */
  class G2O_TYPES_SLAM3D_API EdgeSE3PointXYZUV : public BaseBinaryEdge<2, Vector2D, VertexPointXYZ, VertexSE3 > {
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    EdgeSE3PointXYZUV();
    virtual bool read(std::istream& is);
    virtual bool write(std::ostream& os) const;

    // return the error estimate as a 3-vector
    void computeError();
    // jacobian
    virtual void linearizeOplus();


    virtual void setMeasurement(const Vector2D& m){
      _measurement = m;
    }

    virtual bool setMeasurementData(const double* d){
      Eigen::Map<const Vector2D> v(d);
      _measurement = v;
      return true;
    }

    virtual bool getMeasurementData(double* d) const{
      Eigen::Map<Vector2D> v(d);
      v=_measurement;
      return true;
    }

    virtual int measurementDimension() const {return 3;}

    virtual bool setMeasurementFromState() ;

    virtual double initialEstimatePossible(const OptimizableGraph::VertexSet& from,
             OptimizableGraph::Vertex* to) {
      (void) to;
      return (from.count(_vertices[0]) == 1 ? 1.0 : -1.0);
    }

    virtual void initialEstimate(const OptimizableGraph::VertexSet& from, OptimizableGraph::Vertex* to);

  private:
    Eigen::Matrix<double,2,9,Eigen::ColMajor> J; // jacobian before projection

    virtual bool resolveCaches();
    ParameterCamera* params;
    CacheCamera* cache;
  };

}
#endif