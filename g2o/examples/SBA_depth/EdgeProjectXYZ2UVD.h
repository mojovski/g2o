#ifndef EdgeProjectXYZ2UVD_H
#define EdgeProjectXYZ2UVD_H

/*
NOT READY YET!!
Trying the EdgeSE3PointXYZDepth
first!!
*/

#include "g2o/core/base_vertex.h"
#include "g2o/core/base_binary_edge.h"
#include "g2o/core/base_unary_edge.h"
#include "g2o/types/slam3d/se3_ops.h"
#include "types_sba.h"
#include <Eigen/Geometry>
#include <g2o/types/sba/types_6dof_expmap.h>

namespace g2o {


/*@brief
Projects  a 3d point (xyz) to the image plane with also a depth value 
**/
  class G2O_TYPES_SBA_API EdgeProjectXYZ2UVD : public  g2o::BaseMultiEdge<2, Vector2D>
  {
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    EdgeProjectXYZ2UVD()  {
      resizeParameters(1);
      installParameter(_cam, 0);
    }

    virtual bool read  (std::istream& is);
    virtual bool write (std::ostream& os) const;
    void computeError  ();
    virtual void linearizeOplus ();
    CameraParameters * _cam;
  };


}

#endif