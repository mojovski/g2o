/*
This experiment aims to test the ability of g2o
to start relative odometry process and
to insert some absolute, uncertain measurements
which must help to register the odometry in 
some absolute coordinate frame.
*/

#include <iostream>
#include <cmath>
#include <Eigen/StdVector>
#include "simulator.h"

#include "vertex_se2.h"
#include "vertex_point_xy.h"
#include "edge_se2.h"
#include "edge_se2_pointxy.h"
#include "types_tutorial_slam2d.h"

#include "g2o/core/sparse_optimizer.h"
#include "g2o/core/block_solver.h"
#include "g2o/core/factory.h"
#include "g2o/core/optimization_algorithm_factory.h"
#include "g2o/core/optimization_algorithm_gauss_newton.h"
#include "g2o/solvers/csparse/linear_solver_csparse.h"

using namespace std;
using namespace g2o;
using namespace g2o::tutorial;

int main()
{
  //create some poses SE2(x,y,theta)
  std::vector<SE2> poses;
  poses.push_back(SE2(0,0,0));
  poses.push_back(SE2(0.1,0,0));
  poses.push_back(SE2(0.1,0.1,0));

  //create some odometry
  Eigen::Matrix3d odometry_covariance;
  odometry_covariance.fill(0.);
  odometry_covariance(0, 0) = 0.00100;
  odometry_covariance(1, 1) = 0.00100;
  odometry_covariance(2, 2) = 0.00100;
  Eigen::Matrix3d information = odometry_covariance.inverse();
  Eigen::Matrix2d information2d;
  information2d.fill(0.);
  information2d(0,0)=100;
  information2d(1,1)=100;


  std::vector<Simulator::GridEdge> odometries;
  odometries.push_back(Simulator::GridEdge());
  Simulator::GridEdge& edge=odometries.back();
  edge.from=0;
  edge.to=1;
  edge.trueTransf=poses[0].inverse()*poses[1];
  edge.simulatorTransf=edge.trueTransf;
  edge.information=information;

  odometries.push_back(Simulator::GridEdge());
  Simulator::GridEdge& edge2=odometries.back();
  edge2.from=1;
  edge2.to=2;
  edge2.trueTransf=poses[1].inverse()*poses[2];
  edge2.simulatorTransf=edge2.trueTransf;
  edge2.information=information;

  //Add the absolute positions
  std::vector<Eigen::Vector2d> landmarks;
  landmarks.push_back(Eigen::Vector2d(0,-1));
  landmarks.push_back(Eigen::Vector2d(-1,0));

   /*********************************************************************************
   * creating the optimization problem
   ********************************************************************************/

  typedef BlockSolver< BlockSolverTraits<-1, -1> >  SlamBlockSolver;
  typedef LinearSolverCSparse<SlamBlockSolver::PoseMatrixType> SlamLinearSolver;

  // allocating the optimizer
  SparseOptimizer optimizer;
  SlamLinearSolver* linearSolver = new SlamLinearSolver();
  linearSolver->setBlockOrdering(false);
  SlamBlockSolver* blockSolver = new SlamBlockSolver(linearSolver);
  OptimizationAlgorithmGaussNewton* solver = new OptimizationAlgorithmGaussNewton(blockSolver);

  optimizer.setAlgorithm(solver);

  // add the parameter representing the sensor offset. What is this????
  SE2 sensorOffsetTransf(0,0,0);
  ParameterSE2Offset* sensorOffset = new ParameterSE2Offset;
  sensorOffset->setOffset(sensorOffsetTransf);
  sensorOffset->setId(0);
  optimizer.addParameter(sensorOffset);

  // adding the odometry to the optimizer
  // first adding all the vertices
  cerr << "Optimization: Adding robot poses ... ";
  for (size_t i = 0; i < poses.size(); ++i) {
    const SE2&  p = poses[i];
    VertexSE2* robot =  new VertexSE2;
    robot->setId(i);
    robot->setEstimate(p);
    optimizer.addVertex(robot);
  }
  cerr << "done." << endl;

  // second add the odometry constraints
  cerr << "Optimization: Adding odometry measurements ... ";
  for (size_t i = 0; i < odometries.size(); ++i) {
    const Simulator::GridEdge& simEdge = odometries[i];

    EdgeSE2* odometry = new EdgeSE2;
    odometry->vertices()[0] = optimizer.vertex(simEdge.from);
    odometry->vertices()[1] = optimizer.vertex(simEdge.to);
    odometry->setMeasurement(simEdge.simulatorTransf);
    odometry->setInformation(simEdge.information);
    optimizer.addEdge(odometry);
  }
  cerr << "done." << endl;

  // add the landmark observations
  cerr << "Optimization: add landmark vertices ... \n";
  for (size_t i = 0; i < landmarks.size(); ++i) {
    const Eigen::Vector2d& l = landmarks[i];
    VertexPointXY* landmark = new VertexPointXY;
    landmark->setId(i+poses.size());
    std::cout << "Landmark id: " << i+poses.size() << " added\n";
    landmark->setEstimate(l);
    optimizer.addVertex(landmark);
  }
  cerr << "done." << endl;

  cerr << "Optimization: add landmark observations ... ";
  Eigen::Vector2d measurement_zero(0,0); //the robot stand on the landmark.
  EdgeSE2PointXY* landmarkObservation =  new EdgeSE2PointXY;
  landmarkObservation->vertices()[0] = optimizer.vertex(0);
  landmarkObservation->vertices()[1] = optimizer.vertex(poses.size());
  landmarkObservation->setMeasurement(measurement_zero);//(landmarks[0]);
  landmarkObservation->setInformation(information2d);
  landmarkObservation->setParameterId(0, sensorOffset->id());
  optimizer.addEdge(landmarkObservation);
  std::cout << "Added first " << std::endl << std::flush;

  EdgeSE2PointXY* landmarkObservation2 =  new EdgeSE2PointXY;
  landmarkObservation2->vertices()[0] = optimizer.vertex(2);
  landmarkObservation2->vertices()[1] = optimizer.vertex(poses.size()+1);
  landmarkObservation2->setMeasurement(measurement_zero);
  landmarkObservation2->setInformation(information2d);
  landmarkObservation2->setParameterId(0, sensorOffset->id()); //???????????????
  optimizer.addEdge(landmarkObservation2);

  cerr << "done." << endl;


  /*********************************************************************************
   * optimization
   ********************************************************************************/

  // dump initial state to the disk
  optimizer.save("experiment_before.g2o");

  // prepare and run the optimization
  // fix the first robot pose to account for gauge freedom
  //VertexSE2* firstRobotPose = dynamic_cast<VertexSE2*>(optimizer.vertex(0));
  dynamic_cast<VertexPointXY*>(optimizer.vertex(3))->setFixed(true);
  dynamic_cast<VertexPointXY*>(optimizer.vertex(4))->setFixed(true);
  //firstRobotPose->setFixed(true);
  optimizer.setVerbose(true);

  cerr << "Optimizing" << endl;
  optimizer.initializeOptimization();
  optimizer.optimize(10);
  cerr << "done." << endl;

  optimizer.save("experiment_after.g2o");

  std::cout << dynamic_cast<VertexSE2*>(optimizer.vertex(0)) << std::endl;

  // freeing the graph memory
  optimizer.clear();

  // destroy all the singletons
  Factory::destroy();
  OptimizationAlgorithmFactory::destroy();
  HyperGraphActionLibrary::destroy();

  return 0;
}
