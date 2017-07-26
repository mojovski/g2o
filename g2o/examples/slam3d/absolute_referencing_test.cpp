/* 
This test check whether the absolute referencing feature works with very simple synthetic poses and 3D positions
	
**/

#include <iostream>
#include <cmath>
#include <Eigen/StdVector>
#include <Eigen/Core>

#include <g2o/types/slam3d/vertex_se3.h>
#include <g2o/types/slam3d/vertex_pointxyz.h>
#include <g2o/types/slam3d/edge_se3.h>
#include <g2o/types/slam3d/edge_se3_pointxyz.h>
#include <g2o/types/slam3d/types_slam3d.h>

#include <g2o/core/sparse_optimizer.h>
#include <g2o/core/block_solver.h>
#include <g2o/core/factory.h>
#include <g2o/core/optimization_algorithm_factory.h>
#include <g2o/core/optimization_algorithm_gauss_newton.h>
#include <g2o/solvers/csparse/linear_solver_csparse.h>
#include <g2o/core/eigen_types.h>
#include <g2o/types/slam3d/isometry3d_mappings.h>

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <vector>
#include <math.h>

using namespace std;
using namespace g2o;


int main()
{

	//create some poses SE2(x,y,theta)
	std::vector<SE3Quat> poses;
	Eigen::Quaterniond q(1.0, 0, 0,0);
	for (int i=0; i<5; i++)
	{
		Vector3D t((double)i*0.1, 0, 0);
		SE3Quat pose(q, t);
		poses.push_back(pose);
		cout << "Added pose " << i << ":\n" << pose << endl;
	}
	for (int i=1; i<5; i++)
	{
		Vector3D t(0.4, 0, (double)i*0.1);
		SE3Quat pose(q, t);
		poses.push_back(pose);
		cout << "Added pose " << i << ":\n" << pose << endl;
	}
	EdgeSE3::InformationType information_pose;
	information_pose.setIdentity();
	//information_pose=information_pose*100;
	cout << "information pose: \n" << information_pose << endl;

	//create some landmarks
	std::vector<Eigen::Vector3d> landmarks;
	landmarks.push_back(Eigen::Vector3d(0,0,0));
	cout << "landmark1: " << (*(landmarks.end()-1)).transpose() << endl;

	Eigen::Vector3d l2(0, 0.4, 0);
	//debug

	cout << "landmark2: " << l2.transpose() << endl;
	landmarks.push_back(l2);

	Eigen::Vector3d l3(0, 0.4, 0.4);
	cout << "landmark3: " << l3.transpose() << endl;
	landmarks.push_back(l3);

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
	blockSolver->setSchur(true);
	OptimizationAlgorithmGaussNewton* solver = new OptimizationAlgorithmGaussNewton(blockSolver);


	optimizer.setAlgorithm(solver);


	// Configure the sensor
	// add the parameter representing the sensor offset.
	//--------------------------------
	Eigen::Isometry3d sensorOffsetTransf;
	sensorOffsetTransf.setIdentity();
	ParameterSE3Offset* sensorOffset = new ParameterSE3Offset;
	sensorOffset->setOffset(sensorOffsetTransf);
	sensorOffset->setId(0);
	optimizer.addParameter(sensorOffset);


	// Adding the poses as vertices
	//--------------------------------
	cout << "\nOptimization: Adding robot poses ... \n";
	for (size_t i = 0; i < poses.size(); ++i) {
		const SE3Quat&  p = poses[i];
		VertexSE3* robot =  new VertexSE3;
		robot->setId(i);
		cout << "Adding vertex pose with id " << i << endl;
		robot->setEstimate(p);
		optimizer.addVertex(robot);
	}
	cout << "done ("<< poses.size() << " vertices added)." << endl;


	//----
	// second add the odometry constraints
	cout << "\nOptimization: Adding odometry measurements (edges)... \n";
	for (size_t i = 0; i < poses.size()-1; ++i) {
		const SE3Quat&  p1 = poses[i];
		const SE3Quat&  p2 = poses[i+1];
		SE3Quat p2p1=p1.inverse()*p2; //transforms the points from p2 into p1
		//Isometry3D is declared in g2o_core/eigen_types.h
		Isometry3D rel_pose=internal::fromSE3Quat(p2p1);
		cout << "Rel pose " << i << ": (t,q): (" << p2p1.toVector().transpose() << ")" << endl;

		EdgeSE3* odometry = new EdgeSE3;
		odometry->vertices()[0] = optimizer.vertex(i);
		odometry->vertices()[1] = optimizer.vertex(i+1);
		odometry->setMeasurement(rel_pose);
		odometry->setInformation(information_pose);
		optimizer.addEdge(odometry);
	}
	cout << "done." << endl;


	//----
	// add the landmarks (Ground control points) 

	cout << "\nOptimization: add ground control points (vertices) ... \n";
	for (size_t i = 0; i < landmarks.size(); ++i) {
		const Eigen::Vector3d& l = landmarks[i];
		VertexPointXYZ* landmark = new VertexPointXYZ;
		landmark->setId(i+poses.size());
		std::cout << "Landmark id: " << i+poses.size() << ", xyz: (" << l.transpose() << ") added\n";
		landmark->setEstimate(l);
		optimizer.addVertex(landmark);
	}
	cout << "done." << endl;


	//----
	// add landmark observations
	cout << "Optimization: add landmark observations (edges)... ";
	EdgeSE3PointXYZ::InformationType information_gcp;
	information_gcp.setIdentity();
	//information_gcp=information_gcp*100;

	Eigen::Vector3d measurement_zero(0,0,0); //the robot stand on the landmark.

	EdgeSE3PointXYZ* landmarkObservation =  new EdgeSE3PointXYZ;
	landmarkObservation->vertices()[0] = optimizer.vertex(0);//first pose
	landmarkObservation->vertices()[1] = optimizer.vertex(poses.size()); //connect the very first pose vertex with the first landmark
	landmarkObservation->setMeasurement(measurement_zero); //!!! ADD ALWAYS ZERO, which is the relative measurement of the landmark from the robot pose
	landmarkObservation->setInformation(information_gcp);
	landmarkObservation->setParameterId(0, sensorOffset->id());
	optimizer.addEdge(landmarkObservation);
	std::cout << "\n1. added (" << landmarks[0].transpose() << ")" << endl;
	SE3Quat est=internal::toSE3Quat(dynamic_cast<VertexSE3*>(landmarkObservation->vertices()[0])->estimate());
	cout << "From pose: \n" << est << endl;
	cout << "To PointXYZ:\n" << (dynamic_cast<VertexPointXYZ*>(landmarkObservation->vertices()[1])->estimate().transpose()) << endl;


	EdgeSE3PointXYZ* landmarkObservation2 =  new EdgeSE3PointXYZ;
	landmarkObservation2->vertices()[0] = optimizer.vertex(4); //last pose
	landmarkObservation2->vertices()[1] = optimizer.vertex(poses.size()+1); //second landmark
	landmarkObservation2->setMeasurement(measurement_zero);//!!! ADD ALWAYS ZERO, which is the relative measurement of the landmark from the robot pose
	landmarkObservation2->setInformation(information_gcp);
	landmarkObservation2->setParameterId(0, sensorOffset->id()); //???????????????
	optimizer.addEdge(landmarkObservation2);
	std::cout << "\n2. added (" << landmarks[1].transpose() << "). Done. " << std::endl << std::flush;
	cout << "From pose: \n" << internal::toSE3Quat(dynamic_cast<VertexSE3*>(landmarkObservation2->vertices()[0])->estimate()) << endl;
	cout << "To PointXYZ:\n" << (dynamic_cast<VertexPointXYZ*>(landmarkObservation2->vertices()[1])->estimate().transpose()) << endl;

	EdgeSE3PointXYZ* landmarkObservation3 =  new EdgeSE3PointXYZ;
	landmarkObservation3->vertices()[0] = optimizer.vertex(poses.size()-1);
	landmarkObservation3->vertices()[1] = optimizer.vertex(poses.size()+2);
	landmarkObservation3->setMeasurement(measurement_zero);//!!! ADD ALWAYS ZERO, which is the relative measurement of the landmark from the robot pose
	landmarkObservation3->setInformation(information_gcp);
	landmarkObservation3->setParameterId(0, sensorOffset->id()); //???????????????
	optimizer.addEdge(landmarkObservation3);
	std::cout << "\n3. added third (" << landmarks[2].transpose() << "). Done. " << std::endl << std::flush;
	cout << "From pose: \n" << internal::toSE3Quat(dynamic_cast<VertexSE3*>(landmarkObservation3->vertices()[0])->estimate()) << endl;
	cout << "To PointXYZ:\n" << (dynamic_cast<VertexPointXYZ*>(landmarkObservation3->vertices()[1])->estimate().transpose()) << endl;

	

	/*********************************************************************************
   * optimization
   ********************************************************************************/

  // dump initial state to the disk
  optimizer.save("experiment_absPos3D_before.g2o");

  // prepare and run the optimization
  // fix the first robot pose to account for gauge freedom
  //VertexSE2* firstRobotPose = dynamic_cast<VertexSE2*>(optimizer.vertex(0));
  cout << "Setting vertices with IDs " << poses.size() << ", " << poses.size()+1 << " and " << poses.size()+2 << " to fixed!" << endl;
  dynamic_cast<VertexPointXYZ*>(optimizer.vertex(poses.size()))->setFixed(true);
  dynamic_cast<VertexPointXYZ*>(optimizer.vertex(poses.size()+1))->setFixed(true);
  dynamic_cast<VertexPointXYZ*>(optimizer.vertex(poses.size()+2))->setFixed(true);
  //firstRobotPose->setFixed(true);
  optimizer.setVerbose(true);

  	//----------------------------
	//check the pre-optimization state of the problem
	//----------------------------
	cout << "\n===============================";
	cout << "\n====== Optimization =========\n";
	optimizer.computeActiveErrors();
	SparseOptimizer::EdgeSet::iterator it;
	int ei=0;
	for (it=optimizer.edges().begin(); it!=optimizer.edges().end(); ++it)
	{
		if (ei<poses.size()-1)
		{
			EdgeSE3* e=dynamic_cast<EdgeSE3*>(*it);
			cout << "EdgeSE3 error vector: " << e->error().transpose() << endl;
		} else
		{
			EdgeSE3PointXYZ* e=dynamic_cast<EdgeSE3PointXYZ*>(*it);
			cout << "EdgeSE3PointXYZ error vector: " << e->error().transpose() << endl;
		}
		ei++;
	}



  cerr << "Optimizing" << endl;
  optimizer.initializeOptimization();
  optimizer.optimize(10);
  cerr << "done." << endl;

  cout << "Checking edge errors after optimization:\n";

  optimizer.computeActiveErrors();
  ei=0;
  for (it=optimizer.edges().begin(); it!=optimizer.edges().end(); ++it)
  {
    if (ei<poses.size()-1)
		{
			EdgeSE3* e=dynamic_cast<EdgeSE3*>(*it);
			cout << "EdgeSE3 error vector: " << e->error().transpose() << endl;
		} else
		{
			EdgeSE3PointXYZ* e=dynamic_cast<EdgeSE3PointXYZ*>(*it);
			cout << "EdgeSE3PointXYZ error vector: " << e->error().transpose() << endl;
		}
	ei++;
  }


  optimizer.save("experiment_absPos3D_after.g2o");

  for (int i=0; i<poses.size(); i++)
  {
  	cout << " VertexSE3.pose " << i << ": \n" <<  internal::toSE3Quat(dynamic_cast<VertexSE3*>(optimizer.vertex(i))->estimate()) << endl;
  }

  // freeing the graph memory
  optimizer.clear();


  // destroy all the singletons
  Factory::destroy();

  OptimizationAlgorithmFactory::destroy();

  HyperGraphActionLibrary::destroy();

 //TODO: After all these commands, when the program is cleared up, some constructor fails...
 //Is the problem inside my problem definition?=




}