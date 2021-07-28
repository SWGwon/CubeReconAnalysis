#ifndef __ANALYSISeXAMPLE__
#define __ANALYSISeXAMPLE__

#include <unistd.h>

#include <CubeEvent.hxx>
#include <CubeG4Hit.hxx>
#include <CubeG4Trajectory.hxx>
#include <CubeAlgorithmResult.hxx>
#include <CubeReconTrack.hxx>
#include <CubeReconCluster.hxx>
#include <CubeReconNode.hxx>
#include <CubeHit.hxx>
#include <CubeInfo.hxx>

#include <ToolPrimaryId.hxx>
#include <ToolG4Hits.hxx>
#include <ToolMainTrajectory.hxx>
#include <ToolTrueDirection.hxx>
#include <ToolContained.hxx>
#include <ToolRecon.hxx>

#include "TFile.h"
#include "TTree.h"

int eventNum = 0;
bool SPECIFYEVENT = false;
Cube::Event* event = nullptr;
TFile* inputFile = nullptr;
TTree* inputTree = nullptr;

bool ParseArgs(int argc, char* argv[]);
void PrintSyntax();
void SetBranchAddressInputTree();
void Analysis();
void PrintEventInformation(int eventNum);
void PrintPrimaryInformation(const Cube::Event::G4TrajectoryContainer& trajectories);
void PrintObjectInformation( 
        const Cube::Handle<Cube::ReconObjectContainer>& objects, 
        const Cube::Event::G4TrajectoryContainer& trajectories);
bool isValidObject(const Cube::Handle<Cube::ReconObject>& object);
bool isTPCObject(const Cube::Handle<Cube::ReconObject>& object);
void PrintTrackInformation(
        const Cube::Handle<Cube::ReconTrack>& track,
        const Cube::Event::G4TrajectoryContainer& trajectories);
void PrintClusterInformation(
        const Cube::Handle<Cube::ReconCluster>& cluster,
        const Cube::Event::G4TrajectoryContainer& trajectories);

#endif
