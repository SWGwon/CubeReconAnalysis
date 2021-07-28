#ifndef _EVENTSELECTION_
#define _EVENTSELECTION_

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

#include <TRandom.h>
#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TCanvas.h>

#include <iostream>
#include <sstream>
#include <map>
#include <set>
#include <getopt.h>

#include "TG4Event.h"


/// Return a pointer to the EDepSimEvents Tree found in the geometry file.
TTree*    EDepSimTree();

/// Return a pointer to the i'th entry in the tree.  This is an event.
TG4Event* EDepSimEvent(int i);

/// Return a pointer to the last event read.  
TG4Event* EDepSimEvent();

/// Dump the current event in memory.
void EDepSimDumpEvent();
TG4Event* gEDepSimEvent = NULL;


TVector3 beamDirection(0,0,1);
Cube::Event* event = nullptr;
float trueNeutrinoE;
float recoNeutrinoE;
float trueNu;
float trueNuFromCubeRecon;
float recoNu;
float trueMuonE;
float recoMuonE;
float trueNeutronKE;
float recoNeutronKE;
float leverArm;
float tof;
float trackLength;
float angle;
float eDep;
float neighborDistance;
float trueArm;
float trueTof;
float objectPE[10000];
int numberOfBranches;
int primaryPDG;
int parentPDG;
int trackNum;
int objectPDG;
//0: sig track, 1: sig cluster, 2: bkg track, 3: bkg cluster
int category;
float leptonAngle;
float leptonMomentum;
float Q2;
float Q32;
int numberOfFSNeutron;
int numberOfFSPion;
int numberOfFSChargedPion;
int numberOfFSNeutralPion;
int numberOfFSProton;
float maxAngle;
float maxAngleAllIsolated;
float maxAngleCorrespondingDistance;
float maxDistance;
float vtx[4];
int scatteringCount;
bool isNeutronMatched;
int genieStatus[1000]; 
float genieE[1000]; 
int geniePdg[1000]; 
float genieFSMomentumX[1000];
float genieFSMomentumY[1000];
float genieFSMomentumZ[1000];
int eDepProcess[1000]; 
int eventId; 
float eDepProcessT[1000]; 
int genieN; 
TLorentzVector vertex;
TLorentzVector trueVertex;
std::vector<Cube::Handle<Cube::ReconObject>> muonObjectVector;
Cube::Handle<Cube::ReconObject> earliestObject;
Cube::Handle<Cube::ReconTrack> earliestTrack;
Cube::Handle<Cube::ReconCluster> earliestCluster;
float objectX[1000];
float objectY[1000];
float objectZ[1000];
float isolatedObjectX[1000];
float isolatedObjectY[1000];
float isolatedObjectZ[1000];
float object1X[1000];
float object1Y[1000];
float object1Z[1000];
float object2X[1000];
float object2Y[1000];
float object2Z[1000];
float transverseMomentum;
int numberOfObjectFromPrimaryNeutron[1000];

auto outputFile = std::make_shared<TFile> ("variableOutput.root","RECREATE");
auto outputTree = std::make_shared<TTree> ("tree", "tree");

float muonE;

double EvtVtx[4]; 
double StdHepP4[1000][4]; 
int StdHepStatus[1000]; 
int StdHepPdg[1000]; 
int StdHepN; 

double THRESHOLD = 30;
const double VERTEXACTIVITY = 52;
const double MUONACTIVITY = 42.5;
const double BRANCHACTIVITY = 27;

double NEUTRONMATCH = 0;

void Analysis(Cube::Event* event);
int singleTrack = 0;

void SetOutputTree();
double GetTrackLength(Cube::Handle<Cube::ReconTrack>& earliestTrack);
bool isTPC(const Cube::Handle<Cube::ReconObject>& object);
bool isValidObject(const Cube::Handle<Cube::ReconObject>& object);
bool Tsort(Cube::Handle<Cube::ReconObject> o1, Cube::Handle<Cube::ReconObject> o2);
void InitializeOutputVariables();
bool isNuMuBar();
void FillGenieInformation();
std::vector<Cube::Handle<Cube::ReconObject>> GetMuonObjectVector(Cube::Event::G4TrajectoryContainer trajectories, Cube::Handle<Cube::ReconObjectContainer> objects);
TLorentzVector SetTrueVetex(Cube::Event::G4TrajectoryContainer trajectories);
TLorentzVector SetRecoVertex(Cube::Event::G4TrajectoryContainer trajectories);
int GetNumAssociated(Cube::Handle<Cube::ReconObjectContainer> objects);
void FillPrimaryInformation(Cube::Event::G4TrajectoryContainer trajectories);
Cube::Handle<Cube::ReconObject> GetFirstObjectInTime(Cube::Handle<Cube::ReconObjectContainer> objects);
Cube::Handle<Cube::ReconObject> GetFirstObjectInTimeNoIsolatedCondition(Cube::Handle<Cube::ReconObjectContainer> objects, Cube::Event::G4TrajectoryContainer trajectories);
Cube::Handle<Cube::ReconObject> GetFirstObjectInTimeNoMuon(Cube::Handle<Cube::ReconObjectContainer> objects);
void GetTrueArmTime();
void GetMaxAngleDistance(Cube::Handle<Cube::ReconObjectContainer> objects);
int GetNumberOfBranches(Cube::Handle<Cube::ReconObjectContainer> objects);
float GetNeighborDistance(Cube::Handle<Cube::ReconObjectContainer> objects);
void FillOutputVariables(Cube::Handle<Cube::ReconObjectContainer> objects, Cube::Event::G4TrajectoryContainer trajectories);
void FillEDepSimVariables();
bool isTrueCCQE(Cube::Event::G4TrajectoryContainer trajectories);
bool isInFV(TLorentzVector vertex);
bool isInFV1920width(TLorentzVector vertex);
bool isInFV1600width(TLorentzVector vertex);
bool isInFV1920widthObject(TLorentzVector vertex);
bool isInFV1600widthObject(TLorentzVector vertex);
bool isIn1920width;    
bool isIn1600width;
void FillObjectPosition(Cube::Handle<Cube::ReconObjectContainer> objects);
void FillIsolatedObjectPosition(Cube::Handle<Cube::ReconObjectContainer> objects);
void FillNumberOfObjectFromPrimaryNeutron(Cube::Handle<Cube::ReconObjectContainer> objects, Cube::Event::G4TrajectoryContainer trajectories);
double GetTransverseMomentum(const Cube::Handle<Cube::ReconObject>& object);

#endif
