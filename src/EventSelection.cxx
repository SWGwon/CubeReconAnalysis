#include "EventSelection.hxx"

std::string cubeReconFile;
std::string genieFile;
std::string genieReWeightFile;

float weight[100];

bool ParseArgs(int argc, char* argv[]);
void PrintSyntax();

bool usingMedian = true;

int main(int argc, char** argv) {
    if (!ParseArgs(argc, argv)) return 0;

    SetOutputTree();

    TFile inputFile(cubeReconFile.c_str());
    //TFile inputEDepSimFile("input_edepsim.root");
    TFile inputGenieFile(genieFile.c_str());
    TFile inputGenieReWeightFile(genieReWeightFile.c_str());
    if (!inputGenieFile.IsOpen() || !inputGenieReWeightFile.IsOpen() || !inputFile.IsOpen())
        return 0;
    if (inputGenieFile.TestBit(TFile::kRecovered) || inputGenieReWeightFile.TestBit(TFile::kRecovered) || inputFile.TestBit(TFile::kRecovered))
        return 0;
    TTree* inputChain = (TTree*)inputFile.Get("CubeEvents");
    inputChain->SetBranchAddress("Event", &event);

    TTree* inputGenieTree = (TTree*)inputGenieFile.Get("gRooTracker");
    inputGenieTree->SetBranchAddress("EvtVtx", &EvtVtx);
    inputGenieTree->SetBranchAddress("StdHepP4", &StdHepP4);
    inputGenieTree->SetBranchAddress("StdHepStatus", &StdHepStatus);
    inputGenieTree->SetBranchAddress("StdHepPdg", &StdHepPdg);
    inputGenieTree->SetBranchAddress("StdHepN", &StdHepN);

    TTree* inputGenieReWeightTree = (TTree*)inputGenieReWeightFile.Get("tree");
    inputGenieReWeightTree->SetBranchAddress("eventnum", &eventnum);
    inputGenieReWeightTree->SetBranchAddress("weight", &weight);

    //TTree* gEDepSimTree = (TTree*)inputEDepSimFile.Get("EDepSimEvents");
    //gEDepSimTree->SetBranchAddress("Event",&gEDepSimEvent);

    for (int i = 0; i < inputChain->GetEntries(); i++) {
        eventId = i;
        inputChain->GetEntry(i);
        inputGenieTree->GetEntry(i);
        inputGenieReWeightTree->GetEntry(i);
        std::cout << "weight[0]: " << weight[0] << std::endl;
        //gEDepSimTree->GetEntry(i);
        Analysis(event);
        std::cout << "event: " << i << std::endl;
    }

    outputFile->Write();
    outputFile->Close();

    return 0;
}
//----------------------------------------------------------------------------------------
bool ParseArgs(int argc, char* argv[]) {
    bool status = false;

    int index;
    int iarg = 0;
    const struct option longopts[] =
    {
        {"input-cuberecon",      required_argument, 0, '1'},
        {"input-genie",          required_argument, 0, '2'},
        {"input-genie-reweight", required_argument, 0, '3'},
        {"threshold",            required_argument, 0, 't'},
        {"help",       no_argument,       0, 'h'},
        {0,0,0,0},
    };

    while (iarg != -1) {
        iarg = getopt_long(argc, argv, "1:2:3:t:", longopts, &index);
        switch (iarg) {
            case '1' :
                {
                    cubeReconFile = optarg;
                    break;
                }
            case '2' :
                {
                    genieFile = optarg;
                    break;
                }
            case '3' :
                {
                    genieReWeightFile = optarg;
                    break;
                }
            case 't' :
                {
                    THRESHOLD = std::stod(optarg);
                    break;
                }
            case 'h' :
                {
                    PrintSyntax();
                    break;
                }
        }
    }

    if (!cubeReconFile.empty() && !genieFile.empty() && !genieReWeightFile.empty()) status = true;
    if (!status) {
        if (cubeReconFile.empty())
            std::cout << "cubeRecon file is required" << std::endl;
        if (genieFile.empty())
            std::cout << "genie file is required" << std::endl;
        if (genieReWeightFile.empty())
            std::cout << "genieReWeight file is required" << std::endl;
        PrintSyntax();
    }
    return status;
}
//----------------------------------------------------------------------------------------
void PrintSyntax() {
    std::cout << "./EventSelection\n";
    std::cout << "  --input-cuberecon ${cuberecon file}            (REQUIRED)\n";
    std::cout << "  --input-genie ${genie file}                    (REQUIRED)\n";
    std::cout << "  --input-genie-reweight ${genie reweight file}  (REQUIRED)\n";
    std::cout << "  -t, --threshold ${threshold}                   (OPTIONAL)\n";
    std::cout << "    : edep threshold (pe) for recon object, default: 30pe \n";
    std::cout << "  -h, --help\n";
    std::cout << "    : show this message\n";
    std::cout << std::endl;
}
//----------------------------------------------------------------------------------------
void Analysis(Cube::Event* inEvent) {
    InitializeOutputVariables();

    if (!isNuMuBar())
        return;

    Cube::Event::G4TrajectoryContainer trajectories = inEvent->G4Trajectories;
    Cube::Handle<Cube::ReconObjectContainer> objects = inEvent->GetObjectContainer();
    if (!objects) 
        return;

    for (int i = 0; i < 100; ++i) {
        reWeight[i] = weight[i];
        if (weight[i] != 0) std::cout << weight[i] << std::endl;
    }
    
    //FillObjectPosition(objects);
    FillPrimaryInformation(trajectories);

    if (numberOfFSNeutron == 0)
        return;

    try {
        muonObjectVector = GetMuonObjectVector(trajectories, objects);
    } catch(...) {
        return;
    }

    try {
        trueVertex = SetTrueVetex();
    } catch(...) {
        return;
    }

    try {
        vertex = SetRecoVertex(trajectories);
    } catch(...) {
        return;
    }

    if (!isInFV(vertex))
        return;

    //if (!isTrueCCQE(trajectories)) {
    //    return;
    //}
    //try {
    //    earliestObject = GetFirstObjectInTimeNoMuon(objects);
    //} catch(...) {
    //    return;
    //}

    //single track
    if (GetNumAssociated(objects) != 1) {
        return;
    }

    try {
        earliestObject = GetFirstObjectInTime(objects);
    } catch(...) {
        return;
    }

    FillGenieInformation();
    FillNumberOfObjectFromPrimaryNeutron(objects, trajectories);
    FillOutputVariables(objects, trajectories);
    FillIsolatedObjectPosition(objects);
    //FillEDepSimVariables();

    outputTree->Fill();
    std::cout << "eventId: " << eventId << ", Fill()" << std::endl;
}
//----------------------------------------------------------------------------------------
void SetOutputTree() {
    outputTree->Branch("trueNeutrinoE", &trueNeutrinoE, "trueNeutrinoE/F");
    outputTree->Branch("recoNeutrinoE", &recoNeutrinoE, "recoNeutrinoE/F");
    outputTree->Branch("trueNu", &trueNu, "trueNu/F");
    outputTree->Branch("trueNuFromCubeRecon", &trueNuFromCubeRecon, "trueNuFromCubeRecon/F");
    outputTree->Branch("recoNu", &recoNu, "recoNu/F");
    outputTree->Branch("trueMuonE", &trueMuonE, "trueMuonE/F");
    outputTree->Branch("recoMuonE", &recoMuonE, "recoMuonE/F");
    outputTree->Branch("trueNeutronKE", &trueNeutronKE, "trueNeutronKE/F");
    outputTree->Branch("recoNeutronKE", &recoNeutronKE, "recoNeutronKE/F");
    outputTree->Branch("leverArm", &leverArm, "leverArm/F");
    outputTree->Branch("tof", &tof, "tof/F");
    outputTree->Branch("trackLength", &trackLength, "trackLength/F");
    outputTree->Branch("angle", &angle, "angle/F");
    outputTree->Branch("eDep", &eDep, "eDep/F");
    outputTree->Branch("neighborDistance", &neighborDistance, "neighborDistance/F");
    outputTree->Branch("trackNum", &trackNum, "trackNum/I");
    outputTree->Branch("category", &category, "category/I");
    outputTree->Branch("eventId", &eventId, "eventId/I");
    outputTree->Branch("primaryPDG", &primaryPDG, "primaryPDG/I");
    outputTree->Branch("parentPDG", &parentPDG, "parentPDG/I");
    outputTree->Branch("objectPDG", &objectPDG, "objectPDG/I");
    outputTree->Branch("scatteringCount", &scatteringCount, "scatteringCount/I");
    outputTree->Branch("numberOfBranches", &numberOfBranches, "numberOfBranches/I");
    outputTree->Branch("leptonAngle", &leptonAngle, "leptonAngle/F");
    outputTree->Branch("leptonMomentum", &leptonMomentum, "leptonMomentum/F");
    outputTree->Branch("trueArm", &trueArm, "trueArm/F");
    outputTree->Branch("trueTof", &trueTof, "trueTof/F");
    outputTree->Branch("Q2", &Q2, "Q2/F");
    outputTree->Branch("Q32", &Q32, "Q32/F");
    outputTree->Branch("objectPE", &objectPE, "objectPE[1000]/F");
    outputTree->Branch("numberOfFSNeutron", &numberOfFSNeutron, "numberOfFSNeutron/I");
    outputTree->Branch("numberOfFSPion", &numberOfFSPion, "numberOfFPion/I");
    outputTree->Branch("numberOfFSChargedPion", &numberOfFSChargedPion, "numberOfFChargedPion/I");
    outputTree->Branch("numberOfFSNeutralPion", &numberOfFSNeutralPion, "numberOfFNeutralPion/I");
    outputTree->Branch("numberOfFSProton", &numberOfFSProton, "numberOfFSProton/I");
    outputTree->Branch("maxAngle", &maxAngle);
    outputTree->Branch("maxAngleAllIsolated", &maxAngleAllIsolated);
    outputTree->Branch("maxAngleCorrespondingDistance", &maxAngleCorrespondingDistance);
    outputTree->Branch("maxDistance", &maxDistance);
    outputTree->Branch("isNeutronMatched", &isNeutronMatched);
    outputTree->Branch("vtx", &vtx, "vtx[4]/F");
    outputTree->Branch("genieStatus", &genieStatus, "genieStatus[1000]/I");
    outputTree->Branch("genieE", &genieE, "genieE[1000]/F");
    outputTree->Branch("geniePdg", &geniePdg, "geniePdg[1000]/I");
    outputTree->Branch("eDepProcess", &eDepProcess, "eDepProcess[1000]/I");
    outputTree->Branch("eDepProcessT", &eDepProcessT, "eDepProcessT[1000]/F");
    outputTree->Branch("objectX", &objectX, "objectX[1000]/F");
    outputTree->Branch("objectY", &objectY, "objectY[1000]/F");
    outputTree->Branch("objectZ", &objectZ, "objectZ[1000]/F");
    outputTree->Branch("isolatedObjectX", &isolatedObjectX, "isolatedObjectX[1000]/F");
    outputTree->Branch("isolatedObjectY", &isolatedObjectY, "isolatedObjectY[1000]/F");
    outputTree->Branch("isolatedObjectZ", &isolatedObjectZ, "isolatedObjectZ[1000]/F");
    outputTree->Branch("object1X", &object1X, "object1X[1000]/F");
    outputTree->Branch("object1Y", &object1Y, "object1Y[1000]/F");
    outputTree->Branch("object1Z", &object1Z, "object1Z[1000]/F");
    outputTree->Branch("object2X", &object2X, "object2X[1000]/F");
    outputTree->Branch("object2Y", &object2Y, "object2Y[1000]/F");
    outputTree->Branch("object2Z", &object2Z, "object2Z[1000]/F");
    outputTree->Branch("genieN", &genieN);
    outputTree->Branch("numberOfObjectFromPrimaryNeutron", &numberOfObjectFromPrimaryNeutron, "numberOfObjectFromPrimaryNeutron[1000]/I");
    outputTree->Branch("reWeight", &reWeight, "reWeight[100]/F");

    outputTree->Branch("isIn1920width", &isIn1920width);
    outputTree->Branch("isIn1600width", &isIn1600width);

}
//----------------------------------------------------------------------------------------
double GetTrackLength(Cube::Handle<Cube::ReconTrack>& inTrack) {
    double tempTrackLength = -1;
    Cube::ReconNodeContainer::iterator n = inTrack->GetNodes().begin();
    Cube::Handle<Cube::TrackState> lastState = (*(n++))->GetState();
    while (n != inTrack->GetNodes().end()) {
        Cube::Handle<Cube::TrackState> nodeState = (*(n++))->GetState();
        tempTrackLength += (nodeState->GetPosition().Vect()
                - lastState->GetPosition().Vect()).Mag();
        lastState = nodeState;
    }
    return tempTrackLength;
}
//----------------------------------------------------------------------------------------
bool isTPC(const Cube::Handle<Cube::ReconObject>& object)
{
    Cube::Handle<Cube::HitSelection> inputHits = object->GetHitSelection();

    // Check if this is for the tpc.
    bool isTPC = false;
    for (Cube::HitSelection::iterator h = inputHits->begin();
         h != inputHits->end(); ++h) {
        if (!Cube::Info::IsTPC((*h)->GetIdentifier())) continue;
        isTPC = true;
        break;
    }
    return isTPC;
}
//----------------------------------------------------------------------------------------
bool isValidObject(const Cube::Handle<Cube::ReconObject>& object) {
    if (!object)
        return false;
    if (isTPC(object))
        return false;

    Cube::Handle<Cube::ReconTrack> track = object;
    Cube::Handle<Cube::ReconCluster> cluster = object;
    if (!track && !cluster)
        return false;
    
    int trajId = Cube::Tool::MainTrajectory(*event, *object);
    if (trajId == -1) 
        return false;

    return true;
}
//----------------------------------------------------------------------------------------
bool Tsort(Cube::Handle<Cube::ReconObject> o1, Cube::Handle<Cube::ReconObject> o2) {
    Cube::Handle<Cube::ReconTrack> track1 = o1;
    Cube::Handle<Cube::ReconCluster> cluster1 = o1;
    Cube::Handle<Cube::ReconTrack> track2 = o2;
    Cube::Handle<Cube::ReconCluster> cluster2 = o2;

    if (track1 && track2) {
        return(track1->GetPosition().T() < track2->GetPosition().T());
    }
    if (track1 && cluster2) {
        return(track1->GetPosition().T() < cluster2->GetMedian().T());
    }
    if (cluster1 && track2) {
        return(cluster1->GetMedian().T() < track2->GetPosition().T());
    }
    if (cluster1 && cluster2) {
        return(cluster1->GetMedian().T() < cluster2->GetMedian().T());
    }

    return false;
}
//----------------------------------------------------------------------------------------
void InitializeOutputVariables() {
    trueNeutrinoE = -10;
    recoNeutrinoE = -10;
    trueNu = -10;
    trueNuFromCubeRecon = -10;
    recoNu = -10;
    trueMuonE = -10;
    recoMuonE = -10;
    trueNeutronKE = -10;
    recoNeutronKE = -10;
    leverArm = -10;
    tof = -10;
    trackLength = -10;
    angle = -10;
    eDep = -10;
    trackNum = -10;
    category = -10;
    primaryPDG = -10;
    parentPDG = -10;
    objectPDG = -10;
    numberOfBranches = -10;
    leptonAngle = -10;
    leptonMomentum = -10;
    Q2 = -10;
    Q32 = -10;
    numberOfFSNeutron = 0;
    numberOfFSPion = 0;
    numberOfFSChargedPion = 0;
    numberOfFSNeutralPion = 0;
    numberOfFSProton = 0;
    maxAngle = 0;
    maxAngleCorrespondingDistance = 0;
    maxDistance = 0;
    scatteringCount = -10;
    isNeutronMatched = false;

    for (int i = 0; i < 1000; ++i) {
        genieStatus[i] = -100;
        geniePdg[i] = -100;
        genieE[i] = -100;
        eDepProcess[i] = -100;
        eDepProcessT[i] = -100;
        objectX[i] = 0;
        objectY[i] = 0;
        objectZ[i] = 0;
        object1X[i] = 0;
        object1Y[i] = 0;
        object1Z[i] = 0;
        object2X[i] = 0;
        object2Y[i] = 0;
        object2Z[i] = 0;
        isolatedObjectX[i] = 0;
        isolatedObjectY[i] = 0;
        isolatedObjectZ[i] = 0;
    }

    muonObjectVector.clear();
    vertex.SetX(0);
    vertex.SetY(0);
    vertex.SetZ(0);
    trueVertex.SetX(0);
    trueVertex.SetY(0);
    trueVertex.SetZ(0);

    isIn1920width = false;
    isIn1600width = false;

    transverseMomentum = -10;
}
//----------------------------------------------------------------------------------------
bool isNuMuBar() {
    if (StdHepPdg[0] == -14)
        return true;
    else
        return false;
}
//---------------------------------------------------------------------------------------- 
void FillGenieInformation() {
    genieN = StdHepN;
    for (int i = 0; i < StdHepN; ++i) {
        geniePdg[i] = StdHepPdg[i];
        genieE[i] = StdHepP4[i][3];
        genieStatus[i] = StdHepStatus[i];
        genieFSMomentumX[i] = StdHepP4[i][0];
        genieFSMomentumY[i] = StdHepP4[i][1];
        genieFSMomentumZ[i] = StdHepP4[i][2];
        if (StdHepStatus[i] == 1 && StdHepPdg[i] == -13) {
            TVector3 muonMomentum;
            muonMomentum.SetX(StdHepP4[i][0]);
            muonMomentum.SetY(StdHepP4[i][1]);
            muonMomentum.SetZ(StdHepP4[i][2]);
            //leptonAngle = muonMomentum.Angle(beamDirection);
            leptonMomentum = muonMomentum.Mag();
            Q2 = 4*StdHepP4[i][3]*StdHepP4[0][3]*std::pow(TMath::Sin(leptonAngle/2),2);
            Q32 = std::pow(std::pow(StdHepP4[0][0],2) + std::pow(StdHepP4[0][1],2) 
                    + std::pow(StdHepP4[0][2],2),0.5) - leptonMomentum;
        }
    }
    vtx[0] = EvtVtx[0];
    vtx[1] = EvtVtx[1];
    vtx[2] = EvtVtx[2];
    vtx[3] = EvtVtx[3];
}
//---------------------------------------------------------------------------------------- 
std::vector<Cube::Handle<Cube::ReconObject>> GetMuonObjectVector(
        Cube::Event::G4TrajectoryContainer trajectories, 
        Cube::Handle<Cube::ReconObjectContainer> objects) {
    std::vector<Cube::Handle<Cube::ReconObject>> tempMuonObjectVector;
    for (const auto& o : *objects) {
        if (!isValidObject(o)) {
            continue;
        }
        Cube::Handle<Cube::ReconTrack> track = o;
        Cube::Handle<Cube::ReconCluster> cluster = o;
        int mainTraj = Cube::Tool::MainTrajectory(*event, *o);
        Cube::Handle<Cube::G4Trajectory> traj = trajectories[mainTraj];
        if (traj && traj->GetPDGCode() == -13 && traj->GetParentId() == -1) {
            tempMuonObjectVector.push_back(o);
        }
    }
    if (tempMuonObjectVector.size() == 0) 
        throw std::runtime_error("tempMuonObjectVector.size() == 0");

    //added 20210711
    if (tempMuonObjectVector.size() != 1) 
        throw std::runtime_error("tempMuonObjectVector.size() != 1, broken muon track");

    return tempMuonObjectVector;
}
//---------------------------------------------------------------------------------------- 
void FillIsolatedObjectPosition(Cube::Handle<Cube::ReconObjectContainer> objects) {
    int count = 0;

    for (const auto& m : muonObjectVector) {
        if (count == 999)
            break;
        if (!isValidObject(m))
            continue;
        for (const auto& o : *objects) {
            if (!isValidObject(o)) {
                continue;
            }
            Cube::Handle<Cube::ReconTrack> track = o;
            Cube::Handle<Cube::ReconCluster> cluster = o;
            //if (track) {
                //if (Cube::Tool::AreNeighboringObjects(*m, *track, MUONACTIVITY) || (track->GetPosition().Vect() - vertex.Vect()).Mag() < VERTEXACTIVITY)
                    //continue;
                //isolatedObjectX[count] = track->GetPosition().X();
                //isolatedObjectY[count] = track->GetPosition().Y();
                //isolatedObjectZ[count] = track->GetPosition().Z();
                //count++;
            //}
            if (cluster) {
                if (Cube::Tool::AreNeighboringObjects(*m, *cluster, MUONACTIVITY) || (cluster->GetPosition().Vect() - vertex.Vect()).Mag() < VERTEXACTIVITY)
                    continue;
                isolatedObjectX[count] = cluster->GetPosition().X();
                isolatedObjectY[count] = cluster->GetPosition().Y();
                isolatedObjectZ[count] = cluster->GetPosition().Z();
                count++;
            }
        }
    }
}
//---------------------------------------------------------------------------------------- 
TLorentzVector SetTrueVetex() {
    Cube::Handle<Cube::ReconObject> recoVertex;
    float tempMuonT = 10E+5;
    for (const auto& m : muonObjectVector) {
        if (!isValidObject(m))
            continue;
        Cube::Handle<Cube::ReconTrack> track = m;
        if (track && track->GetPosition().T() < tempMuonT) {
            tempMuonT = track->GetPosition().T();
            vertex = track->GetPosition();
            recoVertex = m;
        }
    }
    if (tempMuonT == 10E+5) {
        throw std::runtime_error("tempMuonT == 10E+5");
    }

    TLorentzVector tempTrueVertex;

    std::vector<Cube::Handle<Cube::G4Hit>> earliestSegs = Cube::Tool::ObjectG4Hits(*event,*recoVertex);
    double earliestTruth = 1E+8;
    for (std::vector<Cube::Handle<Cube::G4Hit>>::iterator t = earliestSegs.begin(); t != earliestSegs.end(); ++t) { 
        if ((*t)->GetStart().T() < earliestTruth) {
            earliestTruth = (*t)->GetStart().T();
            tempTrueVertex.SetX((*t)->GetStart().X());
            tempTrueVertex.SetY((*t)->GetStart().Y());
            tempTrueVertex.SetZ((*t)->GetStart().Z());
            tempTrueVertex.SetT((*t)->GetStart().T());
        }
    }

    return tempTrueVertex;
}
//---------------------------------------------------------------------------------------- 
//select vertex
//earliest muon 'track' position
TLorentzVector SetRecoVertex(Cube::Event::G4TrajectoryContainer trajectories) {
    TLorentzVector tempVertex;
    float tempMuonT = 10E+5;
    for (const auto& m : muonObjectVector) {
        if (!isValidObject(m))
            continue;
        Cube::Handle<Cube::ReconTrack> track = m;
        if (track && track->GetPosition().T() < tempMuonT) {
            int mainTraj = Cube::Tool::MainTrajectory(*event, *track);
            Cube::Handle<Cube::G4Trajectory> traj = trajectories[mainTraj];
            tempMuonT = track->GetPosition().T();
            muonE = traj->GetInitialMomentum().E();
            tempVertex = track->GetPosition();
        }
    }
    if (tempMuonT == 10E+5) {
        throw std::runtime_error("tempMuonT == 10E+5");
    }

    std::cout << "genie vertex: " << EvtVtx[0]*1000 << ", " << EvtVtx[1]*1000 << ", " << EvtVtx[2]*1000 << std::endl;
    std::cout << "selected vertex: " << vertex.X() << ", " << vertex.Y() << ", " << vertex.Z() << std::endl;

    return tempVertex;
}
//---------------------------------------------------------------------------------------- 
//select channel
//vertex activity : 40m
//single track channel
int GetNumAssociated(Cube::Handle<Cube::ReconObjectContainer> objects) {
    int numAssociated = 0;
    for (const auto& o : *objects) {
        if (!isValidObject(o))
            continue;
        Cube::Handle<Cube::ReconTrack> track = o;
        if (track && (track->GetPosition().Vect() - vertex.Vect()).Mag() < VERTEXACTIVITY) {
            numAssociated++;
        }
    }

    return numAssociated;
}
//---------------------------------------------------------------------------------------- 
void FillPrimaryInformation(Cube::Event::G4TrajectoryContainer trajectories) {
    std::vector<int> primaries = Cube::Tool::AllPrimaries(*event);
    for (auto i : primaries) {
        if (trajectories[i]->GetPDGCode() == 2112) {
            numberOfFSNeutron++;
            trueNuFromCubeRecon += (trajectories[i]->GetInitialMomentum().E() - 939.565);
        }
        if (trajectories[i]->GetPDGCode() == 2212) {
            numberOfFSProton++;
            trueNuFromCubeRecon += (trajectories[i]->GetInitialMomentum().E() - 939.565);
        }
        if (trajectories[i]->GetPDGCode() == 111 || std::abs(trajectories[i]->GetPDGCode()) == 211) {
            numberOfFSPion++;
            trueNuFromCubeRecon += (trajectories[i]->GetInitialMomentum().E());
        }
        if (trajectories[i]->GetPDGCode() == 111) {
            numberOfFSNeutralPion++;
        }
        if (std::abs(trajectories[i]->GetPDGCode()) == 211) {
            numberOfFSChargedPion++;
        }
    }
}
//---------------------------------------------------------------------------------------- 
Cube::Handle<Cube::ReconObject> GetFirstObjectInTime(
        Cube::Handle<Cube::ReconObjectContainer> objects) {
    Cube::Handle<Cube::ReconObject> tempEarliestObject;
    float tempEarliestT = 10E+5;

    for (const auto& m : muonObjectVector) {
        if (!isValidObject(m))
            continue;

        for (const auto& o : *objects) {
            if (!isValidObject(o))
                continue;
            Cube::Handle<Cube::ReconTrack> track = o;
            Cube::Handle<Cube::ReconCluster> cluster = o;

            if (track) {
                if (Cube::Tool::AreNeighboringObjects(*m, *track, MUONACTIVITY) || (track->GetPosition().Vect() - vertex.Vect()).Mag() < VERTEXACTIVITY) {
                    continue;
                } else if (usingMedian && track->GetMedian().T() < tempEarliestT && track->GetEDeposit() > THRESHOLD) {
                    tempEarliestT = track->GetMedian().T();
                    tempEarliestObject = o;
                } else if (!usingMedian && track->GetPosition().T() < tempEarliestT && track->GetEDeposit() > THRESHOLD) {
                    tempEarliestT = track->GetPosition().T();
                    tempEarliestObject = o;
                }
            } else if (cluster) {
                if (Cube::Tool::AreNeighboringObjects(*m, *cluster, MUONACTIVITY) || (cluster->GetPosition().Vect() - vertex.Vect()).Mag() < VERTEXACTIVITY) {
                    continue;
                } else if (usingMedian && cluster->GetMedian().T() < tempEarliestT && cluster->GetEDeposit() > THRESHOLD) {
                    tempEarliestT = cluster->GetMedian().T();
                    tempEarliestObject = o;
                } else if (!usingMedian && cluster->GetPosition().T() < tempEarliestT && cluster->GetEDeposit() > THRESHOLD) {
                    tempEarliestT = cluster->GetPosition().T();
                    tempEarliestObject = o;
                }
            }
        }
    }
    if (!isValidObject(tempEarliestObject) || tempEarliestT == 10E+5)
        throw std::runtime_error("invalid tempEarliestObject");

    return tempEarliestObject;
}
//---------------------------------------------------------------------------------------- 
Cube::Handle<Cube::ReconObject> GetFirstObjectInTimeNoMuon(
        Cube::Handle<Cube::ReconObjectContainer> objects) {
    Cube::Handle<Cube::ReconObject> tempEarliestObject;
    float tempEarliestT = 10E+5;
    for (const auto& o : *objects) {
        if (!isValidObject(o))
            continue;
        Cube::Handle<Cube::ReconTrack> track = o;
        Cube::Handle<Cube::ReconCluster> cluster = o;
        if (track) {
            if ((track->GetPosition().Vect() - vertex.Vect()).Mag() < VERTEXACTIVITY) {
                continue;
            } else if (usingMedian && track->GetMedian().T() < tempEarliestT && track->GetEDeposit() > THRESHOLD) {
                tempEarliestT = track->GetMedian().T();
                tempEarliestObject = o;
            } else if (!usingMedian && track->GetPosition().T() < tempEarliestT && track->GetEDeposit() > THRESHOLD) {
                tempEarliestT = track->GetPosition().T();
                tempEarliestObject = o;
            }
        } else if (cluster) {
            if ((cluster->GetPosition().Vect() - vertex.Vect()).Mag() < VERTEXACTIVITY) {
                continue;
            } else if (usingMedian && cluster->GetMedian().T() < tempEarliestT && cluster->GetEDeposit() > THRESHOLD) {
                tempEarliestT = cluster->GetMedian().T();
                tempEarliestObject = o;
            } else if (!usingMedian && cluster->GetPosition().T() < tempEarliestT && cluster->GetEDeposit() > THRESHOLD) {
                tempEarliestT = cluster->GetPosition().T();
                tempEarliestObject = o;
            }
        }
    }
    if (!isValidObject(tempEarliestObject) || tempEarliestT == 10E+5)
        throw std::runtime_error("invalid tempEarliestObject");

    return tempEarliestObject;
}
//---------------------------------------------------------------------------------------- 
Cube::Handle<Cube::ReconObject> GetFirstObjectInTimeNoIsolatedCondition(
        Cube::Handle<Cube::ReconObjectContainer> objects, 
        Cube::Event::G4TrajectoryContainer trajectories) {
    Cube::Handle<Cube::ReconObject> tempEarliestObject;
    float tempEarliestT = 10E+5;
    for (const auto& o : *objects) {
        if (!isValidObject(o))
            continue;
        Cube::Handle<Cube::ReconTrack> track = o;
        Cube::Handle<Cube::ReconCluster> cluster = o;
        if (track) {
            int mainTraj = Cube::Tool::MainTrajectory(*event, *track);
            Cube::Handle<Cube::G4Trajectory> traj = trajectories[mainTraj];
            if (usingMedian && traj->GetPDGCode() != -13 && track->GetMedian().T() < tempEarliestT && track->GetEDeposit() > THRESHOLD) {
                tempEarliestT = track->GetMedian().T();
                tempEarliestObject = o;
            } else if (!usingMedian && traj->GetPDGCode() != -13 && track->GetPosition().T() < tempEarliestT && track->GetEDeposit() > THRESHOLD) {
                tempEarliestT = track->GetPosition().T();
                tempEarliestObject = o;
            }
        } else if (cluster) {
            int mainTraj = Cube::Tool::MainTrajectory(*event, *cluster);
            Cube::Handle<Cube::G4Trajectory> traj = trajectories[mainTraj];
            if (usingMedian && traj->GetPDGCode() != -13 && cluster->GetMedian().T() < tempEarliestT && cluster->GetEDeposit() > THRESHOLD) {
                tempEarliestT = cluster->GetMedian().T();
                tempEarliestObject = o;
            } else if (!usingMedian && traj->GetPDGCode() != -13 && cluster->GetPosition().T() < tempEarliestT && cluster->GetEDeposit() > THRESHOLD) {
                tempEarliestT = cluster->GetPosition().T();
                tempEarliestObject = o;
            }
        }
    }
    if (!isValidObject(tempEarliestObject) || tempEarliestT == 10E+5)
        throw std::runtime_error("invalid tempEarliestObject");

    return earliestObject;
}
//---------------------------------------------------------------------------------------- 
void GetTrueArmTime() {
    std::vector<Cube::Handle<Cube::G4Hit>> earliestSegs = Cube::Tool::ObjectG4Hits(*event,*earliestObject);
    double earliestTruth = 1E+8;
    for (std::vector<Cube::Handle<Cube::G4Hit>>::iterator t = earliestSegs.begin(); t != earliestSegs.end(); ++t) { 
        if ((*t)->GetStart().T() < earliestTruth) {
            earliestTruth = (*t)->GetStart().T();
            trueTof = (*t)->GetStart().T() - trueVertex.T();
            trueArm = ((*t)->GetStart().Vect() - trueVertex.Vect()).Mag();
        }
    }
    //std::cout << "earliestTruth: " << earliestTruth << std::endl;
    //std::cout << "trueVertex.T(): " << trueVertex.T() << std::endl;
}
//---------------------------------------------------------------------------------------- 
void GetMaxAngleDistance(Cube::Handle<Cube::ReconObjectContainer> objects) {
    std::vector<Cube::Handle<Cube::ReconObject>> isolatedObjectVector;
    int isolatedObjectNum = 0;
    for (const auto& m : muonObjectVector) {
        if (!isValidObject(m))
            continue;
        for (const auto& o : *objects) {
            if (!isValidObject(o))
                continue;
            Cube::Handle<Cube::ReconTrack> track = o;
            Cube::Handle<Cube::ReconCluster> cluster = o;
            if (track) {
                if (Cube::Tool::AreNeighboringObjects(*m, *track, MUONACTIVITY)) {
                    continue;
                } else if (track->GetEDeposit() > THRESHOLD) {
                    isolatedObjectVector.push_back(o);
                    isolatedObjectNum++;
                    objectPE[isolatedObjectNum] = track->GetEDeposit();
                }
            } else if (cluster) {
                if (Cube::Tool::AreNeighboringObjects(*m, *cluster, MUONACTIVITY)) {
                    continue;
                } else if (cluster->GetEDeposit() > THRESHOLD) {
                    isolatedObjectVector.push_back(o);
                    isolatedObjectNum++;
                    objectPE[isolatedObjectNum] = cluster->GetEDeposit();
                }
            }
        }
    }
    //time order
    std::sort(isolatedObjectVector.begin(), isolatedObjectVector.end(), Tsort);

    std::vector<float> isolatedObjectAngle;
    std::vector<float> isolatedObjectDistance;

    for (size_t i = 1; i < isolatedObjectVector.size(); ++i) {
        Cube::Handle<Cube::ReconTrack> track1 = isolatedObjectVector.at(i - 1);
        Cube::Handle<Cube::ReconCluster> cluster1 = isolatedObjectVector.at(i - 1);
        Cube::Handle<Cube::ReconTrack> track2 = isolatedObjectVector.at(i);
        Cube::Handle<Cube::ReconCluster> cluster2 = isolatedObjectVector.at(i);
        TVector3 temp1;
        TVector3 temp2;
        if (track1) {
            temp1 = track1->GetPosition().Vect() - vertex.Vect();
        } else {
            temp1 = cluster1->GetPosition().Vect() - vertex.Vect();
        }
        if (track2) {
            temp2 = track2->GetPosition().Vect() - vertex.Vect();
        } else {
            temp2 = cluster2->GetPosition().Vect() - vertex.Vect();
        }
        isolatedObjectAngle.push_back(temp1.Angle(temp2));
        isolatedObjectDistance.push_back(Cube::Tool::DistanceBetweenObjects(*isolatedObjectVector.at(i - 1), *isolatedObjectVector.at(i)));
    }
    for (size_t ii = 0; ii < isolatedObjectAngle.size(); ++ii) {
        if (isolatedObjectAngle.at(ii) > maxAngle) {
            maxAngle = isolatedObjectAngle.at(ii);
            maxAngleCorrespondingDistance = isolatedObjectDistance.at(ii);
        }
    }
    for (auto i : isolatedObjectDistance) {
        if (i > maxDistance)
            maxDistance = i;
    }

    float tempAngle = 0;
    for (size_t i = 0; i < isolatedObjectVector.size(); ++i) {
        for (size_t j = 0; j < isolatedObjectVector.size(); ++j) {
            Cube::Handle<Cube::ReconTrack> track1 = isolatedObjectVector.at(i);
            Cube::Handle<Cube::ReconCluster> cluster1 = isolatedObjectVector.at(i);
            Cube::Handle<Cube::ReconTrack> track2 = isolatedObjectVector.at(j);
            Cube::Handle<Cube::ReconCluster> cluster2 = isolatedObjectVector.at(j);
            TVector3 temp1;
            TVector3 temp2;
            if (track1) {
                temp1 = track1->GetPosition().Vect() - vertex.Vect();
            } else {
                temp1 = cluster1->GetPosition().Vect() - vertex.Vect();
            }
            if (track2) {
                temp2 = track2->GetPosition().Vect() - vertex.Vect();
            } else {
                temp2 = cluster2->GetPosition().Vect() - vertex.Vect();
            }
            if (tempAngle < temp1.Angle(temp2))
                tempAngle = temp1.Angle(temp2);
        }
    }
    maxAngleAllIsolated = tempAngle;
}
//---------------------------------------------------------------------------------------- 
int GetNumberOfBranches(Cube::Handle<Cube::ReconObjectContainer> objects) {
    int tempNumberOfBranches = 0;
    for (auto& o : *objects) {
        if (!isValidObject(o))
            continue;
        bool isMuon = false;
        for (const auto& m : muonObjectVector) {
            if (o == m) {
                isMuon = true;
                break;
            }
        }
        if (isMuon) continue;
        if (o == earliestObject) continue;
        double distance = 10000;
        distance = Cube::Tool::DistanceBetweenObjects(*earliestObject, *o);
        if (distance < BRANCHACTIVITY)
            tempNumberOfBranches++;
    }

    return tempNumberOfBranches;
}
//---------------------------------------------------------------------------------------- 
float GetNeighborDistance(Cube::Handle<Cube::ReconObjectContainer> objects) {
    neighborDistance = 1E+8;
    for (auto& o : *objects) {
        if (!isValidObject(o))
            continue;
        bool isMuon = false;
        for (const auto& m : muonObjectVector) {
            if (o == m) {
                isMuon = true;
                break;
            }
        }
        if (isMuon) continue;
        if (o == earliestObject) continue;
        double distance = Cube::Tool::DistanceBetweenObjects(*earliestObject, *o);
        if (distance < neighborDistance) {
            neighborDistance = distance;
        }
    }
    if (neighborDistance == 1E+8) {
        neighborDistance = -10;
    }

    return neighborDistance;
}
//---------------------------------------------------------------------------------------- 
void FillOutputVariables(
        Cube::Handle<Cube::ReconObjectContainer> objects, 
        Cube::Event::G4TrajectoryContainer trajectories) {
    GetTrueArmTime();
    GetMaxAngleDistance(objects);
    trackNum  = GetNumAssociated(objects);

    earliestTrack = earliestObject;
    earliestCluster = earliestObject;
    Cube::Handle<Cube::G4Trajectory> earliestTraj = trajectories[Cube::Tool::MainTrajectory(*event, *earliestObject)];
    if (isInFV1920width(vertex)) {
        if (earliestTrack && isInFV1920widthObject(earliestTrack->GetPosition()))
            isIn1920width = true;
        if (earliestCluster && isInFV1920widthObject(earliestCluster->GetPosition()))
            isIn1920width = true;
    }
    if (isInFV1600width(vertex)) {
        if (earliestTrack && isInFV1600widthObject(earliestTrack->GetPosition()))
            isIn1600width = true;
        if (earliestCluster && isInFV1600widthObject(earliestCluster->GetPosition()))
            isIn1600width = true;
    }


    //number of branches
    numberOfBranches = GetNumberOfBranches(objects);

    //neighbor distance
    neighborDistance = GetNeighborDistance(objects);

    //tof;
    if (usingMedian) {
        tof = (earliestTrack ? 
                earliestTrack->GetMedian().T() - vertex.T() : 
                earliestCluster->GetMedian().T() - vertex.T());
    } else if (!usingMedian) {
        tof = (earliestTrack ? 
                earliestTrack->GetPosition().T() - vertex.T() : 
                earliestCluster->GetPosition().T() - vertex.T());
    }
    //earliestTrack ? 
            //std::cout << "earliestTrack->GetMedian().T(): " << earliestTrack->GetMedian().T() << std::endl : 
            //std::cout << "earliestCluster->GetMedian().T(): " << earliestCluster->GetMedian().T() << std::endl;
    //std::cout << "vertex.T(): " << vertex.T() << std::endl;
    //if (tof < 0)
    //return;
    
    //leverArm;
    leverArm = (earliestTrack? 
            (earliestTrack->GetPosition().Vect() - vertex.Vect()).Mag() :
            (earliestCluster->GetPosition().Vect() - vertex.Vect()).Mag());

    //angle;
    angle = (earliestTrack?
            TMath::Cos((earliestTrack->GetPosition().Vect() - vertex.Vect()).Angle(beamDirection)) :
            TMath::Cos((earliestCluster->GetPosition().Vect() - vertex.Vect()).Angle(beamDirection)));
    //eDep;
    eDep = (earliestTrack?
            earliestTrack->GetEDeposit() :
            earliestCluster->GetEDeposit());
    //trueMuonE 
    trueMuonE = muonE;
    //recoMuonE
    recoMuonE = muonE * gRandom->Gaus(1, 0.04);
    //recoNeutronKE
    double recoBeta = (leverArm/tof)/300;
    recoNeutronKE = 939.565*(1./std::pow(1.-std::pow(recoBeta,2),0.5)-1.); 
    //trueNeutrinoE
    trueNeutrinoE = StdHepP4[0][3]*1000.;

    //recoNeutrinoE
    recoNeutrinoE = recoNeutronKE + recoMuonE + 40;
    //trueNu
    for (int k = 0; k < StdHepN; k++)
    {
        if (StdHepPdg[k] == -13)
        {
            trueNu = StdHepP4[0][3]*1000. - StdHepP4[k][3]*1000.;
            break;
        }
    }
    //recoNu
    recoNu = recoNeutronKE;
    //trueNeutronKE

    int earliestTrajID = Cube::Tool::MainTrajectory(*event,*earliestObject);
    int earliestPrim = Cube::Tool::PrimaryId(*event,earliestTrajID);

    //parentPDG
    int parentId = earliestTraj->GetParentId();
    int parentPdg = 0;
    if (parentId == -1) {
        parentPdg = earliestTraj->GetPDGCode();
        //std::cout << "parentPdg: " << parentPdg << std::endl;
        //std::cout << "trajectories[earliestPrim]->GetPDGCode()" << trajectories[earliestPrim]->GetPDGCode() << std::endl;
    } else {
        parentPdg = trajectories[parentId]->GetPDGCode();
    }
    objectPDG = earliestTraj->GetPDGCode();
    parentPDG = parentPdg;
    primaryPDG = trajectories[earliestPrim]->GetPDGCode();

    //0: sig track, 1: sig cluster, 2: bkg track, 3: bkg cluster
    if (earliestTrack) {
        trackLength = GetTrackLength(earliestTrack);
        if (trajectories[earliestPrim]->GetPDGCode() == 2112 && parentPdg == 2112) {
            category = 0;
            trueNeutronKE = trajectories[earliestPrim]->GetInitialMomentum().E() - 939.565;
        }
        else {
            category = 2;
        }
    }

    if (earliestCluster) {
        if (trajectories[earliestPrim]->GetPDGCode() == 2112 && parentPdg == 2112) {
            category = 1;
            trueNeutronKE = trajectories[earliestPrim]->GetInitialMomentum().E() - 939.565;
        } else {
            category = 3;
        }
    }

}
//---------------------------------------------------------------------------------------- 
void FillEDepSimVariables() {
    std::vector<int> process = {};
    std::vector<float> processT = {};
    if (numberOfFSNeutron == 1 && (category == 0 || category == 1)) {
        if (!gEDepSimEvent) {
            std::cout << "Event not available" << std::endl;
        }
        TLorentzVector earliestVect = (earliestTrack?
                earliestTrack->GetPosition() :
                earliestCluster->GetPosition());
        isNeutronMatched = false;
        std::vector<int>::iterator it;
        std::vector<double> process4Time;
        std::vector<double>::iterator processit;

        std::vector<int> parentId;
        std::vector<int> trackId;
        bool isFirst = true;
        int tempParentId = 1000;
        int tempTrackId = 1000;

        std::vector<TG4Trajectory>::iterator matchedTrack;
        std::vector<TG4TrajectoryPoint>::iterator matchedPoint;
        for (std::vector<TG4Trajectory>::iterator
                t = gEDepSimEvent->Trajectories.begin();
                t != gEDepSimEvent->Trajectories.end(); ++t) {
            if (!isNeutronMatched && t->GetPDGCode() == 2112) {
                for (std::vector<TG4TrajectoryPoint>::iterator
                        p = t->Points.begin();
                        p != t->Points.end();
                        ++p) {
                    if ((earliestVect.Vect() - p->GetPosition().Vect()).Mag() < NEUTRONMATCH) {
                        matchedTrack = t;
                        matchedPoint = p;
                        tempParentId = t->GetParentId();
                        tempTrackId = t->GetTrackId();
                        isNeutronMatched = true;
                        break;
                    }
                }
            }
        }
        if (tempParentId != 1000) {
            scatteringCount = 0;
            //std::cout << "-------------------------------" << std::endl;
            //std::cout << "eventId: " << gEDepSimEvent->EventId << std::endl;
            //std::cout << "reco vertex: " << vertex.X() << ", " << vertex.Y() << ", " << vertex.Z() << ", " << vertex.T()<< std::endl;
            //std::cout << "true tof: " << trueTof << ", reco tof: " << tof << std::endl;
            //std::cout << "true arm: " << trueArm << ", reco arm: " << leverArm << std::endl;
            //std::cout << "reco neutron position: " << earliestVect.X() << ", " << earliestVect.Y() << ", " << earliestVect.Z() << ", " << earliestVect.T() << std::endl;
            //std::cout << "matched position: " << matchedPoint->GetPosition().X() << ", " << matchedPoint->GetPosition().Y() << ", " << matchedPoint->GetPosition().Z() << ", " << matchedPoint->GetPosition().T() << std::endl;
            //std::cout << "true neutron position: " << trueEarliestPosition.X() << ", " << trueEarliestPosition.Y() << ", " << trueEarliestPosition.Z() << ", " << trueEarliestPosition.T() << std::endl;

            parentId.push_back(tempParentId);
            trackId.push_back(tempTrackId);
            for (std::vector<TG4Trajectory>::iterator
                    t = gEDepSimEvent->Trajectories.end();
                    t != gEDepSimEvent->Trajectories.begin(); --t) {
                if (tempParentId == t->GetTrackId()) {
                    tempParentId = t->GetParentId();
                    tempTrackId = t->GetTrackId();
                    parentId.push_back(tempParentId);
                    trackId.push_back(tempTrackId);
                }
            }
            std::sort(trackId.begin(), trackId.end());

            for (auto id : trackId) {
                for (std::vector<TG4Trajectory>::iterator
                        t = gEDepSimEvent->Trajectories.begin();
                        t != gEDepSimEvent->Trajectories.end(); ++t) {
                    if (id == t->GetTrackId()) {
                        //std::cout << "   Traj " << t->GetTrackId();
                        //std::cout << " ParentId " << t->GetParentId();
                        //std::cout << " " << t->GetName();
                        //std::cout << " Points.size(): " << t->Points.size();
                        //std::cout << std::endl;
                        std::vector<TG4TrajectoryPoint>::iterator firstPoint;
                        for (std::vector<TG4TrajectoryPoint>::iterator
                                p = t->Points.begin();
                                p != t->Points.end();
                                ++p) {
                            //std::cout << "      ";
                            //std::cout << "position: (" << p->GetPosition().X();
                            //std::cout << ", " << p->GetPosition().Y();
                            //std::cout << ", " << p->GetPosition().Z();
                            //std::cout << ", " << p->GetPosition().T();
                            //std::cout << "), Process: " << p->GetProcess();
                            //std::cout << ", distance from reco neutron: " << (earliestVect.Vect() - p->GetPosition().Vect()).Mag();
                            //std::cout << ", distance from true neutron: " << (trueEarliestPosition.Vect() - p->GetPosition().Vect()).Mag();
                            //std::cout << std::endl;

                            process.push_back(p->GetProcess());
                            processT.push_back(p->GetPosition().T());
                            if (p->GetPosition().X() == matchedPoint->GetPosition().X()) {
                                //std::cout << "      ";
                                //std::cout << "ㄴneutron hit matched, break" << std::endl;
                                break;
                            }
                            if (isFirst && p->GetProcess() == 4) {
                                process4Time.push_back(p->GetPosition().T());
                                scatteringCount++;
                                isFirst = false;
                                firstPoint = p;
                                //std::cout << "      ";
                                //std::cout << "ㄴscatteringCount: " << scatteringCount << std::endl;
                            }
                            processit = std::find(process4Time.begin(), process4Time.end(), p->GetPosition().T());
                            if (processit == process4Time.end() && p->GetProcess() == 4) {
                                process4Time.push_back(p->GetPosition().T());
                                scatteringCount++;
                                //std::cout << "      ";
                                //std::cout << "ㄴscatteringCount: " << scatteringCount << std::endl;
                            }
                            //if (p != firstPoint && processit != process4Time.end()) {
                            //std::cout << "      ";
                            //std::cout << "ㄴalready counted, scatteringCount: " << scatteringCount << std::endl;
                            //}
                        }
                    }
                }
            }
            //std::cout << "total scatteringCount: " << scatteringCount << std::endl;
        } else {
            //std::cout << "-------------------------------" << std::endl;
            //std::cout << "eventId: " << gEDepSimEvent->EventId << std::endl;
            //std::cout << "reco neutron position: " << earliestVect.X() << ", " << earliestVect.Y() << ", " << earliestVect.Z() << ", " << earliestVect.T() << std::endl;
            //std::cout << "not matched" << std::endl;
        }
    }
    for (size_t ii = 0; ii < process.size(); ++ii) {
        eDepProcess[ii] = process.at(ii);
        eDepProcessT[ii] = processT.at(ii);
    }
}
//---------------------------------------------------------------------------------------- 
bool isTrueCCQE(Cube::Event::G4TrajectoryContainer trajectories) {
    int FSNeutron = 0;
    int FSProton = 0;
    int FSPion = 0;
    int FSChargedPion = 0;
    int FSNeutralPion = 0;
    std::vector<int> primaries = Cube::Tool::AllPrimaries(*event);
    for (auto i : primaries) {
        if (trajectories[i]->GetPDGCode() == 2112) {
            FSNeutron++;
        }
        if (trajectories[i]->GetPDGCode() == 2212) {
            FSProton++;
        }
        if (trajectories[i]->GetPDGCode() == 111 || std::abs(trajectories[i]->GetPDGCode()) == 211) {
            FSPion++;
        }
        if (trajectories[i]->GetPDGCode() == 111) {
            FSNeutralPion++;
        }
        if (std::abs(trajectories[i]->GetPDGCode()) == 211) {
            FSChargedPion++;
        }
    }

    if (FSNeutron == 1 && FSProton == 0 && FSPion == 0)
        return true;
    else 
        return false;
}
//---------------------------------------------------------------------------------------- 
bool isInFV(TLorentzVector inVertex) {
    if (std::abs(inVertex.X()) > 2400/2 - 2*7.5 || std::abs(inVertex.Y() + 2600) > 2160/2 - 2*7.5 || std::abs(inVertex.Z() - 23400) > 1920/2 - 2*7.5)
        return false;
    else
        return true;
}
//---------------------------------------------------------------------------------------- 
bool isInFV1920width(TLorentzVector inVertex) {
    if (-1200 + 7.5 > inVertex.X() || inVertex.X() > -1200 + 1920 - 7.5 || std::abs(inVertex.Y() + 2600) > 2160/2 - 2*7.5 || std::abs(inVertex.Z() - 23400) > 1920/2 - 2*7.5)
        return false;
    else
        return true;
}
//---------------------------------------------------------------------------------------- 
bool isInFV1600width(TLorentzVector inVertex) {
    if (-1200 + 7.5 > inVertex.X() || inVertex.X() > -1200 + 1600 - 7.5 || std::abs(inVertex.Y() + 2600) > 2160/2 - 2*7.5|| std::abs(inVertex.Z() - 23400) > 1920/2- 2*7.5)
        return false;
    else
        return true;
}
//---------------------------------------------------------------------------------------- 
bool isInFV1920widthObject(TLorentzVector inVertex) {
    if (-1200 > inVertex.X() || inVertex.X() > -1200 + 1920 || std::abs(inVertex.Y() + 2600) > 2160/2|| std::abs(inVertex.Z() - 23400) > 1920/2)
        return false;
    else
        return true;
}
//---------------------------------------------------------------------------------------- 
bool isInFV1600widthObject(TLorentzVector inVertex) {
    if (-1200 > inVertex.X() || inVertex.X() > -1200 + 1600 || std::abs(inVertex.Y() + 2600) > 2160/2|| std::abs(inVertex.Z() - 23400) > 1920/2)
        return false;
    else
        return true;
}
//---------------------------------------------------------------------------------------- 
void FillNumberOfObjectFromPrimaryNeutron(
        Cube::Handle<Cube::ReconObjectContainer> objects, 
        Cube::Event::G4TrajectoryContainer trajectories) {
    std::vector<int> primaries = Cube::Tool::AllPrimaries(*event);
    std::vector<int> primariyNeutrons;
    for (auto i : primaries) {
        if (trajectories[i]->GetPDGCode() == 2112) {
            primariyNeutrons.push_back(i);
        }
    }

    if (numberOfFSNeutron == 1) {
        int index = 0;
        for (auto o : *objects) {
            if (!isValidObject(o)) {
                continue;
            }
            int earliestTrajID = Cube::Tool::MainTrajectory(*event, *o);
            int earliestPrim = Cube::Tool::PrimaryId(*event, earliestTrajID);
            Cube::Handle<Cube::ReconCluster> cluster = o;
            Cube::Handle<Cube::ReconTrack> track = o;
            //if (primariyNeutrons.at(0) == earliestPrim) {
            //std::cout << "primariyNeutrons.at(0) == earliestPrim" << std::endl;
            if (primariyNeutrons.at(0) == earliestPrim)
                numberOfObjectFromPrimaryNeutron[primariyNeutrons.at(0)]++;
            if (cluster) {
                objectX[index] = cluster->GetPosition().X();
                objectY[index] = cluster->GetPosition().Y();
                objectZ[index] = cluster->GetPosition().Z();
                index++;
            }
            /*
               if (track) {
               objectX[index] = track->GetPosition().X();
               objectY[index] = track->GetPosition().Y();
               objectZ[index] = track->GetPosition().Z();
               index++;
               }
               */
            //}
        }
    }
    if (numberOfFSNeutron == 2) {
        //std::cout << "primariyNeutrons.size(): " << primariyNeutrons.size() << std::endl;
        //std::cout << "numberOfFSNeutron: " << numberOfFSNeutron << std::endl;
        //std::cout << "primariyNeutrons.at(0): " << primariyNeutrons.at(0) << std::endl;
        //std::cout << "primariyNeutrons.at(1): " << primariyNeutrons.at(1) << std::endl;
        int index = 0;
        for (auto o : *objects) {
            if (!isValidObject(o)) {
                continue;
            }
            int earliestTrajID = Cube::Tool::MainTrajectory(*event, *o);
            int earliestPrim = Cube::Tool::PrimaryId(*event, earliestTrajID);
            Cube::Handle<Cube::ReconCluster> cluster = o;
            Cube::Handle<Cube::ReconTrack> track = o;
            //if (primariyNeutrons.at(0) == earliestPrim) {
            //std::cout << "primariyNeutrons.at(0) == earliestPrim" << std::endl;
            if (primariyNeutrons.at(0) == earliestPrim)
                numberOfObjectFromPrimaryNeutron[primariyNeutrons.at(0)]++;
            if (cluster) {
                object1X[index] = cluster->GetPosition().X();
                object1Y[index] = cluster->GetPosition().Y();
                object1Z[index] = cluster->GetPosition().Z();
                //std::cout << "index: " << index << std::endl;
                //std::cout << "cluster->GetPosition().X(): " << cluster->GetPosition().X() << std::endl;;
                //std::cout << "object1X[index]: " << object1X[index] << std::endl;
                index++;
            }
            /*
               if (track) {
               object1X[index] = track->GetPosition().X();
               object1Y[index] = track->GetPosition().Y();
               object1Z[index] = track->GetPosition().Z();
               index++;
               }
               */
            //}
            //if (primariyNeutrons.at(1) == earliestPrim) {
            //numberOfObjectFromPrimaryNeutron[primariyNeutrons.at(1)]++;
            if (primariyNeutrons.at(1) == earliestPrim)
                numberOfObjectFromPrimaryNeutron[primariyNeutrons.at(1)]++;
            //std::cout << "primariyNeutrons.at(1) == earliestPrim" << std::endl;
            if (cluster) {
                object2X[index] = cluster->GetPosition().X();
                object2Y[index] = cluster->GetPosition().Y();
                object2Z[index] = cluster->GetPosition().Z();
                index++;
            }
            /*
               if (track) {
               object2X[index] = track->GetPosition().X();
               object2Y[index] = track->GetPosition().Y();
               object2Z[index] = track->GetPosition().Z();
               index++;
               }
               */
            //}
        }
    }
}
//---------------------------------------------------------------------------------------- 
//double GetTransverseMomentum(const Cube::Handle<Cube::ReconObject>& object) {
//    double tempTransverseMomentum = 0;
//    return tempTransverseMomentum;
//}
