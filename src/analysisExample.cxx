#include "analysisExample.hxx"

int main(int argc, char* argv[]) {
    if (!ParseArgs(argc, argv)) return 0;

    SetBranchAddressInputTree();

    Analysis();
}
//----------------------------------------------------------------------------------------
bool ParseArgs(int argc, char* argv[]) {
    bool status = false;
    const char *optstring = "i:n:h";
    char option;

    optind = 1;
    while ( -1 != (option = getopt(argc, argv, optstring)) ) {
        switch (option) {
            case 'i' :
                {
                    std::string fileName = optarg;
                    inputFile = new TFile(fileName.c_str());
                    status = true;
                    break;
                }
            case 'n' : 
                {
                    eventNum = std::stoi(optarg);
                    SPECIFYEVENT = true;
                    break;
                }
            case 'h' : 
                {
                    PrintSyntax();
                    status = false;
                    break;
                }
            default:
                {
                    PrintSyntax();
                    status = false;
                }
        }; // switch (option)
    }; // while ( getopt() ) 
    if (!status) PrintSyntax();
    return status;
}
//----------------------------------------------------------------------------------------
void PrintSyntax() {
    std::cout << "analysisExample" << std::endl;
    std::cout << "  -i ${fileName} (REQUIRED) : name of input CubeRecon root\n";
    std::cout << "  -n ${eventNum} (OPTIONAL) : specify event to look at. \n";
    std::cout << "                              If not specified, loop all events.\n";
    std::cout << "  -h                        : show this message" << std::endl;
}
//----------------------------------------------------------------------------------------
void SetBranchAddressInputTree() {
    std::cout << __func__ << std::endl;
    try {
        inputTree = (TTree*)inputFile->Get("CubeEvents");
    } catch (...) {
        std::runtime_error("can't get tree \"CubeEvents\" from inputFile");
        return;
    }
    inputTree->SetBranchAddress("Event", &event);
}
//----------------------------------------------------------------------------------------
void Analysis() {
    if (SPECIFYEVENT) {
        PrintEventInformation(eventNum);
    } else {
        for (int i = 0; i < inputTree->GetEntries(); ++i) {
            PrintEventInformation(i);
        }
    }
}
//----------------------------------------------------------------------------------------
void PrintEventInformation(int eventNum) {
    inputTree->GetEntry(eventNum);

    std::cout << "evevnt #" << eventNum << std::endl;
    Cube::Event::G4TrajectoryContainer trajectories = event->G4Trajectories;
    Cube::Handle<Cube::ReconObjectContainer> objects = event->GetObjectContainer();
    if (!objects) 
        return;

    PrintPrimaryInformation(trajectories);
    PrintObjectInformation(objects, trajectories);
    std::cout << "-----------------------------------------------------" << std::endl;
}
//----------------------------------------------------------------------------------------
void PrintPrimaryInformation(const Cube::Event::G4TrajectoryContainer& trajectories) {
    std::cout << __func__ << std::endl;

    std::vector<int> primaries = Cube::Tool::AllPrimaries(*event);
    for (auto i : primaries) {
        std::cout << "trajectories[" << i << "]\n";
        std::cout << " |-pdg code: " << trajectories.at(i)->GetPDGCode() << "\n";
        std::cout << " |-energy: " << trajectories.at(i)->GetInitialMomentum().E() << "\n";
        std::cout << " |-position: " << trajectories.at(i)->GetInitialPosition().X() << ", ";
        std::cout <<  trajectories.at(i)->GetInitialPosition().Y() << ", ";
        std::cout <<  trajectories.at(i)->GetInitialPosition().Z() << "\n";
        std::cout << " |-time: " << trajectories.at(i)->GetInitialPosition().T();
        std::cout << std::endl;
    }
}
//----------------------------------------------------------------------------------------
void PrintObjectInformation( 
        const Cube::Handle<Cube::ReconObjectContainer>& objects, 
        const Cube::Event::G4TrajectoryContainer& trajectories) {
    std::cout << __func__ << std::endl;
    int count = 0;
    for (const auto& o : *objects) {
        if (!isValidObject(o))
            continue;
        Cube::Handle<Cube::ReconTrack> track = o;
        Cube::Handle<Cube::ReconCluster> cluster = o;
        if (track) {
            std::cout << "object #" << count << ":";
            PrintTrackInformation(track, trajectories);
            ++count;
        } else if (cluster) {
            std::cout << "object #" << count << ":";
            PrintClusterInformation(cluster, trajectories);
            ++count;
        }
    }
}
//----------------------------------------------------------------------------------------
bool isValidObject(const Cube::Handle<Cube::ReconObject>& object) {
    if (!object)
        return false;
    if (isTPCObject(object))
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
bool isTPCObject(const Cube::Handle<Cube::ReconObject>& object) {
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
void PrintTrackInformation(
        const Cube::Handle<Cube::ReconTrack>& track,
        const Cube::Event::G4TrajectoryContainer& trajectories) {
    //get traj Id corresponding this track
    int mainTraj = Cube::Tool::MainTrajectory(*event, *track);
    std::cout << "track\n";
    std::cout << " |-pdg code: " << trajectories.at(mainTraj)->GetPDGCode() << "\n";
    std::cout << " |-edep: " << track->GetEDeposit() << "\n";
    std::cout << " |-position: " << track->GetPosition().X() << ", ";
    std::cout <<  track->GetPosition().Y() << ", ";
    std::cout <<  track->GetPosition().Z() << "\n";
    std::cout << " |-time: " << track->GetPosition().T();
    std::cout << std::endl;
}
//----------------------------------------------------------------------------------------
void PrintClusterInformation(
        const Cube::Handle<Cube::ReconCluster>& cluster,
        const Cube::Event::G4TrajectoryContainer& trajectories) {
    //get traj Id corresponding this cluster
    int mainTraj = Cube::Tool::MainTrajectory(*event, *cluster);
    std::cout << "cluster\n";
    std::cout << " |-pdg code: " << trajectories.at(mainTraj)->GetPDGCode() << "\n";
    std::cout << " |-edep: " << cluster->GetEDeposit() << "\n";
    std::cout << " |-position: " << cluster->GetPosition().X() << ", ";
    std::cout <<  cluster->GetPosition().Y() << ", ";
    std::cout <<  cluster->GetPosition().Z() << "\n";
    std::cout << " |-time: " << cluster->GetPosition().T();
    std::cout << std::endl;
}
//----------------------------------------------------------------------------------------
