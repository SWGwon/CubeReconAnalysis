# CubeReconAnalysis
This prints basic information of trajectories, objects.
## Entire CubeReconAnalysis folder should be sub-directory of CubeRecon 
Put CubeReconAnalysis into CubeRecon.

At top CMakeLists.txt of CubeRecon, put this line

    add_subdirectory(CubeReconAnalysis)
    
and compile CubeRecon.

## Running
    ./analysisExample 
        -i ${fileName} (REQUIRED) : name of input CubeRecon root
        -n ${eventNum} (OPTIONAL) : specify event to look at.
                                    If not specified, loop all events.
        -h                        : show this message

#### Example

    ./analysisExample -i test.root -n 10

#### Example output
    SetBranchAddressInputTree
    event #10
    PrintTrajectoriesInformation
    trajectories[0]
    |-pdg code: -14
    |-relativistic E: 2534.41MeV
    |-position: 1112.64, -2501.86, 24131.9
    |-time: 1
    |-parentId: -1
    trajectories[1]
    |-pdg code: 2112
    |-relativistic E: 989.128MeV
    |-position: 1112.64, -2501.86, 24131.9
    |-time: 1
    |-parentId: -1
    trajectories[2]
    |-pdg code: 2212
    |-relativistic E: 975.587MeV
    |-position: 1091.69, -2472.14, 24132.5
    |-time: 1.39161
    |-parentId: 1
    PrintObjectInformation
    object #0:cluster
    |-trajectory id: 2
    |-pdg code: 2212
    |-edep: 5080.76pe
    |-position: 1087.47, -2467.23, 24132.5
    |-time: 106.538
    -----------------------------------------------------
