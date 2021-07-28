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
