ChangeLog:

November 14, 2013: 
- DBH: Added notes on solver details and how to run using generic scripts


***************************** HOW TO ***************************

Notes:
- A simple three block (i.e. CCM) case file for testing reaction kinetics
- system/APOLLOControlDict sets up details for the APOLLO Solver.
 
Warning:
- In this version the THV kinetics are unstable and the option should be left using the butler-volmer formulation for the anode.

Setup:
- Environment:
    - Running scripts have been setup for bash environments
    - please add the following to your .bashrc or .bash_profile
        "source $WM_PROJECT_USER_DIR/FCAPOLLO/etc/bashrc"
    - to compile FCAPOLLO, execute ./buildAPOLLO
    - to extract the case files to the $FOAM_RUN, execute ./extractCaseFiles

    * please note that installing FCAPOLLO requires a functioning version of OpenFoam current versions of OpenFOAM can be found at www.openfoam.com, please refer to the documentation there for setup.

How to run:
- Universal Settings:
    - Transport properties are set in each domain under constant
    - OperatingConditions holds the table for the range of voltages, system/controlDict chooses how many points to do and what point to start on
- Serial Runs:
    - fcCleanCase or fcCleanRun (the later if you want to keep the existing mesh)
    - fcCaseSetup (not required a mesh exists and fcCleanRun was used)
    - fcSerialRun (runs and launches a monitor window)
- Parallel Runs:
    - the number of entries in the hosts.lst determines the number of processors
    - fcCleanCase or fcCleanRun (the later if you want to keep the exisiting mesh)
    - fcCaseSetup (not required if a mesh exists and fcCleanRun was used)
    - fcDecompose (note that fcCleanCase or fcCleanRun require this to be re-run after)
    - fcParallel (generally execute this with "&" to move to a background process)
    - fcMonitor (launch a xterm to monitor the case iterative output)
