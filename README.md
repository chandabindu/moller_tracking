# MOLLER-tracking
Reconstruction and analysis code for MOLLER experiments
(Thanks to SBS-offline)

Build prerequisites:

ROOT version 6
cmake version 3.X or higher
Podd version 1.6 and above

**Recommended build procedure**:

Prior to building moller-tracking, it is optional, but highly recommended to set the environment variable MOLLER_REPLAY pointing to the top-level replay directory directory. This enables the convenience features for adding $MOLLER_REPLAY/replay and $MOLLER_REPLAY/scripts to ROOT's macro path, among other things.

Assuming you have a working ROOT and analyzer build and appropriately configured environment for these packages, and assuming that the environment variable ANALYZER points to the top-level installation directory for Podd, the recommended procedure for installing moller-tracking is to create a build directory parallel to the moller-tracking source directory, and an installation directory that is outside of both the build and source directories in a convenient location.

Starting from the directory containing the moller-tracking source directory, do: 

```shell
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=/path/to/desired/installation/directory ..
make -jN install
```
where N is the number of cores available to use for multithread compilation. 

In the top-level install directory, you will find the folders: 

```shell 
bin etc include lib run_replay_here
```

You should add the line 
```
source /path/to/desired/installation/directory/bin/mollerenv.(c)sh 
```
to your login script, where the choice to source mollerenv.sh or mollerenv.csh depends on your preferred shell. This script sets the environment variable MOLLER to the top-level installation directory for moller-tracking and adds $MOLLER/lib to $LD_LIBRARY_PATH (or $DYLD_LIBRARY_PATH on Mac). 

The folder $MOLLER/etc contains a rootlogon.C file that automates the loading of the shared MOLLER library (libsbs.so) when you start the analyzer, and adds $MOLLER/include to ROOT's include path. 

The folder $MOLLER/include contains all the moller-tracking headers and ROOT dictionary files. 

The folder $MOLLER/lib contains the shared library that must be loaded (or linked to any external program) in order to use moller-tracking classes and methods. 

The folder $MOLLER/run_replay_here contains the "magic" .rootrc file that performs several functions. 

If you run ```analyzer``` from $MOLLER/run_replay_here (or any other folder to which $MOLLER/run_replay_here/.rootrc has been copied), then the following steps will happen when you start analyzer/root: 
* **$MOLLER/etc/rootlogon.C** will be executed, automating the loading of libsbs.so and addition of $MOLLER/include to ROOT's include path.
* If $MOLLER_REPLAY was defined prior to the cmake build of moller-tracking, it will add **$MOLLER_REPLAY/replay**, **$MOLLER_REPLAY/scripts**, and **$MOLLER_REPLAY/onlineGUIconfig** to ROOT's macro path. This means that if macro.C exists in any of these folders, then you can do .x macro.C or .L macro.C from any directory containing a copy of $MOLLER/run_replay_here/.rootrc

