
This code requires VTK, CMake and the gts and glib libraries.

If you have gts and glib installed, then simply create a build directory and use cmake to configure the project:

cd Source
mkdir build
cd build
cmake ../vtkSurfaceBooleanOperations
make



If you cannot get gts or glib from your packages maintainer (e.g. fink on OSX or apt-get Debian Linux) the source is provided in the gts directory. The CMakeLists.txt in the gts directory will try to unpack the tar.gz archives and compile and install them locally.

mkdir gts-build
cd gts-build
cmake ../gts

and set the directory names for the project in Source/vtkSurfaceBooleanOperations





