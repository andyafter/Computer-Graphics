
Directional Dipole for Subsurface Scattering
====================

# Instructions #

## Requirements
1. cmake
2. CUDA-Toolkit
3. A computer with NVIDIA graphic card 

## Part 1: Setting up your development environment
The project is in the RayTracer folder.

### Windows

1. Make sure you are running Windows 7/8/10 and that your NVIDIA drivers are
   up-to-date. You will need support for OpenGL 4.0 or better in this course.
2. Install Visual Studio 2013 (**not** 2015).
   * 2010/2012 will also work, if you already have one installed.
   * http://www.seas.upenn.edu/cets/software/msdn/
   * You need C++ support. None of the optional components are necessary.
3. Install [CUDA 7](https://developer.nvidia.com/cuda-downloads?sid=925343).
   * CUDA 7.5 is recommended for its new performance profiling tools.
     However, 7.0 is fine (and is the version on the lab computers).
   * Use the Express installation. If using Custom, make sure you select
     Nsight for Visual Studio.
4. Install [CMake](http://www.cmake.org/download/). (Windows binaries are
   under "Binary distributions.")
5. Install [Git](https://git-scm.com/download/win).

### OS X

1. Make sure you are running OS X 10.9 or newer.
2. Install XCode (available for free from the App Store).
   * On 10.10, this may not actually be necessary. Try running `gcc`
     in a terminal first.
3. Install OS X Unix Command Line Development Tools (if necessary).
4. Install [CUDA 7](https://developer.nvidia.com/cuda-downloads?sid=925343)
   (don't use cask; the CUDA cask is outdated).
   * CUDA 7.5 is recommended for its new performance profiling tools.
   * Make sure you select Nsight.
5. Install [Git](https://git-scm.com/download/mac)
   (or: `brew install git`).
6. Install [CMake](http://www.cmake.org/download/)
   (or: `brew cask install cmake`).

### Linux

Note: to debug CUDA on Linux, you will need an NVIDIA GPU with Compute
Capability 5.0.

1. Install [CUDA 7](https://developer.nvidia.com/cuda-downloads?sid=925343).
   * CUDA 7.5 is recommended for its new performance profiling tools.
   * Make sure you select Nsight.
2. Install Git (`apt-get install git` on Debian/Ubuntu).
3. Install CMake (`apt-get install cmake` on Debian/Ubuntu).

## Part 2: Build & Run ##

* `src/` contains the source code.
* `external/` contains the binaries and headers for GLEW and GLFW.

**CMake note:** Do not change any build settings or add any files to your
project directly (in Visual Studio, Nsight, etc.) Instead, edit the
`src/CMakeLists.txt` file. Any files you add must be added here. If you edit it,
just rebuild your VS/Nsight project to make it update itself.

### Windows

1. In Git Bash, navigate to your cloned project directory.
2. Create a `build` directory: `mkdir build`
   * (This "out-of-source" build makes it easy to delete the `build` directory
     and try again if something goes wrong with the configuration.)
3. Navigate into that directory: `cd build`
4. Open the CMake GUI to configure the project:
   * `cmake-gui ..` or `"C:\Program Files (x86)\cmake\bin\cmake-gui.exe" ..`
     * Don't forget the `..` part!
   * Make sure that the "Source" directory is like
     `.../Project0-CUDA-Getting-Started`.
   * Click *Configure*.  Select your version of Visual Studio, Win64.
     (**NOTE:** you must use Win64, as we don't provide libraries for Win32.)
   * If you see an error like `CUDA_SDK_ROOT_DIR-NOTFOUND`,
     set `CUDA_SDK_ROOT_DIR` to your CUDA install path. This will be something
     like: `C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v7.5`
   * Click *Generate*.
5. If generation was successful, there should now be a Visual Studio solution
   (`.sln`) file in the `build` directory that you just created. Open this.
   (from the command line: `explorer *.sln`)
6. Build. (Note that there are Debug and Release configuration options.)
7. Run. Make sure you run the `cis565_` target (not `ALL_BUILD`) by
   right-clicking it and selecting "Set as StartUp Project".
   * If you have switchable graphics (NVIDIA Optimus), you may need to force
     your program to run with only the NVIDIA card. In NVIDIA Control Panel,
     under "Manage 3D Settings," set "Multi-display/Mixed GPU acceleration"
     to "Single display performance mode".

### OS X & Linux

It is recommended that you use Nsight.

1. Open Nsight. Set the workspace to the one *containing* your cloned repo.
2. *File->Import...->General->Existing Projects Into Workspace*.
   * Select the Project 0 repository as the *root directory*.
3. Select the *cis565-* project in the Project Explorer. From the *Project*
   menu, select *Build All*.
   * For later use, note that you can select various Debug and Release build
     configurations under *Project->Build Configurations->Set Active...*.
4. If you see an error like `CUDA_SDK_ROOT_DIR-NOTFOUND`:
   * In a terminal, navigate to the build directory, then run: `cmake-gui ..`
   * Set `CUDA_SDK_ROOT_DIR` to your CUDA install path.
     This will be something like: `/usr/local/cuda`
   * Click *Configure*, then *Generate*.
5. Right click and *Refresh* the project.
6. From the *Run* menu, *Run*. Select "Local C/C++ Application" and the
   `cis565_` binary.
   
## Modify Scenes
The scene files contains the model information. You can modify the objects, but 
remember the following rules:
1. When adding objects the object number must be in the order:
   1,2,3,4,5...
2. You can only specify sphere and cubes
3. You must specify material BSSRDF, and the value 1 means basic BSSRDF, and 2 means directional dipole method


## Additional Information 
If you are not able to use the code instantly after the above tutorial, 
please contact me at andyafter@gmail.com.

