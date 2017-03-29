#path to cmake
CMEXE="/c/Program Files (x86)/CMake/bin/cmake.exe"
#build type: Release, Debug, RelWithDebInfo, MinSizeRel
CMBT="Release"
#make program
CMMAKE="mingw32-make.exe"
#build executable with pyinstaller
CMBUILDEXE=1
#directory with pyinstaller.exe
CMPYINST_PATH="C:/Python27/Scripts/"
#directory with octave
OCTAVE_PATH="C:/Octave/Octave-4.2.0/"
#build installer with nsis
CMBUILDINSTALLER=1
#directory with nsis
CMNSIS_PATH="C:/Program Files (x86)/NSIS/"
#Install prefix
CMINSTALL_PATH="C:/dev/HybMesh"

"$CMEXE" -G "MSYS Makefiles" .. \
	-DCMAKE_BUILD_TYPE=$CMBT \
	-DCMAKE_MAKE_PROGRAM:PATH=$CMMAKE \
	-DUSE_NSIS:BOOL=$CMBUILDINSTALLER \
	-DUSE_PYINSTALLER:BOOL=$CMBUILDEXE \
	-DPYINST_HINT_PATH:PATH="$CMPYINST_PATH" \
	-DNSIS_HINT_PATH:PATH="$CMNSIS_PATH" \
	-DCMAKE_INSTALL_PREFIX:PATH="$CMINSTALL_PATH" \
	-DOCTAVE_HOME:PATH="$OCTAVE_PATH"

