export CPPFLAGS=""
export CFLAGS="-arch i386 -arch x86_64 -O2"
export CXXFLAGS="-arch i386 -arch x86_64 -O2"
export LDFLAGS=""
export LIBS="-arch i386 -arch x86_64"
export CC="/usr/bin/gcc-4.0"
export CXX="/usr/bin/g++-4.0"
export CPP="/usr/bin/cpp-4.0"
export CXXCPP="/usr/bin/cpp-4.0"
bash configure --disable-dependency-tracking 
