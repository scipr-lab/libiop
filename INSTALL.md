# Installation instructions

Libiop depends on the following libraries:

* [Boost](https://www.boost.org/users/download/) - For program options within instrumentation
* [CMake](https://cmake.org/install/)
* [Google Test (GTest)](http://github.com/google/googletest) - For testing and benchmarking
* [libff](https://github.com/scipr-lab/libff) - For implementations of fields with smooth multiplicative subgroups
* [libfqfft](https://github.com/scipr-lab/libfqfft) - For FFTs onto smooth subgroups of multiplicative fields
* [libsodium](https://download.libsodium.org/doc/installation/) - For blake2b and randomness sampling

Google Test, libff, and libfqfft are setup using git submodules,
and their dependencies can all be installed from most package managers directly.
Boost, CMake, and libsodium can similarly be installed directly from most package managers.
We include the commands to install these dependencies for Ubuntu and Fedora below.

### Dependencies for Ubuntu:
* On Ubuntu 18.04 LTS:

```bash
    $ sudo apt-get install build-essential cmake git libgmp3-dev libprocps-dev libboost-all-dev libssl-dev libsodium-dev
```

* On Ubuntu 16.04 LTS:

```bash
    $ sudo apt-get install build-essential cmake git libgmp3-dev libprocps4-dev libboost-all-dev libssl-dev libsodium-dev
```

### Dependencies for Fedora:

* Fedora 29:

```bash
    $ sudo yum install gcc-c++ cmake make git gmp-devel procps-ng-devel boost-devel libsodium-devel
```


### Dependencies Mac OS X

On Mac OS X, install GMP from MacPorts (`port install gmp`).

MacPorts does not write its libraries into standard system folders, so you
might need to explicitly provide the paths to the header files and libraries by
appending `CXXFLAGS=-I/opt/local/include LDFLAGS=-L/opt/local/lib` to the line
above.

### Setting up git submodules and building

Fetch dependencies from their GitHub repos:

```bash
    $ git submodule init && git submodule update
```

Create a `build` folder:

```bash
    $ mkdir build; cd build
```

To create the Makefile:

```bash
    $ cmake ..
```

Then, to compile the library, tests, and profiling harness, run this within the `libiop/build`
directory:

```bash
    $ make
```

Build flags include:
| Flag | Value | Description |
| ---- | ----- | ----------- |
| NDEBUG | false | Enables debug mode. |
| WITH_PROCPS | ON | Enables `libprocps`, which is by default turned off since it is not supported on some systems such as MacOS. |
| -GNinja |  | Builds with `Ninja` instead. |
