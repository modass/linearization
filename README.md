# Linearization
[![Test](https://github.com/modass/linearization/actions/workflows/test.yml/badge.svg)](https://github.com/modass/linearization/actions/workflows/test.yml)

This is the project repository for the linerarization project using Carleman linearization.

## Dependencies

The project has the following dependencies which need to be available:

* [HyPro](github.com/hypro/hypro), a library for state set representations and linear hybrid systems reachability analysis.

There are additional dependencies, which will be downloaded and built via the build system, e.g., [mc++](https://github.com/coin-or/MCpp).

## Setup

* Clone the repository, e.g., `git clone git@github.com:modass/linearization.git`
* Make sure you have all dependencies available
* Create a build folder (or let your IDE do it for you)
* cd to your build folder and run cmake (`cmake ..`)
* You now should be able to build the targets configured in this project, e.g., `make tool` (for the tool) or `make linearization` (for the library)

The whole process can be run with the following script:
```bash
git clone git@github.com:modass/linearization.git &&\
cd linearization && mkdir build && cd build &&\
cmake .. && make tool
```