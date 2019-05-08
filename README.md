# EIRENE

## Source code

The source code of EIRENE is [hosted on JuGit](https://jugit.fz-juelich.de/eirene/eirene).

## Monte Carlo transport solver
- multi species
- nonlinear (neutral-neutral)
- time dependent

## Curiosities

Albert Einstein:
"Everything should be made as simple as possible, but not simpler."

Shaw’s Principle:
"Build a system that even a fool can use, and only a fool will want to use it."

Science:
"No agreement between experiment and theory validates a theory (no matter how 
many). But a single discrepancy invalidates a theory"

## Documentation

A manual can be found on http://www.eirene.de/eirene.pdf.

## Atomic and molecular databases

Following atomic and molecular databases are usually used within EIRENE:
- [Amjuel](http://www.eirene.de/amjuel.pdf)
- [Hydhel](http://www.eirene.de/hydhel.pdf)
- [Methane](http://www.eirene.de/methane.pdf)
- [H2vibr](http://www.eirene.de/h2vibr.pdf)

## Installation

EIRENE is served with a CMake config file that allows to control which [interface](src/interfaces) (`EIRENE_INTERFACE`) and [user-routines](src/user-routines) (`EIRENE_USER-ROUTINES`) are compiled with the code. Options to activate trace output (`TRACE`) and usage of the Message Passing Interface (`MPI`) are available.

Typical targets of the generated makefile are `EIRENE` for the EIRENE library to be linked into plasma codes, `eirene` for a standalone version of EIRENE, and `doc` for a Doxygen documentation (by far not complete).

Getting started:
```bash
cd eirene
mkdir buildRelease
cd buildRelease
FC=gfortran cmake ../src
make -j EIRENE
```
Above lines generate an EIRENE library in `eirene/libRelease` with the "Dummy" interface and "default" user-routines, where `j` denotes to use as many threads as available for compilation. If you like to compile for SOLPS-ITER use e.g.
```bash
FC=gfortran cmake ../src -DEIRENE_INTERFACE=SOLPS-ITER -DEIRENE_USER-ROUTINES=iter
```
instead. The variable values correspond to the folder name without prefix (`couple_`, `user_`).

**Attention**

[Interfaces](src/interfaces/) and [user-routines](src/user-routines/) may not be up to date. Use the routines served with the plasma code repository if you are not sure if the EIRENE repository contains the interface and user-routines you want.

## Contributing

EIRENE is an open source code and we are very happy to accept contributions. Please refer to the [contributing guide](CONTRIBUTING.md) for more details.

## Legal aspects

The EIRENE code primary application domain is linear kinetic transport, mainly 
to study interaction of neutral gas and radiation with magnetized plasmas.

The EIRENE code is a continuously updated "moving target", developed at 
Forschungszentrum Jülich GmbH (FZJ) mainly for own scientific applications.

The EIRENE team at FZJ is not a "code provider". We do not offer technical and 
3rd party support at the level which would be necessary for such a complex 
numerical tool.

For this reason no download of EIRENE as ready-to-use software is foreseen.

On the other hand, EIRENE is an open source code, which we are happy to provide 
and jointly employ within projects and collaborations of mutual interest.