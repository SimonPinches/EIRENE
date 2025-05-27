# EIRENE

## Legal aspects

The EIRENE code is distributed under the dual "CC BY-NC-ND 4.0" and Forschungszentrum Juelich GmbH (FZJ) developed dedicated "EPL" licence. The user is bound by at least one of the licences at the own choice. The EPL requires the mandatory registration at www.eirene.de website and recognition of the FZ Juelich origin of the EIRENE code, however it is much less restrictive for research-oriented developers. The EPL allows as well a creation of a new developers community so as providing the limited commercial services (for details please check the EPL.md file or www.eirene.de/EPL)

The EIRENE public licence (EPL) is contained in the same directory (EIRENE repository root) in the "EPL.md" file. The licence is versioned together with the code and the manual, which is ensured by Git hooks.

The "CC BY-NC-ND 4.0" licence is also saved to the EIRENE root (as markdown and converted PDF) files based on the orginal licence from https://creativecommons.org/licenses/by-nc-nd/4.0/ converted to markdown format by https://github.com/idleberg/Creative-Commons-Markdown.

The EIRENE code primary application domain is linear kinetic transport, mainly
to study interaction of neutral gas and radiation with magnetized plasmas.

The EIRENE code is a continuously updated "moving target", developed at
Forschungszentrum Jülich GmbH (FZJ) mainly for own scientific applications.

We cannot offer technical support beyond the documentation. On the other hand, EIRENE is an open source code, which we are happy to provide
and jointly employ within research projects and collaborations of mutual interest.



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

A manual can be found on http://www.eirene.de/

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


