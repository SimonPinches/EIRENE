# Rules for the Changelog

- Add new entries at the top of the file
- Include Author, branch and last commit sha-1 hash
- Include the last EIRENE release version the development is based on
- The commit sha-1 hash is obviously not known at the time of the commit, so best practice is to add an additional commit just including an update to the CHANGELOG.md file. This commit can then reference the previous commits sha-1 hash.
- Include list of changes with short description
- If test cases were updated describe the necessity and the amount by which results have changed

---

---

# Changelog of EIRENE repository:

---

---
## D. Harting - DMH/fix_ci_version_check - b6e4b2df6d71f4195faa1e23fa1dd79a0f56a5c8
**Changes**

- Stope the CI pipeline if the compliance_check failed
- Corrected logic for version check in compliance_check
   - compliance_check failed when e.g. minor version number decreased and patch number decreased. This now fixed
   - Check not only the last commit if version number increased but check the last two commits which changed the first line of version.txt if the version number increased

## H. Leggate - HJL/add_changelog_and_coding_rules_version_check - 5cc860bae983754e7aaabb8797f8b340e80e2ac7

#### Based on Vx.y.z(undefined before this commit - use instead develop - 375c8f50e4e5641fed3654c9b59a41e429f590f9)

New_EIRENE_Version=1.0.10

**Changes**

This adds further checks on versions at various points in the code and adds a coding rules document

- Added coding_rules.md imported from ModCR
- Added check on Version in coding_rules.md to pre-commit git hook
- Added check for New_EIRENE_Version in CHANGELOG.md in ci_version_check.sh
- Added short guide to merging in CONTRIBUTING.md
- Changed OpenMP Sample Case version to add in 2 exclusions from 2D-poly

---

## D. Harting - DMH/streamline_ci - 54b400477165ea16b146ca28b0361cf4ff916a48

#### Based on on Vx.y.z (undefined before this commit - use instead develop - 375c8f50e4e5641fed3654c9b59a41e429f590f9)

**Changes**

- All the ITER test cases were building the whole EIRENE library again

  - Build now only once the ITER EIRENE executable at preparation stage to speedup the CI pipeline

  - Reuse the ITER EIRENE executable at run stage

  - No changes in referece outputs of ITER test cases


- Reduced number of MC particles in EMC3 test cases to speedup CI pipeline

  - Reference outputs of EMC3 test cases were updated due to new results with reduced number of MC particles


- Removed -fopenmp compile switch from non openmp test cases

  - Removed -fopenmp switch from CMakeList.txt of EIRENE repository and from the Makefiles in the sample case repository

  - Updated test case ITER_1573_scaling due to small changes from removal of -fopenmp switch

  - Removed ouput files fort.111 and fort.113 generated if CHECKBIN environment variable is set during compile time as these files are now inconsistent with updated test cases and not used during CI checks

---

## X. Bonnin - feature/SOLPS_push_whitespace-comments - 66c03ee2edeeb58def8b6a4ef7b6aefd02ae797b

#### Based on Vx.y.z(undefined before this commit - use instead develop - d3033c4626f220acfb060aecc8c29c523b8f4d4c)

**Changes**

This commit is meant is decrease the distance between the SOLPS Eirene branch and the MsV reference develop branch. It is not meant to contain any functional changes.

- Updated version number to 1.0.9.
- Streamlined detection of graphical library dependencies for compilation.
- Some formatting corrections to avoid stars in output.
- Added a more informative error message when EWALL=0 is set.
- Added informative message from CMake as to what compiler version is being used.
- Added flushing of buffers on error exit to make sure all output is printed.
- Bug fix to OUTAU array size to prevent a CI test failure.
- Corrected generalized interface UPTUSR routine to do the same as its counterparts.
- Added SAVE and THREADPRIVATE statements for local arrays in SUMOSTRA, EIRMOD_COUTAU and EIRMOD_LOCATE routines.
- Alignment bug fix in EIRMOD_COMXS to ensure RP%POLY is allocated.
- Keeping function and routine name declarations on a single line to facilitate code searches.
- Changed CI rules to always produce artifacts.
- Cleaned up debugging output.
- Corrected explanation of SRCML array (folllowing bug fix from Niels Horsten).
- Corrected default value of DTIMV in Manual.
- Added ISPCOPT variable description in Manual.
- Corrected typos.
- Added clarification comments (mostly from Detlev Reiter).
- Right-justifying line labels to satisfy the NAG compiler preprocessor.
- Some aethestic changes in the positioning of empty lines (in both the code source and the output).
- Some minor additional output to match that from the SOLPS branch.
- Removing trailing whitespaces as per the CONTRIBUTING rules.
- Removing tabs as indentation characters.
- Removed some superfluous commented out code lines.
- Removing some superfluous parentheses.
- Truncating long comment lines to avoid compiler warnings.
- Enforced proper indentation when found misaligned.
- Added some spaces to improve code readability.
- Changes of character case to match the SOLPS branch.
- Removed duplicate prb89.dat file.
- Added local user source files to .gitignore list.

---

## D. Harting - DMH/develop_fix_compare_scripts - 0d9077f21ce0c99835d44719821f7e00ad5073f3

#### Based on Vx.y.z (undefined before this commit - use instead develop - d3033c4626f220acfb060aecc8c29c523b8f4d4c)

**Changes**

- Only changes to CI logic, no code changes
  
  - Removed fort.37 from reference outputs of ITER cases. Output to fort.37 was abandoned and is now included in the standard output of EIRENE

- Removed dependency of OPENMP cases on previous cases
  
  - OPENMP library and test case preparation downloaded all previous artifacts. This is not necessary as the OPENMP cases should be fresh compilations

- Removed some cells in the 2D-D_polygon_openmp  test cases from the comparison as they contain small numbers and relative change is large due to statistical noise

---

## D. Harting - DMH/develop_fix_compare_scripts - 7cf909853e2065c8213459946d13e63a9e9a427d

#### Based on Vx.y.z (undefined before this commit - use instead develop - d3033c4626f220acfb060aecc8c29c523b8f4d4c)

**Description**

- CI checks (compare scrips) were switched off and CI signaled 'green' even though the output has changed ([issue #82](https://jugit.fz-juelich.de/eirene/eirene/-/issues/82))
- The checks were switched off by commits 40a86ae38b0eeb860f7f4553a15ade890749f3f3 (2021-10-14) and 28ba4aa745452a67b75d6836c1e60aefa061af20 (2022-12-14)
- Fixed compare scripts to flag differences to reference output
- Updated standard test cases
- **IMPORTANT!!! NEW TEST CASES WERE NOT YET VALIDATED !!!**

**Changes in EIRENE Database since compare error was introduced in  40a86ae38b0eeb860f7f4553a15ade890749f3f3**

- amjuel.tex
  
  - Changes in coefficients:
    - Amjuel reaction 2.2.9  e + H_2  -> 2e + H_2^+
    - Amjuel Reaction 2.3.6A0  C^+ + e -> C
    - Amjuel Reaction 2.6A0  e + C  -> C^+   + 2e (ADAS 93 -> ADAS 96)
    - Amjuel Reaction 2.8A0  e + O  -> O^+   + 2e
  - Changes to T1MIN, T1MAX, N2MIN, N2MAX
    - Amjuel Reaction 2.7A0r  e + N  -> N^+   + 2e
    - Amjuel Reaction 2.6A0 C  + e -> C^+   + 2e
  - Changes to N2MIN, N2MAX
    - Amjuel reaction 2.1.5JH  e + H -> H^+ + 2e
    - Amjuel reaction 2.1.5o  e + H -> H^+ + 2e   Ly-opaque
    - Amjuel Reaction 2.1.8o H^+ + e -> H(1s)  Ly-opaque
    - Amjuel Reaction 2.2.h2c H_2 + e ->  ....
    - Amjuel Reaction 2.2.h2r  H_2 + e ->  ....
    - Amjuel Reaction 2.3.9a  e + He(1s^21S) -> He^+(1s) + 2e
    - Amjuel Reaction 2.3.13a  e+ He^+(1s) ->   He(1s^21S)
  - Changes to T1MIN
    - Amjuel reaction 7.2.3a    p + H^{-} ->  H + H (for cold H^-)
    - Amjuel reaction 7.2.3b   p + H^{-} ->  H + H^+ + 2e (for cold H^-)
  - Changes to P2MIN, P2MAX
    - Amjuel Reaction 2.1.5   H + e -> H^+ +2e,  Ratio H^+/H(1)
    - Amjuel Reaction 2.1.5a  H + e -> H^+ + 2e , Ratio H(3)/H(1)
    - Amjuel Reaction Reaction 2.1.5b  H + e -> H^+ + 2e, Ratio H(2)/H(1)
    - Amjuel eaction 2.1.5c  H + e -> H^+ + 2e, Ratio H(4)/H(1)
    - Amjuel Reaction 2.1.5d  H + e -> H^+ + 2e, Ratio H(5)/H(1)
    - Amjuel Reaction 2.1.5e  H + e -> H^+ + 2e, Ratio H(6)/H(1)
    - Amjuel Reaction 2.1.5tot  H + e -> H^+ + 2e, Ratio H(tot)/H(1)
    - Amjuel Reaction 2.1.5de   H + e -> H^+ + 2e
    - Amjuel Reaction 2.1.5o    H + e -> H^+ + 2e  Ly-opaque
    - Amjuel Reaction 2.1.8   H^+ + e -> H(1s), Ratio $H(1)/H^+
    - Amjuel Reaction 2.1.8a  H^+ + e -> H(1s), Ratio H(3)/H^+
    - Amjuel Reaction 2.1.8b  H^+ + e -> H(1s) ,  Ratio H(2)/H^+
    - Amjuel Reaction 2.1.8c  H^+ + e -> H(1s) ,  Ratio H(4)/H^+
    - Amjuel Reaction 2.1.8d  H^+ + e -> H(1s) ,  Ratio H(5)/H^+
    - Amjuel Reaction 2.1.8e  H^+ + e -> H(1s) ,  Ratio H(6)/H^+
    - Amjuel Reaction 2.1.8tot  H^+ + e -> H(1s) , Ratio H(tot)/H^+
    - Amjuel Reaction 2.1.8de H^+  + e -> H(1s) , + 13.6 eV
    - Amjuel Reaction 2.1.8o H^+ + e -> H(1s) , 13.6 eV
    - Plus others......

- h2vibr.tex
  
  - H.1 :  Fits for sigma(E)
    - Nearly all Eth changed
  - Changes in coefficients 
    - Reaction 2.4l1T
  - A lot of reactions changed T1MIN T1MAX
  - New Reactions added

- hydhel.tex
  
  - Changes of Eth
    - Reaction 3.2.2    p + H_2(v=0) -> p + H_2(v > 0)

- methane.tex
  
  - Changes in coefficients
    - Reaction 3.2      p + C -> C^+ + H 

**Differences of new test cases generated**

- fort.37 generated by eirmod.infcop was abandoned (ecdedfcbb9d2f9a105462ead067205e38b5bc421).
  
  - The output to fort.37 was moved to the general EIRENE output file (IUNOUT)
  - Removed fort.37 from reference cases

- Differences in standard output file of EIRENE
  
  - The maximum number of ATOMIC and molecular reaction has increased by one (NREAC)
    - Affected test cases: 1D-H_slab, 1D_cylindric, 1D-elliptic, 2D-D_slab, 2D-D_slab_finite_cylinder_11x11,  2D-D_triang, 2D_cylindric, 2D_elliptic, 3D_slab_11x11x11,  ALTS2,  ALTS2_bits, cylinder, cylinder_triang_no_octree,  cylinder_triang_octree
    - Test cases not affected:  2D-D_polygon, 2D-D_polygon_NLPLG
  - Some test cases are signalling IEEE_DIVIDE_BY_ZERO
    - 2D-D_polygon, 2D-D_polygon_NLPLG
  - Some test cases are **NOT** signalling any more IEEE_INVALID_FLAG
    - cylinder
  - CREACD changed significantly
    - 2D-D_polygon,  2D-D_slab,  2D-D_slab_finite_cylinder_11x11,  2D-D_triang,  3D_slab_11x11x11,  3D-D_tetra,  3D-D_tetra
  - TABEI1 (AVR) changed significantly 
    - 2D-D_polygon
  - EELEI1 (AVR) changes slightly 
    - 2D-D_polygon_NLPLG 
  - More warnings from SLREAC extrapolation
    - 2D-D_slab_finite_cylinder_11x11, 2D-D_triang, 3D_slab_11x11x11,  3D-D_tetra, 

- **All EMC3 test cases are identical**
  
  - They have limited precision in their ASCII outputs

- Even though these testcases show differences, their global balances in standard output of EIRENE show no differences -> **OK to update these test cases**
  
  - 1D_cylindric, 1D_elliptic, 1D-H_slab, 2D_cylindric, cylinder, ITER_1573, ITER_1573_def_lines, ITER_1573_dens_model, ITER_1573_scaling, ALT2S, 2D-D_polygon_NLPLG, 2D_elliptic, ALT2S_bits
  
  - ITER test cases show differences due to neutral-neutral collisions in the IN-TALLIES

- These test cases show differences with the new EIRENE database, but they show no differences in their global balances with the old EIRENE database -> **OK to update test cases**
  
  - 2D-D_slab_finite_cylinder_11x11, 3D-D_tetra, 2D-D_polygon, 3D-D_slab_11x11x11

- These test cases show differences with the new and the old EIRENE Database. **Reasson for changes needs still to be undestood, not OK to update**
  
  - 2D-D_triang, cylinder_triang_no_octree, cylinder_triang_octree
  
  - 2D-D_slab: Shows only some small differences in global balances with old database !!!

---
