#!/usr/bin/env bash
#
# A script to check that the version number in version.txt is
# the same in the source, manual and license
# Called by the EIRENE CI pipeline with no arguments.  The hook should
# exit with non-zero status after issuing an appropriate message if
# it wants to cause the pipeline to fail.
#

# Read eirene version from version.txt file
ver=$(head -n 1 version.txt)

echo "Running check for consistent EIRENE version" $ver

# Get the version text from  modules/eirmod_parmmod.f90
src_file="src/modules/eirmod_parmmod.f"
src_ver=$(awk 'BEGIN { FS = "=" }/EIRENE_VERSION_STRING/{ print $2 }' $src_file )
src_ver=${src_ver:1:-1}

# Get the version text in EPL.md to the value in version.txt
epl_file="EPL.md"
epl_ver=$(awk '/Version:/{ print $2 }' $epl_file )

# Get the version text in eirene.tex
man_file="Manual/eirene.tex"
man_ver=$(awk '/Manual version/{print $4 }' $man_file )
man_ver=${man_ver::-2}

if [ "$ver" != "$src_ver" ] || [ "$ver" != "$epl_ver" ]  || [ "$ver" != "$man_ver" ]; then

    echo "The version numbers in version.txt($ver), eirmod_parmmod.f($man_ver), eirene.tex($man_ver) and EPL.md($epl_ver) do not match."

    cat <<\EOF

Make sure the EIRENE git hooks are enabled in your development environment with

  git config core.hooksPath hooks

EOF
    exit 1
fi
