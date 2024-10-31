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

# Check that the version is not smaller than the previous maj, min or patch versions
oldver=$(git show HEAD^:version.txt | head -n 1)
echo "Found in previous commit EIRENE version: " $oldver
newversplit=( ${ver//./ } )
oldversplit=( ${oldver//./ } )
if [ ${newversplit[0]} -lt ${oldversplit[0]} ]; then
    printf "Test failed.\nThe new major version (${newversplit[0]}) is lower than the previous major version(${oldversplit[0]}).\n"
    exit 1
else
    if [ ${newversplit[1]} -lt ${oldversplit[1]} ]; then
	printf "Test failed.\nThe new minor version(${newversplit[1]}) is lower than the previous minor version(${oldversplit[1]}).\n"
	exit 1
    else
	if [ ${newversplit[2]} -lt ${oldversplit[2]} ]; then
	    printf "Test failed.\nThe new patch version(${newversplit[2]}) is lower than the previous patch version(${oldversplit[2]}).\n"
	    exit 1
	fi
    fi
fi

# Get the version text from  modules/eirmod_parmmod.f90
src_file="src/modules/eirmod_parmmod.f"
src_ver=$(awk 'BEGIN { FS = "=" }/EIRENE_VERSION_STRING/{ print $2 }' $src_file )
src_ver=${src_ver:1:-1}

# Get the version text in EPL.md to the value in version.txt
epl_file="EPL.md"
epl_ver=$(awk '/Version:/{print $2}' $epl_file )

# Get the version text in eirene.tex
man_file="Manual/eirene.tex"
man_ver=$(awk '/Manual version/{print $4}' $man_file )
man_ver=${man_ver::-2}

# Get the version text in coding_rules.md
cod_file="coding_rules.md"
cod_ver=$(awk '/Version: /{print $2}' $cod_file )

# Get the version text in change_log.md
cha_file="CHANGELOG.md"
cha_ver=$(awk -F'[/=]' '/New_EIRENE_Version/{print $2}' $cha_file)

if [ "$ver" != "$src_ver" ] || [ "$ver" != "$epl_ver" ]  || [ "$ver" != "$man_ver" ] || [ "$ver" != "$cod_ver" ] || [ "$ver" != "$cha_ver" ]; then

    printf "The version numbers in version.txt($ver), eirmod_parmmod.f($man_ver), eirene.tex($man_ver), EPL.md($epl_ver), coding_rules.md($cod_ver) and CHANGELOG.md($cha_ver) do not match.\n"

    cat <<\EOF

Make sure the EIRENE git hooks are enabled in your development environment with

  git config core.hooksPath hooks

EOF
    exit 1
fi
