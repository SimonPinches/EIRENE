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

# get the last two commits which changed the first line of file version.txt
tmp=($(git blame -L1,1 -- version.txt))
last_changed_commit1=${tmp[0]}
tmp=($(git blame -L1,1 -- version.txt "$last_changed_commit1^"))
last_changed_commit2=${tmp[0]}

# if last change to version.txt is not the current version number throw error
test_ver=$(git show $last_changed_commit1:version.txt | head -n 1)
if [ $ver != $test_ver ]; then
    printf "Test failed\nCurrent version number and last chagned version number do not match $ver $test_ver\n"
    exit 1
fi

# Check that the version is not smaller than the previous maj, min or patch versions
oldver=$(git show $last_changed_commit2:version.txt | head -n 1)
echo "Found previous EIRENE version: " $oldver
newversplit=( ${ver//./ } )
oldversplit=( ${oldver//./ } )
if [ ${newversplit[0]} -lt ${oldversplit[0]} ]; then
    printf "Test failed.\nThe new major version (${newversplit[0]}) is lower than the previous major version(${oldversplit[0]}).\n"
    exit 1
elif [ ${newversplit[0]} -eq ${oldversplit[0]} ]; then
    if [ ${newversplit[1]} -lt ${oldversplit[1]} ]; then
	printf "Test failed.\nThe new minor version(${newversplit[1]}) is lower than the previous minor version(${oldversplit[1]}).\n"
	exit 1
    elif [ ${newversplit[1]} -eq ${oldversplit[1]} ]; then
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
# just get latest version number
tmp=($(awk -F'[/=]' '/New_EIRENE_Version/{print $2}' $cha_file))
cha_ver=${tmp[0]}

if [ "$ver" != "$src_ver" ] || [ "$ver" != "$epl_ver" ]  || [ "$ver" != "$man_ver" ] || [ "$ver" != "$cod_ver" ] || [ "$ver" != "$cha_ver" ]; then

    printf "The version numbers in version.txt($ver), eirmod_parmmod.f($man_ver), eirene.tex($man_ver), EPL.md($epl_ver), coding_rules.md($cod_ver) and CHANGELOG.md($cha_ver) do not match.\n"

    cat <<\EOF

Make sure the EIRENE git hooks are enabled in your development environment with

  git config core.hooksPath hooks

EOF
    exit 1
fi
