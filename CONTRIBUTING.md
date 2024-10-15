# Contribute to EIRENE

Thank you for your interest in contributing to EIRENE. This guide details how to contribute to EIRENE in a way that is easy for everyone.

Looking for something to work on? Look for open issues.

## Golden rule

__One branch, one feature.__

Contributors can speed up the integration of new features or bug fixes by focussing on one topic per branch / merge request. This allows the reviewer to get a rapid overview on the changes which simplifies the process of merging contributions into our development branch.

## Contribution flow

When contributing to EIRENE, your merge request is subject to review by our maintainers.

When maintainers are reading through a merge request they may request further details or improvements from the contributor.

Once the maintainer is satisfied with the contribution it can be merged into the development branch.

## Testing

Especially, new features require simple test cases that can be run automatically with as little input as possible. Therefore, contributors of new features are strongly encurraged to also contribute a new sample case to eirene/EIRENE-sample-cases>.

## Client side git hooks

The EIRENE repository contains git hooks to ensure certain rules are followed. You should instruct git to use these hooks with this command

git config core.hooksPath hooks

Not doing this may result in code being rejected by the Gitlab server. At present only a pre-commit hook is
present, this can be run at any time by running the script hooks/pre-commit but will always run with a new commit.

## Creating a merge request

Once a new feature has been developed and tested including the running of the Gitlab CI, a merge request should be made to add the new or changed code to the EIRENE develop branch. There are several steps that must be carried out to complete this process.

1) The version number after the merge must be decided. If there are no actual code changes this can be the same as the previous version or otherwise can be increased in line with semantic versioning.

2) The version in version.txt should be updated to this value. If git hooks are set up (see contributing.md) then the pre-commit hook will ensure that the version numbers are the same in the following files.

3) This code should then be committed and a note taken of the git hash.

4) At this point an entry must be added to CHANGELOG.md in line with the requirements given in that file. The git hash of the previous commit should be included in the CHANGELOG.md entry and a final commit should be made.

5) The code should then be pushed to the EIRENE repository and a merge request to develop made.