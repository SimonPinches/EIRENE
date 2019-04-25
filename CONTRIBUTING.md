# Contribute to EIRENE code documentation

Thank you for your interest in contributing to the EIRENE code documentation. This guide details how to contribute in a way that is easy for everyone.

Looking for something to work on? Look for open issues.

## Golden rule

__One branch, one sort of changes.__

Contributors can speed up the integration of corrections or improvements by focussing on one topic per branch / merge request. This allows the reviewer to get a rapid overview on the changes which simplifies the process of merging contributions into our master branch.

## Style guide

* Keep number of Underful/Overful boxes and Warnings as small as possible.
* Use the `cleveref` package for in-document references, e.g. `\cref{.}`, `\Cref{.}`.
* Use the `acronym` package for acronyms, i.e. define new acronyms via `\acro{}{}` at the end of `eirene.tex` and use any acronym via `\ac`, `\acp`, etc.
* Use the `siunitx` package for typesetting units, e.g. `\si{\ampere}`, `\SI{1e2}{\ampere}`.
* Use the `fmtcount` package for typesetting ordinals, e.g. `\ordinalnum{1}`.
* Remember `e.g.` and `i.e.` will trigger full stop space behind, use `e.g.\` and `i.e.\` instead.
* The space between initials and surname should be filled with `~`. It prevents the full stop space and keeps initials and surname in one line.

## Contribution flow

When contributing to the EIRENE code documentation, your merge request is subject to review by our maintainers.

When maintainers are reading through a merge request they may request further details from the contributor.

Once the maintainer is satisfied with the contribution it can be merged into the master branch.