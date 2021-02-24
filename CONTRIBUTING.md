# Contributing to the [MOAB][homepage] repository

All contributions to configuration, source code, and documentation is accepted through the 
pull request process in [Bitbucket][https://bitbucket.org/fathomteam/moab]. However, if you
plan to submit a large contribution to MOAB, we recommend discussing the changes and motivations
through the use of issue pages, or [developer mailing list][moab-dev@mcs.anl.gov]. 

Please note we have a [code of conduct](CODE_OF_CONDUCT.md), and so please follow it in all 
your interactions with the project team and members.

## Pull Request Process

1. The pull request target should be "master". The "develop" branch will be used by the continuous 
   integration system to ensure source stability and verification.
2. Update the README.md with details of changes to the interface, this includes new environment 
   variables, exposed ports, useful file locations and container parameters. If appropriate, document
   all source artifacts with Doxygen-based comments in order to accurate reflect in the documentation
   pages.
3. The continous integration (CI) builds in the develop branch need to pass all tests strictly, and 
   the nightly [Buildbot][buildbot] system should be green before the merge request can be processed. 
4. Increase the version numbers in any examples files and the README.md to the new version that this
   Pull Request would represent. The versioning scheme we use is [SemVer](http://semver.org/).
5. You may merge the Pull Request in once you have the sign-off of two other developers, or if you 
   do not have permission to do that, you may request the second reviewer to merge it for you. 

MOAB is distributed under the LGPL version-3 license (see LICENSE). The act of submitting a pull 
request will be understood as an affirmation of the following:

  Developer's Certificate of Origin 1.1

  By making a contribution to this project, I certify that:

  (a) The contribution was created in whole or in part by me and I
      have the right to submit it under the open source license
      indicated in the file; or

  (b) The contribution is based upon previous work that, to the best
      of my knowledge, is covered under an appropriate open source
      license and I have the right under that license to submit that
      work with modifications, whether created in whole or in part
      by me, under the same open source license (unless I am
      permitted to submit under a different license), as indicated
      in the file; or

  (c) The contribution was provided directly to me by some other
      person who certified (a), (b) or (c) and I have not modified
      it.

  (d) I understand and agree that this project and the contribution
      are public and that a record of the contribution (including all
      personal information I submit with it, including my sign-off) is
      maintained indefinitely and may be redistributed consistent with
      this project or the open source license(s) involved.

[homepage]: https://sigma.mcs.anl.gov/moab-library
[buildbot]: http://gnep.mcs.anl.gov:8010


