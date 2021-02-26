# Contributing to the [MOAB][homepage] repository

All contributions to configuration, source code, and documentation is accepted through the pull request process in [Bitbucket](https://bitbucket.org/fathomteam/moab). If you plan to submit a large contribution to MOAB, we do recommend discussing the changes and motivations through the use of issue pages, or in the [developer mailing list](mailto:moab-dev@mcs.anl.gov).

Please note we have a [code of conduct](CODE_OF_CONDUCT.md), and so please follow it in all your interactions with the project team and members.

## Pull Request (PR) Process

1. When submitting a PR containing your contribution, please ensure that the destination branch is `master`.
2. All PR changes will be merged to the `develop` branch by maintainers in order to verify that tests in the continuous integration (CI) systems ([Buildbot][buildbot], [CircleCI][circleci]) pass cleanly. This process ensures source stability, and minimizes regressions. The contributor are encouraged to additionally request review of code formatting, test coverage and memory leaks specifically on their changeset.
3. If the PR contains changes to the build system (autotools or CMake), which affect code configuration in certain architectures, ensure that the PR description explicitly states it. If appropriate, changes to README.md should be propagated about the updated build process.
4. All public API changes in source artifacts should include Doxygen-based comments to ensure that the [code documentation][moabdocs] is current for `master` branch. If appropriate, request a version number change during PR submission in line with the [SemVer](http://semver.org/) scheme.
5. Please address comments and concerns from developers in a timely manner when reviewers **Request changes** to a PR changeset.
6. You may request the PR be merged to `master` once you have received the **Approved-by** status from two reviewers.

MOAB is distributed under the LGPL version-3 license (see [LICENSE](LICENSE)). The act of submitting a pull request will be understood as an affirmation of the following:

<ins>*Developer's Certificate of Origin: version 1.1*</ins>

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
[circleci]: https://app.circleci.com/pipelines/bitbucket/fathomteam/moab
[moabdocs]: http://ftp.mcs.anl.gov/pub/fathom/moab-docs/index.html

