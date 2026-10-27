# Roadmap: GitHub presence: personal repos → WignerQuantumOptics, CI green, issues published

Legend: **[A]** you do it (GitHub UI, `gh`, `git push`, anything network or credential, per Wigner policy) ·
**[C]** I prepare it locally (edits, workflow files, batch scripts) for you to review and run ·
✔ = checkpoint that must hold before moving on.

Starting state (2026-10-08):
- Both repos are under `vukics`, with nothing unpushed. Both have an untracked `CLAUDE.md`; cppqed also has a modified `.gitignore` (my local-settings entry).
- cppqed has no CI, only `jekyll-gh-pages.yml`, which is still tracked even though the last commit message says it was removed. Its executable `jc_qjmc.cc` is gitignored.
- The library has CI (gcc-14/clang-18 × Debug/Release on `ubuntu-latest`) and a docs workflow. (Correction 2026-10-21: it has no Zenodo DOI; Zenodo was never used. The claim came from the library CLAUDE.md.)
- The forks popl, pcg-cpp and Xoshiro-cpp are not needed (verified).

---

## Phase 0: Decisions and a clean local state (no GitHub changes) — done 2026-10-20
1. Decide:
   - ~~the ledger renaming~~ done 2026-10-15: F1, F1b, F2…F11 in both CLAUDE.md files
   - ~~whether to keep Xoshiro~~ decided 2026-10-19: **keep it.** It is a second PRNG with independent streams, so it stays a dependency (C-07) and the RNG tests must cover it alongside `pcg64`.
   - ~~the Pages URL strategy~~ decided 2026-10-19: **`wignerquantumoptics.github.io/<repo>`**, no custom domain.
   - ~~whether the v2-era open issues get closed as obsolete~~ decided 2026-10-19: **close all v2-era open issues and projects as obsolete.** Most issues are already closed. Before closing the projects, skim their items for ideas worth bringing to v3 (Phase 4, step 3).
2. ~~**[C]** Commit the two `CLAUDE.md` files (after the corrections from the review: E2 and the library layout section) and the `.gitignore` change.~~ Done 2026-10-20: E2 applied to cppqed CLAUDE.md, D11 to the library CLAUDE.md. ~~Decide whether `jc_qjmc.cc` gets un-ignored now (C-52/C-74)~~ done 2026-10-20: un-ignored and committed; that's needed for cppqed CI anyway.
3. ~~**[A]** Push both repos.~~ Done 2026-10-20 (the final Phase 0 commits still need one more push in each repo).
4. ~~**[A]** Write down what the transfer could break. For each repo, record:~~ Done 2026-10-20: none of these exist, per András.
   - Pages source and branch
   - Actions secrets and variables (names only)
   - webhooks and installed apps
   - branch protection
   - collaborators
5. ~~**[A]** Dump the existing issues and PRs of both repos (JSON, the command I gave earlier) for deduplication.~~ Done 2026-10-20: neither repo has open issues or PRs, so there is nothing to dedupe against. Also dump the v2-era projects with their items (JSON), so I can skim them for v3 ideas. *Projects dumped 2026-10-20 into `dev/review-2026-10/dumps/` (kept untracked).*

✔ Both repos clean and pushed; inventory of settings and issues saved locally.

## Phase 1: Check and configure the existing WignerQuantumOptics organization (nothing is created) — done 2026-10-22
1. ~~**[A]** Check these organization settings and adjust only where they differ:~~ Done 2026-10-20.
   - require 2FA: already on
   - base member permission: Read (changed from Write; no member lost write access)
   - Actions policy: allow GitHub-owned actions plus a selected list (changed from "all"). Selected: `snok/install-poetry@*` (Wigner_Time), `peaceiris/actions-gh-pages@*` (library docs). Verified by re-running Wigner_Time's tests workflow. SHA pinning not required for now (our workflows use tags; revisit in Phase 5).
   - default `GITHUB_TOKEN`: read-only, already; Actions may not approve PRs
   - Pages allowed for public repos: already
2. ~~**[A]** Define the native **issue types**: Bug, Feature, Task, plus custom Physics, Refactor and Docs if you want them. Types are per organization, so these replace the `type:` label family from the scheme.~~ Done 2026-10-20: Task, Bug, Feature and Decision already existed (set up for Wigner_Time); Physics, Refactor and Docs added. The organization also already has issue fields Priority (Urgent/High/Medium/Low), Effort (High/Medium/Low), Start date and Target date. Mapping in the inventory, §2.
3. ~~**[A]** Create the organization **Project** "C++QED v3"~~ Done 2026-10-21:
   - private since 2026-10-22 (the restricted-visibility issue fields cannot be added to a public Project); one view: Board by Status; further views (by milestone, critical, by repository) can be added in the browser later if needed (views cannot be created via `gh` or the API)
   - no "Roadmap item" field (it would duplicate the milestones)
   - built-in workflows reviewed and kept; auto-add is switched on after the transfer (Phase 2, step 6)
   - Priority and Effort: decided 2026-10-21, the organization issue fields (as Wigner_Time), set on each issue in Phase 4; added to the Project as fields 2026-10-22 (`createProjectV2IssueField`)
4. ~~If you want to double-check, ask IBF whether institute-branded code under the group's organization needs anything.~~ Done: IBF had no comments (András, 2026-10-26).

✔ Org ready; a test repo created in it shows the expected Actions and Pages behaviour (then delete it).

## Phase 2: Transfer — done 2026-10-21
1. ~~**[A]** Transfer `cppqed-multidiagonal` first (it's the dependency), then `cppqed`. In each repo: Settings → Danger Zone → Transfer.~~ Done 2026-10-21; pushing to the new URLs works. The Phase 1 test-repo checkpoint was skipped in favour of transferring directly.
2. **[A]** Never create a repo named `cppqed` or `cppqed-multidiagonal` under `vukics` again; that would break the redirects.
3. ~~**[A]** Point local clones at the new URLs with `git remote set-url origin git@github.com:WignerQuantumOptics/<repo>.git`~~, done for both local checkouts 2026-10-21; cluster checkouts done too (András, 2026-10-26) (they work through the redirect until then).
4. ~~**[A]** Zenodo: enable the library repo under the organization in Zenodo's GitHub settings.~~ Not applicable: Zenodo was never connected (corrected 2026-10-21).
5. ~~**[A]** Archive the three forks under `vukics` (archive, don't delete: old `.gitmodules` reference them).~~ Done 2026-10-21.
6. ~~**[A]** In the "C++QED v3" Project, switch on the built-in "Auto-add to project" workflow for the transferred repos (filter `is:issue is:open`).~~ Done 2026-10-21 for cppqed only: the plan allows one auto-add workflow. Pull requests deliberately not auto-added (direct commits to the default branch; "Fixes #N" closes issues). For the library: issue templates `bug.yml` and `other.yml` with `projects:` pointing at the Project, blank issues disabled (library commit after 6dc2567; the Project is number 1). The Phase 4 script adds its issues itself.

✔ Old URLs redirect (`git ls-remote https://github.com/vukics/cppqed` works); issues, releases and stars are present; Actions still run on push in the library.

## Phase 3: URL and metadata fix-up (one commit per repo) — done 2026-10-21
1. ~~**[C]** Update every hard-coded URL:~~ Done 2026-10-21. `.gitmodules` intentionally still points at the archived forks under `vukics` (until C-07).
   - README badges and links
   - `HOMEPAGE_URL` in both CMakeLists
   - `.gitmodules` (or skip, if C-07 FetchContent lands in Phase 5)
   - Doxyfile and mainpage, the docs workflows
   - both CLAUDE.md files, the developer-guide link
2. ~~**[C]** Add `CITATION.cff` to both repos, so GitHub shows "Cite this repository" (and Zenodo would pick up the metadata, if it is connected later). Optional, but it's the expected norm for research software.~~ Done 2026-10-21: minimal files (author, affiliation, license, URLs; the library also has version 0.2.0 and its tag date). No ORCID or paper reference yet.
3. ~~**[A]** Re-enable Pages under the organization; the sites live at `wignerquantumoptics.github.io/cppqed` and `wignerquantumoptics.github.io/cppqed-multidiagonal` (no custom domain). Check that the developer guide and the library's Doxygen docs build and are reachable.~~ Done 2026-10-21: Pages carried over with the transfer (cppqed: GitHub Actions; library: branch `gh-pages`); both sites load. The library's Doxygen site uses the stock theme → M-20.

✔ No `vukics/` URLs left except intentional historical references; docs sites live at their new URLs.

## Phase 4: Publish the issues — done 2026-10-22
1. **[C]** Turn the inventory into batch scripts that stay idempotent on re-run:
   1. labels: `severity:`, `area:`, `state:`, `origin:` (no `type:` labels: issue types and the Priority/Effort issue fields replace them)
   2. milestones
   3. epics
   4. leaf issues, with bodies that cite paths at a fixed commit SHA
   5. sub-issue links and blocked-by links (including cross-repo, e.g. cppqed C-01 blocked by library M-01)
   6. adding everything to the Project
2. ~~**[A]** Review and run the scripts one repo at a time: library first, so cppqed issues can link to real library numbers.~~ Done 2026-10-22 (run by Claude after dry runs reviewed by András).
   - Library done 2026-10-22: 17 labels, milestones v0.2.1 and v0.3.0, issues #2–#20 created and #1 (M-04) updated, all in Project 1 with type, Priority and Effort; #19 blocked by #14. Scripts kept outside the repo (decided 2026-10-22).
   - cppqed in progress 2026-10-22: 26 labels, milestones M0–M6, epics #14–#20, E-1 children #21–#25; library issues attached to their cppqed epics. #13 (developer-guide typos) closed as fixed in f1f20c85. Master (C-10…C-14 → #26–#30) and QJMC/elements physics (C-20…C-31 → #31–#42) filed. Cluster runs (C-40…C-47 → #43–#50) and build/tests/cleanup (C-50…C-55 → #51–#56) filed; C-56 fixed directly. Roadmap/known problems, renames, mdspan and carry-overs → #57–#86. cppqed complete: 73 issues, all with type, Priority/Effort, Project 1, epics as parents.
3. v2-era issues and projects:
   1. ~~**[C]** Skim the project dumps from Phase 0 for ideas worth carrying over to v3; propose them as inventory additions (`origin:v2-migration`, usually `type:discussion`) for you to accept or drop.~~ Done 2026-10-20: accepted as C-90…C-96 and C-76, plus notes in C-14, C-51, C-53, C-67.
   2. **[A]** ~~Close the remaining v2-era open issues~~ (done 2026-10-20). ~~Close projects #3 "C++QED v3dev" and #4 "C++QED Milestone 13" as obsolete; delete the empty, untitled projects #5 and #6.~~ Done 2026-10-22: #3 and #4 closed (their ideas are in cppqed #80–#86, #73 and #64), #5 and #6 deleted.
4. ~~**[A]** Dump again; **[C]** diff the result against the inventory.~~ Done 2026-10-22: verification queries after each batch (types, fields, Project, milestones, parents, blocked-by) and idempotent re-runs that skipped every filed item.

✔ Every inventory item maps to exactly one GitHub issue; the epics show their children; the Project views (by milestone, by severity, cross-repo blockers) are populated.

## Phase 5: CI working (each change closes the matching issues via "Fixes #N") — done 2026-10-26
**Library** (already has CI; make it strict):
1. ~~**[C]** Fix:~~ Done 2026-10-24 (library #3, #4 closed; the warnings under #10).
   - the Clang warnings (M-10)
   - the const `vectorize`/`matricize` overloads (M-02)
   - the `conj` test gap (M-03): the whole #4 checklist, incl. rank-3 composition and constructJC validation
2. ~~**[C]** Then add:~~ Done 2026-10-24 (#10 closed). The consumer job found a real bug: the installed package exported `cppqed::cppqed-multidiagonal`, not the documented `cppqed::multidiagonal` → library #21, fixed. Actions pinned at checkout v4.4.0 and actions-gh-pages v4.1.0 (latest v4; a major bump is a separate decision).
   - `-DCMAKE_COMPILE_WARNING_AS_ERROR=ON`
   - an ASan/UBSan job
   - a consumer job that installs the library and consumes it via `find_package`, plus one via FetchContent
   - `apt-get update` in `docs.yml`
   - pin `runs-on: ubuntu-24.04` instead of `ubuntu-latest`
   - pin actions by SHA
3. ~~**[C]** Move the headers to `include/multidiagonal/` (M-01).~~ Done 2026-10-24 (#2 closed, f8492fd). ~~Release v0.3.0 **[A tag/release]**~~: released 2026-10-26 by Claude (versioning delegated by András): tag v0.3.0 at fd468e8 plus a GitHub release; milestones rearranged (v0.2.1 → v0.3.1 for the open doc/header items, the open v0.3.0 items → v0.4.0); cppqed pinned to the release commit. Re-versioned 2026-10-27: renaming the misnamed `twoTimesImagPartOf` (library #22) makes the doc/header sweep breaking, so v0.3.1 → v0.4.0 and the previous v0.4.0 → v0.5.0. v0.4.0 released 2026-10-27 (#6, #7, #8, #12, #13, #22; tag at 24d651c), cppqed pinned to it.

**cppqed** (no CI yet):
4. ~~**[C]** Track `jc_qjmc.cc` (C-52/C-74). Add `ci.yml`:~~ Done 2026-10-24 (ddb9df71; green baseline on gcc-13/14, clang-18 × Debug/Release). Also fixed: the version banner still printed the `vukics` URL.
   - ubuntu-24.04, gcc-13/14 and clang-18, × Debug/Release
   - apt packages for Boost and Eigen
   - build all targets
   - a short `jc_qjmc` smoke run
   This first version runs on the *current* code (submodules), which gives a green baseline before anything changes.
5. ~~**[C]** The split itself (E-1, plus the structural renames from E-6 in the same batch):~~ Done 2026-10-24, decided with András: downloads allowed for local verification; all four structural renames (C-70, C-71, C-73, C-74 → #68, #69, #71, #72). Order: regression reference from the pre-switch code (e661220c), renames (9f9ea671), switch (f5d4be20). #21, #22, #25 closed; #14 stays open for #24 (SliceIterator location), #16 for C-08 and library #11. A CI job builds fully disconnected from local source copies (the cluster recipe).
   - FetchContent for the library and the upstream popl/pcg (C-01, C-07)
   - the framework Traits.h (C-02/C-73)
   - the operator ports (C-03)
   - a regression job comparing `jc_qjmc` output against a stored reference (C-05)
   CI must be green before merging.
6. ~~**[C]** Fix the ODR problems (C-50), so a Catch2 test target can be added (C-51). From then on, CI runs `ctest`.~~ Done 2026-10-24 (#51, #23 closed; #52 stays open for the physics tests). New finding: Qbit σz has the opposite sign of spin::sz → #87 (convention decision, `[!mayfail]` test).
7. ~~**[A]** Optional, when CI is stable: a ruleset on `v3/dev` / `main` that requires the CI checks.~~ 2026-10-26: rulesets that block force-pushes and deletion on both default branches (no bypass). Required checks and pull requests come later, tracked in #88.

✔ Both repos green on the default branch; cppqed builds against the library at a pinned SHA; the README CI badges show passing.

## Phase 6: Steady state
- Every fix references its issue, and the regression and physics tests grow with each fix.
- From here the issue order follows the milestones: M1 (correct QJMC: C-20 tolerances, C-21 H_eff, C-23 SliceIterator, C-40/41 run safety) → M2 (Master restored, Bloch test) → M3 (propagators).
- Ledger rule (written into the library CLAUDE.md 2026-10-15): the F-series is closed; new findings = GitHub issue numbers, cited in test names and commits. Added to cppqed's CLAUDE.md (Testing) 2026-10-24, together with its test suite.

---

### Dependencies at a glance
Phase 0 → 1 → 2 → 3, then 4 and 5 can overlap. File the issues first, so the CI and split commits can close them. Within Phase 5: library M-01/M-02 → cppqed split → cppqed tests.
