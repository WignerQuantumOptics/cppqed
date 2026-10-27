# CLAUDE.md — C++QED v3

C++QED is a framework for simulating open quantum dynamics: quantum-jump Monte Carlo
(QJMC) trajectories, the Lindblad master equation, and ensembles. It targets
**production runs on HPC/HTC clusters**: you compile one executable per physical system
and run it as thousands of job-array tasks.

- It is deliberately *not* an interactive prototyping tool (that is QuTiP's niche).
- v3 rewrites v2's OO/template-metaprogramming design in a C++20 functional style
  (concepts, ranges, Boost.Hana) while keeping v2's performance core: compile-time
  resolution of the tensor-product structure.
- Development resumed in 2026 after a pause in 2023. The first end-to-end executable ran
  in April 2026.

License: BSL-1.0. Developer Guide: `docs/`, published at https://wignerquantumoptics.github.io/cppqed/.

## Repository relationships

- **The array and operator foundations live in a separate repo, `cppqed-multidiagonal`**
  (`Traits.h`, `MultiArray.h`, `MultiArrayComplex.h`, `MultiDiagonal.h`), included as
  `<multidiagonal/X.h>`. The switch-over (epic E-1) was done in 2026-10.
  - Changes to those foundations go to the **library**. If a task needs one, stop and say
    so. That repo has its own CLAUDE.md.
  - `CPPQEDutils/FrameworkTraits.h` holds only framework extras on top of
    `<multidiagonal/Traits.h>`. Never duplicate symbols that the library already defines
    (`passByValue_v`, `ReferenceMF`, …). Include it directly where its symbols are used.
  - `SliceIterator.h` exists only here, in `CPPQEDutils/`; where it ends up is open (C-04).
- Source dependencies come through **FetchContent pinned to full commit SHAs**
  (`CPPQEDutils/CMakeLists.txt`): cppqed-multidiagonal, and the upstream popl, pcg-cpp and
  Xoshiro-cpp. There are no submodules any more.
  - An installed cppqed-multidiagonal ≥ 0.4 is preferred when found (`FIND_PACKAGE_ARGS`).
  - Co-development against a sibling checkout:
    `-DFETCHCONTENT_SOURCE_DIR_CPPQED-MULTIDIAGONAL=../cppqed-multidiagonal`.
  - Bumping a pin: update the SHA (and its date comment); the jc_qjmc regression test
    must still pass.
- Cluster compute nodes have no outbound network, so builds there must work offline with
  system Boost and Eigen: `FETCHCONTENT_FULLY_DISCONNECTED=ON` plus
  `FETCHCONTENT_SOURCE_DIR_<NAME>` for each dependency, pointing at copies of `_deps/*-src`
  from a configured build. Don't introduce dependencies that break offline builds.

## Layout

```
CPPQEDutils/      physics-agnostic scaffolding: ODE.h (Boost.Odeint, RKCK54), Trajectory.h (run()),
                  StochasticTrajectory.h (Ensemble, AverageTrajectories), Pars.h (popl + JSON), Random.h,
                  SliceIterator.h, FrameworkTraits.h
CPPQEDcore/
  quantumdata/        StateVector, DensityOperator, LazyDensityOperator (LDO)
  quantumoperator/    SparseMatrix (MultiDiagonal comes from the library)
  quantumtrajectory/  QuantumTrajectory.h, Master.h, QuantumJumpMonteCarlo.h — header-only
  quantumdynamics/    Hamiltonian, Liouvillian, ExactPropagator, ExpectationValues,
                      Composite.h (Broadcaster); namespace quantumdynamics
CPPQEDelements/
  subsystems/     Mode, Qbit, Spin
  interactions/   JaynesCummings.h
examples/         jc_qjmc.cc, the reference executable (CMake builds every *.cc here)
tests/            regression/: jc_qjmc against a stored reference (CTest, Python 3)
```

The reference executable is `examples/jc_qjmc.cc`, target `jc_qjmc` (qubit on axis 0,
mode on axis 1). It shows the whole user-facing API in about 40 lines.

## Build

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Debug -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
cmake --build build --target <executable>
```

- `compile_commands.json` is needed by clangd. Development happens in VS Code on WSL2
  Ubuntu.
- Compile as C++23, which the multidiagonal headers require.
- Expect long template-heavy builds. Build single targets rather than `all` while
  iterating.

## Design principles: preserve these in every change

- **Compile time first.** Everything known at compile time is resolved at compile time:
  - RANK
  - retained/dummy axes (`retainedAxes<i...>` NTTPs)
  - slice ranks
  - Hamiltonian composition (Hana sequences)
  - ODE/streaming modes (template parameters of `run()`)

  Don't turn any of these into runtime values. Hilbert-space *dimensions* are runtime
  (`Extents<R>`) for now.
- **Functional, not OO.** Use free functions constrained by concepts, and no base-class
  hierarchies.
  - A SubSystem is not an object. `make(pars)` returns a plain tuple
    `(cutoff, hamiltonian, liouvillian, expectationValues, descr)`.
- **Views vs owners.** Views are passed by value; owners are move-only. Watch for
  dangling spans/iterators after moves and reallocations.
- **The H/i convention.**
  - Hamiltonian functionals *add* `H/i · ψ` to `dψ/dt`.
  - H is generally **non-Hermitian**: in QJMC and Master contexts it is
    `H_eff = H_phys − (i/2) Σ J†J`.
  - Three time-dependence levels: `H(psi,dpsidt)`, `H(t−t0,psi,dpsidt)` (interaction
    picture), and `H(t,psi,dpsidt,t0)`.
- **Exact propagators.** These are per-subsystem diagonal propagators U that are **not
  necessarily unitary**. For Mode, `U = diag(1, e^{−zt}, …)` with
  `z = κ(2n_th+1) − iδ`.
  - `applyPropagator` *replaces* ψ by Uψ, whereas Hamiltonians *add*.
- **Continuous interaction-picture reset.** This is numerically essential:
  - After every ODE step, the state is transformed back to the Schrödinger picture and
    `t0 ← t`.
  - Without the reset, matrix elements grow like `e^{NκT}`, which is `e^{10000}` for
    N=100, κ=10, T=10.
  - With it, growth stays bounded by `e^{NκΔt}`. Never remove or defer the reset.
- **The interaction picture belongs to the trajectory, not to operators.** Operators
  (MultiDiagonal) carry no frequencies.
  - v2's `dynamic_cast`-based picture-aware coupling (`expINKX(particle::Ptr, …)`) is
    deliberately not reproduced.
  - Where phase tabulation is needed, use a separate `(alpha, freqs)` MultiDiagonal pair
    in the physics layer.
- **Master equation reuses state-vector code.**
  - It broadcasts `H_eff/i` over the rows of ρ, then applies `twoTimesRealPartOfSelf`
    (X + X†) to obtain the full Lindblad equation minus the recycling terms.
  - The recycling terms `J ρ J†` come from the superoperator loop.
  - ρ is re-symmetrized and renormalized after each step.
- **LDO abstracts over the state type.** Expectation values take a
  `lazy_density_operator`, so the same functional works for a pure state, a mixed state,
  and an ensemble average.
- **Composition.**
  - `Broadcaster<RANK, ra...>` lifts rank-|ra| descriptors to RANK.
  - Hamiltonians and ExpectationValues are recursive Hana sequences, so nesting
    composites works to any depth.
  - Liouvillians are flat `std::vector<Lindblad<RANK>>`; nest them by broadcast plus
    concatenation.
  - ExactPropagators are not composed: each leaf contributes one, broadcast to the full
    rank.
- **RNG.** Use `pcg64` with independent stream ordinals for ensembles. **Never use
  `std::mt19937` with seed+i.**
- **Parameters.**
  - `Pars` structs register through `add(mod, op, tjc, "Title", _("name","desc",default,binding), …)`.
  - `mod` suffixes namespace the options (`--cutoffm`, `--deltaq`).
  - `noJSON` excludes initial conditions from the system descriptor.
- **run() is the production interface.**
  - Options: `--T`/`--NDt`, `--Dt`/`--dc`, `--o`, `--sdf`, `--initialFileName`,
    `--streamSwitch`, `--precision`.
  - Runs continue automatically from `<o>.state`, and every output file documents itself.
  - Keep state files binary-exact; continuation requires full precision.

## Current state and known problems

**The master equation does not compile**, so there are no current Master results at all.
Finding IDs (B…) refer to `dev/review-2026-10/findings.md`.

- `Master.h` (and `Schrödinger.h`) still use the API removed in 2bc5d3f8 (2024-04-30):
  `QuantumSystemDynamics.h`, `getHa/getLi/getEx`, `hamiltonian_ns::broadcast`. Nothing
  includes them, so the build doesn't notice (B4). Port Master to the `(ha, li, ev[, ex])`
  API that QJMC uses.
- Once it compiles, three known defects each corrupt it independently (a fourth, library
  finding F5, `twoTimesRealPartOfSelf` quadrupling the diagonal, arrived fixed with the
  switch to the library):
  - B3: `conj(m.rho)` after every step
  - B5: Mode superoperators get √rate instead of rate
  - B7: the `SliceIterator` first-slice offset
- **Confirmation:** a driven-qubit Bloch-equation test comparing three results:
  - Master
  - QJMC ensemble
  - the analytic optical-Bloch steady state

Other known problems:
- `UnaryDiagonalPropagator` operator bodies are unimplemented. This is the highest
  priority, because it blocks exact propagation for Qbit (`propagator()` is defined in
  Qbit.cc, but `qbit::make` doesn't use it yet), Mode (where
  `propagator()` needs reinstating), and the future MDoF subsystem.
- QJMC `stepwise` is broken: it has stale `q.ex`/`q.qsd` references.
  - Open question before reinstating it: stepwise evolves under the physical H and
    accounts for decay through rates. It probably needs only the **unitary part** of the
    propagator, `diag(e^{iδnτ})`. Using Mode's decaying U there would double-count
    dissipation. v2's old code holds the intended answer.
- `makeEnsemble` is commented out, so ensemble runs are disabled.
- Autostop is **broken**: `AutostopHandlerGeneric` assumes a numeric TDP, but TDPs are
  labelled Hana tuples.
  - Callers must pass `trajectory::observerNoOp` explicitly. The default argument is dead
    code because template deduction ignores defaults.
- `quantumdynamics::concatenate` fails with non-constexpr arguments. Concatenate Liouvillians
  manually.
- `SliceIterator` (B7): the first slice drops the parent view's offset, and every
  complete loop dereferences one past the end of the offsets. Top-level QJMC is not
  affected; Master, nested Broadcasters, and `partialTrace` of sliced LDOs are.
- Stubs:
  - `compose()` in Hamiltonian.h (interactions across composite boundaries)
  - `broadcastFor`
  - `Placeholder` slicing beyond RANK==2
  - `InitializeEnsembleFromArrayOnlyArchive`
- In-trajectory log in DC mode doesn't work for `dc > 1`.

## Roadmap: migration from v2, not new features

1. Interaction-picture propagators through Mode, Qbit, and MDoF.
2. Motional degree of freedom (v2 `Particle`/`ParticleCavity`).
   - Momentum-space representation, `Spatial` grids.
   - The kinetic propagator `exp(−i ω_rec k_j² τ)` becomes a `UnaryDiagonalPropagator`.
   - `expINKX(nK, dim)` becomes a single-diagonal MultiDiagonal at offset `+nK`, and
     `cos`/`sin` follow from it.
   - The 80-line `mfComposition` table is replaced by `operator|`.
3. Restore autostop (`hana::for_each` over the TDP, or flatten it to
   `valarray<double>`).
4. Re-verify the master equation (Bloch test above).
5. Documentation and tutorials.

Later or under discussion, so **don't start these unasked**:
- renames: acronym casing (C-72), `qbit` → `qubit` (C-76). The structural renames are done.
- compile-time leading subsystem dimensions
- the `OBSERVER` concept unification in Trajectory.h
- HDF5 serialization

## Testing

`ctest --test-dir build` runs everything; CI runs it on every push.
- `tests/unit/` (target `cppqedTests`, Catch2 v3 as in the library, SHA-pinned):
  - dense-matrix checks of the subsystem operators
  - `testHeaders.cc`, compiled twice as two translation units, so a non-inline
    definition in a framework header fails the link (#51). Keep framework headers
    ODR-clean: `inline` or templates only.
- `tests/regression/`: `jc_qjmc` against a stored reference (`compareTrajectory.py`).
  Regenerate the reference only for an intended change of results, and say so in the
  commit message.
- Findings ledger: new findings get a GitHub issue; cite `#N` in the test name and the
  commit message. `[!mayfail]` marks a test that documents an open issue, and goes when
  the fix lands.

The framework-level suite is still thin. Test physics against analytic results:
- Mode: coherent-state evolution, thermal state
- Qbit: Rabi oscillations
- Jaynes–Cummings: vacuum Rabi splitting
- Master: trace and Hermiticity preservation
- QJMC: ensemble convergence to Master

Operators are truncated, so compare against the truncated dense result, never against
textbook commutators.

## Style and conventions

- camelCase. Keep established physics abbreviations (`dcomp`, `noInit`). Use `FWD(x)`.
  Use `slv<"label">(v)` for labelled output values.
- Doxygen `///` comments. Each file starts with a copyright plus BSL-1.0 line.
- Use physics language in docs and comments ("quanta moved", "phase belongs to the
  diagonal"), not group-theory vocabulary. κ may be set to 1 as the natural unit.
- Be honest about what works. Don't advertise broken features (e.g. autostop) in docs or
  examples.
- Physics facts to keep straight:
  - Generators (a, a†, n, Kerr, two-photon drive, JC coupling, e^{inkx}) are
    multidiagonal.
  - The unitaries they generate (displacement, squeezing) are **not**. The developer
    guide (§7.9) and README state this correctly; keep it that way.

## Working rules (Wigner Research Centre institutional policy)

- **git, gh, glab:** Claude runs them with András's existing logins and remotes: fetch, pull,
  push, issues, labels, milestones, Project items, and CI runs. Ask first before:
  - force-push or rebases of pushed history
  - adding new remotes
  - changing repository or organization settings
  - deleting issues or other GitHub objects
- **Other network commands** (downloads, configure steps that download) and **package
  installs:** only on explicit request.
- **Out-of-tree files and credentials:** commands that delete or move files outside this
  repo, or that create, change, or reveal keys or tokens (e.g. `gh auth refresh`): show the
  command, explain it, and let András run it.

  Building, testing, and editing inside the repo are fine.
- **Destructive local operations:** ask for confirmation first. This covers `rm -r`,
  `git reset --hard`, and `git clean`.
- **Simulation and research data:** process locally only. No telemetry, and no external
  APIs or services unless explicitly requested.
- **Secrets:** never hardcode them.
- **System safety:** never disable antivirus, firewalls, updates, or the VPN, and never
  change network configuration.
