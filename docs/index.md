---
layout: default
title: Developer Guide
nav_order: 1
---

# C++QED v3 — Developer Onboarding Guide

*April 2026*

---

## 1. Project Overview

C++QED is a framework for simulating open quantum dynamics, originally developed by András Vukics. Version 3 (v3) is a rewrite of v2, transitioning from a heavily templated object-oriented style to a modern C++20-aided functional style, while preserving the core performance advantage: compile-time resolution of quantum system structure.

The framework is organized into three main libraries:

- **CPPQEDutils** — generic scaffolding (arrays, ODE engine, trajectory driver, stochastic ensemble)
- **CPPQEDcore** — quantum-specific layers (`quantumdata`, `quantumoperator`, `quantumtrajectory`, `quantumdynamics`)
- **CPPQEDelements** — concrete physical systems (SubSystems) and interactions

### 1.1 Design Intent: Production, Not Prototyping

C++QED is explicitly **not** designed for interactive prototyping (Jupyter notebooks, REPL sessions). It occupies a different niche from QuTiP or QuantumOptics.jl: **production-scale simulations** run as compiled executables, submitted as batch jobs on HPC/HTC infrastructure, often as large parameter studies with thousands of independent runs.

|  | QuTiP / QuantumOptics.jl | C++QED v3 |
|---|---|---|
| **Interface** | Interactive notebook / REPL | Compiled executable |
| **Workflow** | Prototype → explore → iterate | Define system → compile → deploy |
| **Execution** | Single workstation, interactive | Command-line, batch jobs, HPC/HTC clusters |
| **Parallelism** | Typically manual | QJMC ensemble trivially parallelizable as independent jobs |
| **Parameter sweep** | Python/Julia loop | Cluster job array, one executable invocation per parameter point |
| **State** | In-memory, session-based | Checkpointed, resumable, serialized |
| **Output** | Inline plots | Structured data files, streamable |

This intent is expressed at every level of the design:
- The **highest-level user interface is an API for writing executables** — `Pars` structs + `run()`, not a scripting interface
- The **command-line parameter machinery** (`Pars.h`) is designed for `--option value` invocation, making job array submission trivial
- The **`run()` function** handles checkpointing, resumption, autostopping (currently broken, see §7.1), and structured logging — all essential for long HPC jobs
- The **compile-time tensor product structure** means you pay the compilation cost once, then run efficiently many times — the right trade-off for production use
- **QJMC ensemble parallelism** in HTC: the natural model is submitting N independent single-trajectory jobs (`nTraj=0`) and averaging results in post-processing, rather than running N trajectories in a single process

C++QED targets the regime where interactive frameworks break down: Hilbert spaces with millions of dimensions, parameter studies requiring thousands of runs, and trajectories requiring days of compute time on cluster nodes.

---

## 2. Core Design Principles

### 2.1 Compile-Time Philosophy

The single most important design principle: everything known at compile time must be processed at compile time. Since one top-level executable simulates one fixed physical system, the entire tensor product structure — number of subsystems, Hilbert space dimensions, axis assignments — is fully known at compile time. This means:

- Retained/dummy axis sets are non-type template parameters (NTTPs), not runtime values
- Slice ranks are compile-time constants
- Hamiltonian composition is a compile-time Hana sequence, not a runtime list
- ODE mode (DC vs DT, T vs NDT) is a template parameter of `run()`

This is not micro-optimization. For Hilbert spaces with millions of dimensions (a practical target), compile-time index arithmetic is what makes the problem tractable. Python/QuTiP frameworks cannot reach this regime efficiently due to runtime reshape/transpose overhead.

> **Note:** Hilbert space *dimensions* are currently runtime values (`Extents<R>`), even for systems whose dimension is structurally fixed (qubits, qutrits). A design for promoting leading dimensions to compile-time NTTPs is under consideration — see section 10.

### 2.2 Functional Style

v3 replaces class hierarchies with free functions constrained by C++20 concepts. There is no Composite base class — composition is done functionally. Operations are algorithms parameterized by `retainedAxes`, acting on plain `MultiArrayView`s.

### 2.3 Value vs Reference Semantics

`MultiArrayView` is a non-owning view (like `std::span`), passed by value. `MultiArray` owns its data and is move-only. `LazyDensityOperator` follows the same pattern. This is reflected in the `passByValue_v` trait.

---

## 3. CPPQEDutils — Scaffolding Layer

### 3.1 MultiArray.h

The fundamental data structure. A rank-R strided array where the mapping from multi-index to flat storage is:

```
position = offset + sum_{k=0}^{R-1} stride_k * index_k
```

| Type | Description |
|---|---|
| `MultiArrayView<T,R>` | Non-owning view over `std::span<T>`. Passed by value. Pointer-like const semantics. |
| `MultiArray<T,R>` | Owns data via `std::vector<T>`. Move-only. Value-like const semantics. Derived from `MultiArrayConstView<T,R>`. |

Key operations: `operator()(multi-index)`, `Placeholder`-based slicing for `RANK==2`, JSON serialization, Boost.Serialization support.

**Known issues:**
- `MultiArray::load()` (Boost.Serialization) does not re-seat the inherited `dataView` span after deserializing `data_` into a new buffer — dangling span bug. Move construction is safe since `std::vector` move preserves the buffer address; the bug is specific to `load()`.
- `assignTo` naming is counterintuitive (copies source into `*this`)
- `zeroInit` uses `T(0.)` — should be `T{}` for generality
- `Placeholder` slicing only implemented for `RANK==2`, needs variadic generalization
- TODO comment identifies a cleaner `MultiArrayBase<StorageType>` design

### 3.2 SliceIterator.h

The core mechanism for acting on subsystems. Given a rank-R array and a compile-time set of retained axes (e.g. `retainedAxes<0,2>` in a rank-4 array):

- `filterIn<retainedAxes>(extents/strides)` — projects onto retained axes
- `filterOut<retainedAxes>(extents/strides)` — projects onto dummy axes
- `calculateSlicesOffsets` — precomputes all flat storage offsets for dummy index combinations (paid once at initialization)
- `sliceRange<retainedAxes>(mav)` — returns a lazy range of `MultiArrayView`s of rank `|retainedAxes|`
- `SliceRangeOwning` / `SliceRangeReferencing` — owns or references the offset vector respectively
- `transpose` — special case where `|retainedAxes| == RANK`

**Note on `SliceRangeOwning`:** it stores a `std::vector<size_t>` of offsets and initializes `b_` (a `SliceIterator`) with `offsets_.begin()` in the member initializer list. This is safe as long as the class stays move-only (as it is now): a defaulted move transfers the vector's buffer, so `b_` keeps pointing at valid storage — the same reason `MultiArray` move construction is safe. Making the class copyable with a defaulted copy would leave `b_` pointing into the source's vector.

**Open items:**
- `broadcastFor` — stub, intended as the high-level functional primitive for applying a function across matching slices of multiple arrays simultaneously
- `MultiIndexRange` — attempt to make `incrementMultiIndex` a proper `std::ranges::range`; currently satisfies `std::input_iterator` but not `std::ranges::range`

> Any future extension carrying compile-time extent information through `sliceRange` must respect the transposition constraint described in section 10.4 — the static prefix of the result depends on the *order* of `retainedAxes`, not just their set.

### 3.3 ODE.h

Adaptive ODE engine wrapping Boost.Odeint. Provides the `step()` free function required by the `adaptive_steppable` concept.

| Type | Description |
|---|---|
| `ode::system<T,State>` | C++20 concept: callable `(stateIn, stateOut, time) -> void` |
| `ode::controlled_stepper<T,State>` | C++20 concept wrapping Boost.Odeint controlled stepper |
| `ode::Base<CES>` | CRTP-free base providing `step()`, `getDtDid()`, `logIntro()`, `logOutro()`, `stateIO()` |
| `ODE_EngineBoost<StateType>` | Concrete instantiation with RKCK54 stepper |
| `StepperDescriptor<S>` | Extensibility point: template variable specializations register stepper names |

Notable: `nextDtTryCorrectionFactor` patches a subtle numerical issue where tiny completion steps would corrupt the adaptive stepsize controller.

### 3.4 Trajectory.h

Generic time-evolution driver, entirely decoupled from physics. Concept hierarchy:

```
time_keeper            has getTime()
  adaptive_time_keeper   + getDtDid()
    adaptive_steppable     + step(t, deltaT) -> LogTree
      uniform_step           + advance(), temporalDataPoint()
        adaptive               uniform_step + adaptive_steppable
```

The `run()` function template is the primary user-facing entry point. All mode selection (`RLT`, `SFT`) is via template parameters — zero branching overhead in the hot loop.

#### 3.4.1 User-Facing Capabilities of `run()`

**Run length** — mutually exclusive; in Dt mode `NDt` takes precedence over `T` (in dc mode `NDt` is ignored):

| Parameter | Effect |
|---|---|
| `--T value` | Run until absolute simulation time T |
| `--NDt N` | Run for exactly N streaming intervals — avoids floating-point accumulation (`0.1×10 ≠ 1.0`), preferred for reproducible runs |

**Streaming frequency** — mutually exclusive; a nonzero `dc` takes precedence over `Dt`. Since `--dc` defaults to 10, Dt mode (and with it `--NDt`) needs an explicit `--dc 0`:

| Parameter | Effect |
|---|---|
| `--Dt value` | Stream every Dt simulated time units |
| `--dc N` | Stream every N ODE steps — denser output when dynamics are fast, sparse when slow; adaptive output density |

**Output destination:**

| Parameter | Effect |
|---|---|
| `--o filename` | Stream to file (appending); automatic continuation detection |
| (empty) | Stream to stdout |

**State serialization / checkpointing:**

| Parameter | Effect |
|---|---|
| `--sdf N` | Save state every N streaming events |
| `--initialFileName f` | Load initial state from file `f` |
| (automatic) | If output file exists and is non-empty, `.state` file is read and run resumes — no special flag needed |

State files are uncompressed binary archives (bzip2 compression was removed in May 2026). Each archive is self-contained with `SerializationMetadata` (typeID, trajectoryID, rank, protocolVersion) — prevents silent corruption from loading a wrong state file. `ARRAY_ONLY` mode allows loading raw array data (e.g. from NumPy) without full trajectory metadata.

**Output control** (`--streamSwitch`, 4-bit field, default `1011`). The value is parsed as a `std::bitset` string, so bit 0 is the **rightmost** character: the default `1011` means intro on, outro on, in-trajectory log off, initial-state serialization on.

| Bit | Effect |
|---|---|
| 0 | Write introduction: full command line, version, system parameters as JSON, column key |
| 1 | Write outro: ODE step counts, failed steps, derivative calls |
| 2 | Write in-trajectory log (ODE step events) |
| 3 | Serialize initial state before first step |

**Output precision:** `--precision N` — number of significant digits.

**Autostopping** — detects convergence to steady state via sliding window. **Currently broken** with the labelled Hana-tuple TDPs (see §7.1); pass `trajectory::observerNoOp` instead:

| Parameter | Effect |
|---|---|
| `--autoStopRepetition N` | Stop when TDP converges: same value N consecutive times within tolerance (0 = disabled) |
| `--autoStopEpsilonRel` | Relative tolerance for convergence |
| `--autoStopEpsilonAbs` | Absolute floor — components below this are not compared |

#### 3.4.2 HPC/HTC Production Workflow

A typical production run: compile once, submit a job array where each job varies `--o jobN.out` and some physical parameters. Each job streams data to its own file and checkpoints state to `jobN.out.state`. If a job is interrupted, resubmit with identical arguments — continuation is automatic. Once autostopping is restored (§7.1), it will terminate the job cleanly at steady state without a wall-time kill. For QJMC parameter studies, each job runs with `nTraj=0` (single trajectory); post-processing averages results across jobs.

The introduction written by `streamSwitch[0]` records the full command line and all system parameters — every output file is fully self-documenting for reproducibility.

#### 3.4.3 Internal Structure Notes

- `advance(traj, deltaT)` drives the trajectory for exactly `deltaT` by repeatedly calling `step()` — the loop handles the final partial step correctly
- `CommentingStream` prefixes log lines with `#`, interleaving data and logs in the same file while keeping them distinguishable
- `returnStreamedArray=true` accumulates all streamed TDPs into a returned `DataStream` — enables programmatic use of `run()` beyond pure file I/O
- The `observer` callback receives each TDP at streaming time — `AutostopHandlerGeneric` is the intended default (currently broken, §7.1); custom observers can be passed for e.g. online analysis

**Open items:**
- In-trajectory log in DC mode doesn't work for `dc > 1`
- TODO at line 349: unify observer, datastreamer, TDP calculation, `DataStream` return, and state streaming into a single `OBSERVER` concept
- `Read/WriteState` and `readViaSStream/writeViaSStream` will be partially superseded by the planned serialization modernization: array data moves to HDF5, small binary state (ODE/RNG) stays binary, metadata moves to JSON/HDF5 attributes — binary precision for continuation is a hard requirement that rules out JSON for state data

### 3.5 StochasticTrajectory.h

Physics-agnostic, lives in CPPQEDutils. Defines the `stochastic` concept and the `Ensemble` class.

The `stochastic` concept extends `uniform_step` with:
- `EnsembleAverageResult` — type of the accumulated ensemble average (e.g. `DensityOperator<RANK>`)
- `EnsembleAverageElement` — type of a single trajectory's contribution (e.g. `const StateVector<RANK>&`)
- `averaged(t)` — returns the trajectory's contribution to the average

`Ensemble<SingleTrajectory, TDP_Calculator>` is itself a `uniform_step` — it can be passed directly to `run()` without special-casing. Serial evolution is a plain loop over `trajs`. Recursive design: since `Ensemble` satisfies `stochastic`, it can be an element of a larger `Ensemble`.

`AverageTrajectories<SingleTrajectory>` is a customization point (traits class) for how trajectories are averaged. The default naive implementation assumes `EnsembleAverageElement` is additive and divisible by `double`. The `EnsembleAverageResult&` first argument is a pre-allocated output buffer — specializations use it to avoid repeated allocation (e.g. accumulating `|psi><psi|` dyads into a pre-allocated `DensityOperator`).

**Lazy ensemble averaging:** between streaming steps the buffer is never touched. At each streaming point, `temporalDataPoint` calls `tdpCalculator(getTime(e), averaged(e))`, which accumulates into the buffer and immediately passes it as `LDO<DensityOperator,RANK>` to `calculateExpectationValues`. The same `ExpectationValues` functional works identically for a single pure state, a single mixed state, or an ensemble average — `LDO` abstracts over all cases.

---

## 4. CPPQEDcore — Quantum Layer

### 4.1 Folder Structure

| Folder | Contents |
|---|---|
| `quantumdata` | Quantum state representations: `StateVector`, `DensityOperator`, `LazyDensityOperator` |
| `quantumoperator` | Sparse operators over tensor-product Hilbert spaces: `MultiDiagonal`, `SparseMatrix` |
| `quantumtrajectory` | Trajectory types satisfying `adaptive_steppable`: Master, QJMC (`Schrödinger.h` is a stale 2023 draft, see §7.3) |
| `quantumdynamics` | Physical system descriptors: `Hamiltonian`, `Liouvillian`, `ExactPropagator`, `ExpectationValues` |

### 4.2 quantumdata: LazyDensityOperator.h

The `LazyDensityOperator` (LDO) abstraction provides a unified interface for quantum averages from either a pure state (`StateVector`) or a mixed state (`DensityOperator`), without constructing the full density matrix in the pure case.

Implementation:

```cpp
template <template<size_t> typename TAG, size_t RANK>
struct LDO : MultiArrayConstView<dcomp, multiArrayRank_v<TAG<RANK>>>
```

- `TAG = StateVector`: underlying rank is `RANK`, indexing `_(psi,i,j)` computes `psi(i)*conj(psi(j))` on the fly
- `TAG = DensityOperator`: underlying rank is `2*RANK`, indexing `_(rho,i,j)` reads `rho(concatenate(i,j))` from storage
- `lazy_density_operator` concept checks that `_(t,i,j)->dcomp` and `_(t,i)->double` are well-formed

`partialTrace`: for `StateVector`, slices over dummy axes of the rank-R array. For `DensityOperator`, uses `extendedAxes` to slice only the diagonal in the dummy indices. Both accumulate via `fold_left_first` over a transformed `sliceRange`. Laziness is resolved entirely at compile time — no virtual dispatch, no type erasure.

### 4.3 quantumoperator: MultiDiagonal (from cppqed-multidiagonal)

`MultiDiagonal<RANK>` comes from the standalone library [cppqed-multidiagonal](https://github.com/WignerQuantumOptics/cppqed-multidiagonal), included as `<multidiagonal/MultiDiagonal.h>` and fetched at a pinned commit (§9). It is a sparse matrix representation storing an arbitrary set of diagonals, designed to replace the v2 `Tridiagonal` class. The key design motivation: `Tridiagonal` is not closed under composition (product of two tridiagonals is pentadiagonal in general), whereas `MultiDiagonal` is closed — composition produces new diagonals at all pairwise sums of offsets.

Data structure:

```cpp
using Offsets  = std::array<ptrdiff_t,RANK>;  // signed offset per axis: + upper, − lower, 0 main
using Diagonal = MultiArray<dcomp,RANK>;      // the actual diagonal values
std::map<Offsets, Diagonal> diagonals;        // ordered: deterministic accumulation order
```

The key of a diagonal is the number of quanta it moves along each axis. A diagonal at offset $d$ on an axis of dimension $N$ has length $N-|d|$ along that axis, stored as `d[k] = A(k+d⁻, k+d⁺)`. Subsystem operators are built directly into the map, e.g. `res.diagonals.emplace(Offsets{+1}, multiarray::fromStorage(v))` for the ladder operator $a$; σ⁺ = |1⟩⟨0| is the **lower** diagonal, `Offsets{-1}`.

Operations:
- `operator()` — application with accumulate semantics: `dpsidt += A psi`, iterating over each diagonal's valid index range
- `operator*` — direct product: concatenates the `Offsets` arrays, takes `directProduct` of the `Diagonal` arrays
- `operator|` — composition, rank-general: offsets add per axis; colliding offset pairs accumulate into the same diagonal
- `operator+=` / `-=` — addition: new diagonals may appear; existing ones accumulate element-wise
- `hermitianConjugate()` — negates all offsets and conjugates the diagonal values
- `twoTimesRealPartOf` / `twoITimesImagPartOf` — A + A† = 2 Re A (Hermitian) and A − A† = 2i Im A (anti-Hermitian)
- `multidiagonal::identity(dim)` — constructs identity operator
- `calculateAndCheckDimensions` — derives and validates Hilbert space dimension from diagonal extents+offsets

**Design decision: `MultiDiagonal` carries no frequency information.**

In v2, `Tridiagonal` fused static operator coefficients with interaction-picture phase rates (`freqs_`), and `propagate(t)` multiplied each element by `exp(freqs * t)` in place. This was necessary because `Tridiagonal` objects were passed between subsystems (via `dynamic_cast` on `particle::Ptr`) and had to carry picture information with them — the canonical example being `expINKX(particle::Ptr, nK)`, which would embed the particle's propagator frequencies into the coupling operator.

In v3 this cross-subsystem coupling is eliminated: the interaction picture is a property of the trajectory, not of individual operators. `ExactPropagator` / `UnaryDiagonalPropagator` handle all picture transformations globally at the trajectory level. `MultiDiagonal` is therefore a purely structural, time-independent object — the right abstraction boundary, and essential for the standalone library (§7.9).

When interaction-picture phase tabulation is genuinely needed (e.g. for the Particle subsystem, where $e^{inkx}$ couples momentum states with frequency differences $\omega_\text{rec}(k_j^2 - (k_j+n)^2)$), the correct representation is a **separate** `MultiDiagonal<RANK>` carrying the frequency values alongside the coefficient `MultiDiagonal<RANK>`:

```cpp
MultiDiagonal<RANK> alpha;  // static complex coefficients
MultiDiagonal<RANK> freqs;  // same sparsity structure, complex frequency values
```

with `operator()` evaluating `alpha(idx) * exp(freqs(idx) * t)` lazily. Composition of a `(alpha, freqs)` pair follows naturally: `alpha` composes via `operator|`, and `freqs` composes additively ($\delta_\text{result} = \delta_a + \delta_b$, from $e^{\delta_a t} \cdot e^{\delta_b t} = e^{(\delta_a+\delta_b)t}$). This pair lives one level up from `MultiDiagonal`, in the physics layer (CPPQEDelements), not in the mathematical library.

### 4.4 quantumtrajectory layer

All trajectory types satisfy `adaptive_steppable` and compose `ode::Base` with physics from `quantumdynamics`.

#### 4.4.1 QuantumTrajectory.h

Shared infrastructure:

- `initialTimeStep(h, d)` — estimates a good initial ODE timestep as 0.1 / max|H|psi>| by probing the Hamiltonian with a uniform state
- `TDP_DensityOperator<RANK,EV>` — the `tdpCalculator` used by `Ensemble`: wraps an `ExpectationValues` functional and calls `calculateExpectationValues` with `LDO<DensityOperator,RANK>`. This is the bridge between the ensemble averaging machinery and the `structure` layer.

#### 4.4.2 Master.h

Solves the Lindblad master equation directly on a `DensityOperator<RANK>`.

Key design insight — **reuse of state-vector operations for density matrix evolution, exploiting the H/i convention:**

- **Deterministic part** of Lindblad: applies `hamiltonian_ns::broadcast<compileTimeOrdinals<RANK>>` to rows of `ρ` (treating each row as a state vector), giving $\frac{H_\text{eff}}{i}\rho$ in `drhodt`, then calls `twoTimesRealPartOfSelf(drhodt)` which computes `X + X†` in place, yielding:

$$-i[H_\text{phys},\rho] - \frac{1}{2}\sum_k\left(J_k^\dagger J_k \rho + \rho J_k^\dagger J_k\right)$$

This is the full Lindblad equation **except** the recycling terms $J_k\rho J_k^\dagger$.

- **Recycling / jump part**: the $J_k\rho J_k^\dagger$ terms are handled separately by the `applySuperoperator` loop over each `Lindblad`. The decomposition exactly mirrors the QJMC split between coherent evolution and quantum jumps.
- **Exact propagator** `ρ → UρU†`: calls `exact_propagator_ns::broadcast` twice, with a `hermitianConjugateSelf` in between. $U$ is not necessarily unitary (e.g. for Mode, $U(t) = \text{diag}(1, e^{-zt}, \ldots)$ with complex $z$) — the same algebra works regardless.
- `rowIterationOffsets_` is precomputed once in the constructor via `calculateSlicesOffsets<compileTimeOrdinals<RANK>>` and shared by all three operations.

After each ODE step, `ρ` is symmetrized (`twoTimesRealPartOfSelf`) and renormalized to counteract numerical drift from Hermiticity and trace.

#### 4.4.3 QuantumJumpMonteCarlo.h

Implements single quantum-jump Monte Carlo trajectories. Templated on `qjmc::Algorithm` with two specializations.

`QuantumJumpMonteCarloBase<RANK,HA,EV,OE,RandomEngine>` holds all shared state: `ha`, `li`, `ev`, `psi`, `oe` (ODE engine), `re` (random engine), logging, serialization, `calculateRates()`, `performJump()`, `averaged()` (returns `const StateVector<RANK>&`).

```
EnsembleAverageElement = const StateVector<RANK>&
EnsembleAverageResult  = DensityOperator<RANK>
```

**Choice of random engine — critical for ensemble correctness:**

For parallel QJMC ensemble simulations, `std::mt19937` is **not recommended** — Mersenne Twister is designed for a single sequential stream, and seeding multiple instances with nearby seeds (e.g. `seed+i`) does not guarantee statistical independence between trajectories. Correlated random sequences can introduce systematic bias in ensemble averages.

The correct choice is **`pcg64`** (Permuted Congruential Generator), which supports provably independent parallel streams via a stream parameter. `randomutils::EngineWithParameters<pcg64>` encapsulates this: each trajectory gets a different `prngStream` ordinal, guaranteeing independence. `incrementForNextStream(p)` increments the stream ordinal between trajectory initializations in `makeEnsemble`. `Xoshiro256PlusPlus` is another high-quality option with good parallel properties.

`std::mt19937` has no `Pars` specialization in `Random.h` — `pcg64` is the intended production engine.

**`Algorithm::integrating`** (QuTiP-style):
- Evolves under non-Hermitian $H_\text{eff}/i$ (norm decays continuously, encoding integrated jump probability)
- Samples random threshold `normAt_`; when `normSqr(psi) < normAt_`, a jump occurs
- Bisects to find the exact jump time via `bisect_` struct (bracketing interval `[t0,t1]` with norms `[n0,n1]`)
- More accurate — jump time is exact to within `normTol`

**`Algorithm::stepwise`** (Comp. Phys. Comm. 238:88 2019):
- Evolves under physical Hermitian H, renormalizes after each step
- Computes rates explicitly, samples jump probabilistically as `sampleRandom()/getDtDid() < totalRate`
- Adjusts timestep so `dp = totalRate * dtTry ≈ dpLimit`
- Simpler but introduces discretization error proportional to `dpLimit`
- **Currently broken** — stale references to `q.ex` and `q.qsd` from an earlier refactor

`AverageTrajectories` specialization for QJMC:
```cpp
rho = dyad(trajs[0].psi, trajs[0].psi);
for each remaining traj: addTo(traj.psi, rho);
rho /= nTraj;
```
Accumulates `|ψk⟩⟨ψk|` directly into the pre-allocated `DensityOperator` buffer — no intermediate density matrices allocated.

### 4.5 structure layer (→ `quantumdynamics`)

All structure elements are concept-constrained callables, not class instances. Design consistency across the layer:
- **Hamiltonian, ExpectationValues**: compile-time Hana sequences — structure fixed at compile time
- **Liouvillian**: runtime `std::vector` — number of jump operators may vary; `std::join_view` limitations make Hana composition impractical
- **ExactPropagator**: per-subsystem by nature, no sequence needed — `broadcast()` with `retainedAxes` handles subsystem application

#### 4.5.1 Hamiltonian.h

Recursive compile-time Hana sequences of functionals. Three levels of time dependence:

| Concept | Signature |
|---|---|
| `time_independent_functional` | `H(psi, dpsidt)` |
| `one_time_dependent_functional` | `H(t-t0, psi, dpsidt)` — interaction picture |
| `two_time_dependent_functional` | `H(t, psi, dpsidt, t0)` — full Schrödinger picture |

**Important:** the functionals apply $H/i$ to the state vector, not $H$ itself. Furthermore $H$ is **not necessarily Hermitian** — in the QJMC and Master contexts it is the full effective non-Hermitian Hamiltonian $H_\text{eff} = H_\text{phys} - \frac{i}{2}\sum_k J_k^\dagger J_k$, which includes the no-jump backaction from the dissipator. The $H/i$ convention is what makes the algebra click consistently across all trajectory types (pure state, master equation, QJMC).

Similarly, the exact propagator $U$ is **not necessarily unitary** — for Mode with complex frequency $z = \kappa(2n_\text{th}+1) - i\delta$, $U(t) = \text{diag}(1, e^{-zt}, \ldots)$ decays as well as rotates.

`applyHamiltonian` recurses through Hana tuples with `hana::for_each`, dispatching via `if constexpr`. `NoOp` is a zero-overhead sentinel. `vectorize()` constructs an explicit sparse matrix by probing with unit vectors.

**Open item:** `compose()` for interaction-type elements (`H_12` acting on a bipartite system) is declared but not implemented.

#### 4.5.2 ExactPropagator.h

Represents the free evolution operator U in the interaction picture — diagonal in the chosen basis, per-subsystem, and not necessarily unitary. `applyPropagator` **replaces** `|psi>` with `U|psi>` (unlike Hamiltonian which **adds** to `d|psi>/dt`). No recursive concept needed: each subsystem contributes exactly one diagonal propagator.

`UnaryDiagonalPropagator`: caches `std::valarray<dcomp>` diagonal, recomputes only when `t` or `t0` changes. Operator bodies are declared but not yet implemented — **high priority open item** (see section 7.2).

`broadcast<retainedAxes>()` applies a propagator to a subsystem by iterating `sliceRange`.

**The continuous interaction-picture reset — a key numerical method:**

When an `ExactPropagator` is supplied, the framework assumes the Hamiltonian is defined in the interaction picture defined by $H_\text{free}$ (the exactly exponentiated part), while expectation values, jump operators, and rates are defined in the Schrödinger picture. After each ODE step from $t_0$ to $t$, `applyPropagator` transforms the state back to the Schrödinger picture via $U(t,t_0)$, and then $t_0 \leftarrow t$ is updated. This **continuously resets the interaction picture origin** to the current time.

This is not mere bookkeeping — it is the mechanism that prevents numerical blow-up in the non-unitary case. Without resetting, the interaction-picture operator matrix elements accumulate exponential factors $e^{\pm nz\tau}$ over the full simulation time $T$. For a cavity mode with decay $\kappa$, Fock cutoff $N$, and simulation time $T$, the ratio between largest and smallest matrix elements would grow as $e^{N\kappa T}$ — for $N=100$, $\kappa=10$, $T=10$ this is $e^{10000}$, far beyond double precision.

With continuous reset, the same ratio is bounded by $e^{N\kappa\Delta t}$ where $\Delta t$ is the adaptive ODE stepsize — implicitly controlled by the error tolerance. The representation stays well-conditioned throughout, consistent with the physical Fock space truncation assumption. For large Hilbert spaces with significant decay — exactly the regime C++QED targets — this is likely the difference between a correct simulation and a numerically corrupted one.

#### 4.5.3 Liouvillian.h

Dissipative evolution via Lindblad operators. Each `Lindblad<RANK>` bundles: `label` (string), `jump` (`Jump<RANK>` variant), `rate` (`Rate<RANK>` variant), `superoperator` (`Superoperator<RANK>` variant). All dispatch via `std::visit` with `overload`.

Utilities:
- `rateFromJump` — derives rate as `||J|psi>||^2`, avoiding separate specification
- `superoperatorFromJump` — derives full superoperator action on `rho` from jump operator via `sliceRange` row iteration
- `concatenate` — fold-based composition of multiple `Liouvillian`s into one vector

#### 4.5.4 ExpectationValues.h

Follows the same recursive Hana-sequence pattern as Hamiltonian. Functionals take `lazy_density_operator` (not `StateVectorConstView` directly), so they work identically for pure and mixed states.

| Concept | Signature |
|---|---|
| `time_independent_functional` | `ev(rho) -> temporal_data_point` |
| `time_dependent_functional` | `ev(t, rho) -> temporal_data_point` |

`NoOp` carries a label and returns `hana::make_tuple()` as an empty `temporal_data_point` — systems with no observables participate correctly in data streaming without special-casing.

---

## 5. Parameter Machinery — Pars.h / Pars.cc

The parameter system is built on top of [popl](https://github.com/badaix/popl) and extends it significantly. It plays a central role in the highest-level user interface — every `Pars` struct in the framework uses it.

### 5.1 What popl provides

Basic command-line parsing: register options with names, descriptions, defaults, and bindings, then call `parse()`. C++QED wraps this in `optionParser()` (which adds `--help` and `--version` switches) and `parse()` (which also records the full command line in `parsedCommandLine` for logging).

### 5.2 What the Framework Adds

**Hierarchical namespacing via `mod`:**

Every `Pars` constructor accepts a `mod` string suffix appended to every option name. So `mode::Pars` constructed with `mod="1"` registers `--cutoff1`, `--delta1`, etc. Two modes in a composite system use `mod="1"` and `mod="2"`, giving `--cutoff1`, `--cutoff2` — no name collisions, no special handling needed.

**Automatic JSON serialization via `to_json_converter`:**

`to_json_converter` is a `std::list<std::function<void(json::object&)>>` accumulating closures — one per parameter, each capturing the option name and a reference to the bound variable. Calling `jsonize()` walks the list and builds a `json::object` of current parameter values. Parameters registered with `noJSON` are excluded — e.g. initial conditions in `mode::Pars` that don't belong in the system descriptor JSON.

**The `_()` helper for tuple-based registration:**

```cpp
_("cutoff", "Fock space cutoff", 10, cutoff)
```
returns a 4-tuple. A `noJSON` variant returns a 5-tuple. These are unpacked by `add_dispatch` via a `multilambda` fold — titles (plain strings), normal parameters, and `noJSON` parameters are all handled uniformly in a single variadic call.

**The result — a single `add(...)` call does three things at once:**
- Registers all options with popl (with `mod` namespacing)
- Wires up JSON serialization for all non-`noJSON` parameters
- Inserts section titles in the `--help` output

### 5.3 Example: mode::Pars

```cpp
struct Pars : ::parameters::JSONizable
{
  size_t cutoff, initFock;
  double delta, omegaKerr, kappa, nTh;
  dcomp eta, init;

  Pars(popl::OptionParser& op, std::string mod="")
  {
    using namespace ::parameters;
    add(mod, op, tjc, "Mode",
        _("cutoff",    "Fock space cutoff",              10,       cutoff),
        _("delta",     "detuning",                       -10.,     delta),
        _("omegaKerr", "Kerr constant",                  0.,       omegaKerr),
        _("eta",       "drive amplitude",                dcomp(0), eta),
        _("kappa",     "decay rate",                     10.,      kappa),
        _("nTh",       "thermal photon number",          0.,       nTh),
        _("initFock",  "Fock state initial condition",   0,        initFock, noJSON),
        _("init",      "Coherent state initial condition", dcomp(0), init,   noJSON)
    );
  }
};
```

- `tjc` is inherited from `JSONizable` — it accumulates the serialization closures
- `"Mode"` is the section title in `--help` output
- `initFock` and `init` use `noJSON` — they control initial conditions but are not part of the system descriptor passed to `make()`
- `make(const Pars& p)` calls `p.jsonize()` to produce the `json::object descr` passed into the tuple

### 5.4 Hierarchical Composition

Since every SubSystem `Pars` struct takes `(popl::OptionParser&, std::string mod)`, and `trajectory::Pars<BASE>` / `ode::Pars<BASE>` take `(popl::OptionParser&)` and chain constructors via the `BASE` template parameter, a composite simulation's full parameter set is assembled by inheritance plus members:

```cpp
using TrajPars = trajectory::Pars<qjmc::Pars<qjmc::Algorithm::integrating, RandomEngine>>;

struct MyPars : TrajPars
{
  mode::Pars pm;
  qbit::Pars pq;
  MyPars(popl::OptionParser& op)
    : TrajPars{op}, pm{op,"m"}, pq{op,"q"} {}
};
```

A single `parse(op, argc, argv)` call then populates every parameter across all subsystems, with `--cutoffm`, `--deltaq` etc. namespaced automatically.

---

## 6. CPPQEDelements — Physical Systems and Interactions

Concrete physical systems and interactions built on top of CPPQEDcore. Unary systems are called **SubSystems**; interactions are binary (or higher) functionals coupling them. Each SubSystem follows a consistent pattern:

- Operators — as `MultiDiagonal` or direct lambdas
- `hamiltonian(...)` — returns a `hamiltonian<RANK>`-compatible callable or `MultiDiagonal`, applying $H_\text{eff}/i$
- `Liouvillian<RANK>` built from `Lindblad` structs
- `expectationValues` lambda over `lazy_density_operator`
- `make(...)` — returns a plain tuple `(cutoff, hamiltonian, liouvillian, expectationValues, descr)` — a SubSystem is not an object, it's a tuple of constituents
- `Pars` struct for command-line/JSON parameter parsing
- Initial state constructors

### 6.1 Mode.h / Mode.cc — Harmonic Oscillator (Bosonic Mode)

Quantum harmonic oscillator truncated to Fock space of dimension `cutoff`. The canonical SubSystem.

Operators via `MultiDiagonal` — closure under composition pays off immediately:
```cpp
aOp(cutoff)    // annihilation: superdiagonal with sqrt(n+1)
aDagOp(cutoff) // creation: hermitianConjugateOf(aOp)
nOp(cutoff)    // number: aDagOp | aOp  — composition via operator|
```

Hamiltonian ($H_\text{eff}/i$, complex frequency $z = \kappa(2n_\text{th}+1) - i\delta$):
```cpp
ham -= z * nOp(cutoff);                     // detuning + decay (complex frequency)
ham += (omegaKerr/1.i) * (aDag|aDag|a|a);  // Kerr nonlinearity a†a†aa
ham += twoITimesImagPartOf(eta*aDag);        // coherent drive H = i(η a† − η* a), so H/i = η a† − η* a
```
The Kerr term `aDag|aDag|a|a` is four chained compositions — impossible with `Tridiagonal`, trivial with `MultiDiagonal`.

Lindblads: `photonLoss(kappa, nTh)`, `photonGain(kappa, nTh)` — each bundles jump, rate, and superoperator.

`expectationValues`: computes photon number $\langle\hat{n}\rangle$, $\langle\hat{n}^2\rangle$, and $\langle a\rangle$ via `_(rho,n)` and `_(rho,n,n-1)` — works identically for pure and mixed states via LDO.

Initial states: `coherent(alpha, cutoff)` (with Stirling approximation for high Fock states), `fock(n, dim)`.

**Open item:** `UnaryDiagonalPropagator<> propagator(dcomp z)` is declared but commented out — blocked by unimplemented `UnaryDiagonalPropagator` operator bodies. In v2, Mode had exact propagation for the free oscillator evolution $U(t) = \text{diag}(1, e^{-zt}, e^{-2zt}, \ldots)$, removing fast oscillations from the ODE and allowing much larger timesteps. High priority to restore.

### 6.2 Qbit.h / Qbit.cc — Two-Level System

Two-level system (qubit/spin-½). Hilbert space dimension is always 2 — operators implemented as direct lambdas rather than `MultiDiagonal` (overkill for 2×2).

Hamiltonian as composed lambdas:
```cpp
diagonalH(z)      // dpsidt(1) -= z*psi(1)
offDiagonalH(eta) // off-diagonal drive: H = i(η σ+ − η* σ-)
fullH(z, eta)     // composition by lambda capture
```

Lindblads: `loss(gamma_m)`, `gain(gamma_p)`, dephasing (commented out).

`expectationValues`: population of ground state and polarization $\langle\sigma_+\rangle = \rho_{01}$ (note: Mode's ladder-operator average is $\langle a\rangle$, the opposite convention).

`UnaryDiagonalPropagator<> propagator(dcomp z)` **is implemented** in `Qbit.cc` (but not yet used by `qbit::make`):
```cpp
return {"diagonalPropagator", 2, [=] (double t, auto& d) {d[0]=1; d[1]=exp(-z*t);}};
```
giving $U(t) = \text{diag}(1, e^{-zt})$ — but cannot be used yet because `UnaryDiagonalPropagator::operator()` bodies are unimplemented stubs. This makes `UnaryDiagonalPropagator` the highest-priority open item — it would unlock exact propagation for both Qbit and Mode.

`multidiagonal` sub-namespace defines `splus`, `sminus`, `sx`, `sy`, `sz` as `MultiDiagonal<1>` — with a TODO to define a leaner `Sigma` sparse type.

### 6.3 Composite.h — Broadcaster

`Composite.h` lives in `CPPQEDcore/quantumdynamics/`. `Broadcaster<RANK, ra...>` lifts a rank-`|ra|` system descriptor to rank-`RANK` by iterating `sliceRange<retainedAxes>` over the dummy axes. Template parameters `ra...` specify which axes the subsystem occupies in the composite system. This is the key composition primitive for building composite systems from SubSystems.

```cpp
Broadcaster<2,0> bQ{psi.extents};  // qubit on axis 0
Broadcaster<2,1> bM{psi.extents};  // mode on axis 1
```

`Broadcaster` handles all three descriptor types via overloaded `operator()`:

- **Hamiltonian / ExpectationValues** — templated `operator()` with `if constexpr` dispatch: zips `sliceRange<retainedAxes>` over `psi`+`dpsidt` (Hamiltonian) or uses `partialTrace` (ExpectationValues)
- **`Lindblad<RRANK>` → `Lindblad<RANK>`** — dedicated `operator()`:
  - `jump`: iterates `sliceRange<retainedAxes>` over `psi`, applies `applyJump` to each slice
  - `rate`: uses `partialTrace` with `calculateRate` as the per-slice function
  - `superoperator`: lazily populates `matrixOffsets` on first call using the `extendedAxes` diagonal slice pattern, then zips `sliceRange<extendedAxes<retainedAxes,RANK>>` over both `rho` and `drhodt`
- **`Liouvillian<RRANK>` → `Liouvillian<RANK>`** — transforms each `Lindblad` in the vector via the above

Status: largely complete. May have minor issues to resolve before compiling cleanly.

### 6.4 JaynesCummings.h — Cavity-Atom Interaction

Interaction between a cavity mode and a two-level atom: $H_{JC} = i\left(g^* a^\dagger \sigma_- - g\, a \sigma_+\right)$ (the same phase convention as the drive terms of Mode and Qbit).

```cpp
template <size_t i, size_t j>
auto hamiltonian(size_t cutoff, dcomp g)
```

Template parameters `i`, `j` are the qubit levels coupled at compile time — covering both standard JC and multi-level generalizations without runtime cost.

Implementation uses `Placeholder` slicing: `psi(i,_)` extracts the rank-1 cavity slice at qubit level `i`, then applies the rank-1 `MultiDiagonal` to it.

The interaction is a `time_independent_functional<2>` — slots directly into a Hana tuple with Mode and Qbit Hamiltonians for `applyHamiltonian`.

**Open items:**
- `label()` for the lambda is commented out — ADL doesn't pick it up for lambda types
- A purely binary `MultiDiagonal` design is noted as an alternative — would require `compose()` in `Hamiltonian.h`
- Generalizing `Placeholder` slicing beyond `RANK==2` needed for higher-rank composite systems

### 6.5 Particle.h / ParticleCavity.h — Motional Degree of Freedom

The motional degree of freedom (MDoF) subsystem was present in C++QED from its earliest versions, originating in PhD work on mobile atoms in cavity fields. It is a distinguishing feature of the framework — competitors rarely provide first-class support for spatially extended quantum systems. Up to 2D motion can be represented as two independent MDoF subsystems. These files are reviewed from v2 and have not yet been migrated to v3.

#### Physics and representation

A particle moving in a periodic potential (e.g. a standing-wave optical lattice or cavity field) is most naturally represented in **momentum space**. The kinetic Hamiltonian $H_\text{kin} = \omega_\text{rec}\, \hat{k}^2$ is diagonal in momentum eigenstates $|k_j\rangle$, so it is exactly exponentiated:

$$U_\text{kin}(t,t_0) = \exp\!\bigl(-i\,\omega_\text{rec}\,k_j^2\,(t-t_0)\bigr)$$

This makes the interaction picture essential: without it, the ODE timestep would be constrained by $1/(\omega_\text{rec} k_\text{max}^2)$, which for `fin=6` (dimension $2^6 = 64$, $k_\text{max} = 32$ in units of $\Delta k$) is already $k_\text{max}^2 = 1024$ times smaller than the recoil timescale $1/\omega_\text{rec}$. The exact propagator removes the kinetic term entirely from the ODE.

The optical potential introduces $\sin(nkx)$ or $\cos(nkx)$ terms. In momentum space, $e^{inkx}$ acts as a pure shift: $e^{inkx}|k_j\rangle = |k_j + n\rangle$. This is a `MultiDiagonal<1>` with a single off-diagonal band of ones at offset $+n$.

#### `particle::Spatial`

A pure data helper: for a given finesse `fin` (dimension $= 2^\text{fin}$, i.e. `1<<fin`) and lattice spacing $\Delta k$, it precomputes the coordinate grids

$$x_j = j\,\Delta x - x_\text{max}, \quad k_j = j\,\Delta k - k_\text{max}$$

stored as `DArray<1>` arrays. The two grids are related by the discrete FFT, encoding the canonical commutation $[X,K]=i$.

#### `particle::Exact` → `UnaryDiagonalPropagator`

In v3 this becomes a `UnaryDiagonalPropagator` with `updateDiagonal` lambda:

```cpp
[=](double tau, auto& d) { for (size_t j=0; j<dim; j++) d[j] = exp(-1i*omrec*k(j)*k(j)*tau); }
```

This is one of the clearest motivations for completing the `UnaryDiagonalPropagator` implementation (§7.2): the particle subsystem's entire numerical advantage depends on it.

#### `expINKX`, `cosNKX`, `sinNKX` — coupling operators

In v3:
- `expINKX(nK, dim)` → `MultiDiagonal<1>` with a single diagonal at signed offset `+nK`, filled with ones
- `cosNKX(nK, dim)` → `(expINKX(nK,dim) + expINKX(-nK,dim)) / 2`
- `sinNKX(nK, dim)` → `(expINKX(nK,dim) - expINKX(-nK,dim)) / (2i)`

No frequency information is embedded in these operators — that is the responsibility of the `UnaryDiagonalPropagator` at the trajectory level.

#### `mfComposition` — eliminated by `MultiDiagonal`

In v2, `mfComposition` was an 80-line switch table computing $f^*(n_1 kx)\cdot g(n_2 kx)$ for all combinations of $\{\sin, \cos, e^+, e^-\}^2$, restricted to the case $|n_1|=|n_2|$ (where the result is tridiagonal). This restriction and the entire table are artifacts of `Tridiagonal`'s lack of closure under composition. In v3, the product of any two `MultiDiagonal` operators is computed by `operator|` — the result has diagonals at all pairwise offset sums, correctly and automatically, with no restriction on $|n_1|$ vs $|n_2|$.

#### Initial conditions

Three families:
- **`wavePacket`**: Gaussian $\psi(x)\propto e^{-(x-x_0)^2/(4\sigma^2)+ik_0 x}$, FFT'd to momentum space if `kFlag=true`
- **`hoState`**: harmonic oscillator eigenstate using Hermite polynomials (Boost.Math)
- **`init()`**: dispatcher selecting wavepacket vs HO eigenstate based on `hoInitn`

#### `ParticleCavity` interactions — design note

In v2, `expINKX(particle::Ptr, nK)` detected (via `dynamic_cast<const particle::Exact*>`) whether the particle carried an exact propagator, and if so, embedded the corresponding frequency differences into the `Tridiagonal`'s own frequency data. This cross-subsystem coupling is **not reproduced in v3**, and deliberately so. The mechanism introduced a hidden dependency: operator correctness depended on a runtime property of a partner object, with no compile-time enforcement. The v3 architecture places the interaction picture entirely at the trajectory level, so coupling Hamiltonians are written once for the interaction picture and remain correct regardless. The price is that the picture is no longer a self-adapting property of individual operators — it is an explicit, visible property of the trajectory construction. This is the right trade-off.

The `ParticleAlongCavity` and `ParticleOrthogonalToCavity` interactions migrate to v3 as rank-2 `hamiltonian<2>` functionals built from `Broadcaster`-lifted rank-1 operators. The `isSpecialH_` flag and the `addContribution_v` override both become unnecessary — `MultiDiagonal::operator|` handles arbitrary mode function products correctly and generally.

### 6.6 Example: examples/jc_qjmc.cc

A complete, **verified working** highest-level executable (first successfully compiled and run April 2026, after a 4-year development hiatus, on Windows 11 via WSL2/Ubuntu 24.04). Demonstrates the full composition pattern.

```cpp
// axis 0 = qubit, axis 1 = mode

#include "Composite.h"
#include "JaynesCummings.h"
#include "Mode.h"
#include "Qbit.h"
#include "QuantumJumpMonteCarlo.h"
#include "Pars.h"
#include "pcg_random.hpp"

using namespace cppqedutils;
using namespace quantumdata;
using namespace quantumtrajectory;
using namespace quantumdynamics;

// pcg64 is the correct engine for parallel QJMC ensemble simulations:
// std::mt19937 is NOT recommended — nearby seeds do not guarantee independent streams.
// pcg64 supports provably independent parallel streams via the prngStream parameter.
using RandomEngine = pcg64;

int main(int argc, const char* argv[])
{
  auto op = optionParser();

  // All parameters registered into a single op — no custom Pars struct needed
  cppqedutils::trajectory::Pars<qjmc::Pars<qjmc::Algorithm::integrating, RandomEngine>> pt{op};
  qbit::Pars pq{op, "q"};
  mode::Pars pm{op, "m"};
  double g;
  add(op, parameters::_("g", "Jaynes-Cummings coupling strength", 1., g));

  parse(op, argc, argv);

  // Initial state: rank-2 tensor product psi_ij = psiQ_i * psiM_j
  StateVector<2> psi{ qbit::init(pq) * mode::init(pm) };

  // Broadcasters: lift rank-1 descriptors to rank-2
  Broadcaster<2,0> bQ{psi.extents};  // qubit on axis 0
  Broadcaster<2,1> bM{psi.extents};  // mode on axis 1

  // Unpack subsystem descriptors
  auto [cutoffQ, haQ, liQ, evQ, descrQ] = qbit::make(pq);
  auto [cutoffM, haM, liM, evM, descrM] = mode::make(pm);

  // Compose: Hamiltonian as Hana tuple, Liouvillian by manual concatenation
  // (quantumdynamics::concatenate has a constexpr bug with non-constexpr arguments)
  auto ha = hana::make_tuple( bQ(haQ), bM(haM),
                               jaynescummings::hamiltonian<0,1>(pm.cutoff, g) );
  auto li = bQ(liQ);
  for (auto& l : bM(liM)) li.push_back(std::move(l));
  auto ev = hana::make_tuple( bQ(evQ), bM(evM) );

  auto descr = json::object{
    {"system","Jaynes-Cummings"}, {"qubit",descrQ}, {"mode",descrM}, {"g",g}
  };

  // observerNoOp must be passed explicitly — template default argument deduction
  // does not work for AH, so the default in quantumtrajectory::run is dead code
  quantumtrajectory::run<qjmc::Algorithm::integrating, ODE_EngineBoost>(
    ha, li, ev, descr, std::move(psi), pt, trajectory::observerNoOp);
}
```

**Key observations:**
- The entire composite system is assembled in ~10 lines of substantive code
- `Pars` instances registered directly into `op` — no custom struct needed
- `operator*` builds the rank-2 tensor product initial state from rank-1 subsystem states
- `Broadcaster` does all axis-embedding; `JaynesCummings` slots in as a third Hana element
- `quantumdynamics::concatenate` bypassed due to constexpr bug — manual `push_back` used instead
- `trajectory::observerNoOp` must be passed explicitly (template deduction limitation)
- `pcg64` used throughout for parallel-safe PRNG with independent streams

**Known issues encountered during first compilation (April 2026):**
- `qjmc::Pars` requires two template arguments: `<Algorithm, RandomEngine>`
- `quantumdynamics::concatenate` fails with non-constexpr arguments — see bug in section 7.1
- `quantumtrajectory::run` default observer argument is dead code — must pass explicitly
- `AutostopHandlerGeneric` incompatible with Hana tuple TDPs — bypassed via `observerNoOp`

---

## 7. Open Development Threads

### 7.1 Known Bugs

- **`AutostopHandlerGeneric`** in `Trajectory.h`: uses `std::ranges::fill`, arithmetic, and `abs()` directly on `TDP` — assumes `TDP` is a homogeneous numeric range like `std::valarray<dcomp>`. Broken when `TDP` is a Hana tuple of labelled `slv_ns::_` values (the current representation). This is a half-finished update: the TDP type was upgraded from `std::valarray` to a labelled Hana tuple, but `AutostopHandlerGeneric` was not updated to match. Workaround: pass `observerNoOp` explicitly to `run()`. Fix: either make `AutostopHandlerGeneric` work with Hana tuples via `hana::for_each`, or flatten TDP to `std::valarray<double>` before passing to it.
- **`quantumtrajectory::run` default observer**: the default argument `AH&& observer = trajectory::observerNoOp` is dead code — C++ template deduction does not consider default arguments, so `AH` cannot be deduced and the caller must always pass the observer explicitly. The default should be removed to avoid confusion.

### 7.2 High Priority: Unimplemented

- **`UnaryDiagonalPropagator` operator bodies** — declared in `ExactPropagator.h` but not implemented. Blocks exact propagation for `Qbit` (`propagator()` defined, but not yet used by `make()`), `Mode` (needs `propagator()` reinstated), and the future Particle subsystem (the primary numerical motivation for MDoF in the framework). Significant numerical advantage: removes fast oscillations from the ODE, allowing much larger timesteps. For Particle with `fin=6`, the kinetic frequencies span a factor of $k_\text{max}^2 = 1024$ — without exact propagation the ODE timestep is constrained accordingly.
- **QJMC stepwise algorithm** — stale references to `q.ex` and `q.qsd` from earlier refactor, currently broken.
- **`makeEnsemble`** in `QuantumJumpMonteCarlo.h` — commented out, ensemble runs currently disabled.

### 7.3 Other Unimplemented / Stub

- `broadcastFor` in `SliceIterator.h` — generic functional primitive for applying `f` across matching slices of multiple arrays
- `compose()` in `Hamiltonian.h` — interaction Hamiltonian composition for binary/composite systems
- `Placeholder` slicing generalization beyond `RANK==2`
- `MultiDiagonal operator|` composition — only rank-1 implemented
- `InitializeEnsembleFromArrayOnlyArchive` for QJMC — just throws; design question unresolved
- Schrödinger equation trajectory — not yet implemented (`Schrödinger.h` is a 2023 draft that still includes the removed `QuantumSystemDynamics.h`)
- **Particle / MDoF subsystem** — v2 files reviewed (§6.5), migration not yet started

### 7.4 Design Improvements Under Consideration

- `MultiArrayBase<StorageType>` refactor: single template parameterized on storage type instead of inheritance
- `MultiIndexRange` as a proper `std::ranges::range`
- Unification of observer/streamer/state concerns in `Trajectory.h` into a single `OBSERVER` concept
- JSON interface to numpy (may make some Boost.Serialization `Read/WriteState` machinery redundant)
- `Sigma` sparse type for two-level operators, leaner than `MultiDiagonal<1>` for `Qbit`
- Compile-time leading dimensions for `StateVector` / `DensityOperator` — see section 10

### 7.5 Naming Changes and Consistency

- **`structure/` → `quantumdynamics/`** (done 2026-10, together with the switch to cppqed-multidiagonal): the old name was opaque; `quantumdynamics` parallels the other CPPQEDcore folders (`quantumdata`, `quantumoperator`, `quantumtrajectory`) and accurately describes the contents
- **`Free` → `SubSystem`** (directory `CPPQEDelements/subsystems/`, done 2026-10): "Free" has a specific field-theory meaning (quadratic Hamiltonian, exactly solvable) that doesn't hold for all unary systems (e.g. a driven qubit is not strictly free); `SubSystem` is more neutral and descriptive of the software role
- **Acronym casing**: camelcase is the chosen convention throughout, but acronyms are not handled consistently. Recommended rule: **capitalise only the first letter of an acronym** in camelcase contexts (`Qjmc`, `Ode`, `Ldo`, `Tdp`). Cosmetic — defer until the codebase is otherwise stable.

### 7.6 Planned: Test Suite

A clean, comprehensive test suite is a prerequisite for both confident development and the SciPost paper.

#### 7.6.1 v2 Test Suite Architecture

The v2 test suite (`Testing/`) used a Python driver (`testdriver.py`, ~700 lines) orchestrating CTest, with configuration-file-driven test descriptions. The following sub-directories existed:

- **`compilefail/`** — tests that should fail to compile *with a specific expected error message*. The `CompileTarget` class invokes CMake to build a target, expects a non-zero return code, and checks that the expected string appears in stdout or stderr. Useful for verifying SFINAE/concept constraints, deleted operations, and invalid template instantiations.
- **`boost/`** — Boost.Test unit tests, compiled and run as standard CTest tests.
- **`run/`** — runs simulation scripts and checks they succeed for all runmodes (`single`, `master`, `ensemble`). Verifies basic executability rather than physics correctness.
- **`continue/`** — runs a simulation, then resumes it from the saved state file with identical arguments, checking that continuation works and output is consistent.
- **`expected/`** — compares simulation output (trajectory data files and binary state files) against stored reference files. The `Verifier` class supports both full trajectory comparison (`verify=full`) and final-outcome-only comparison (`verify=outcome`). Changes in simulation output trigger test failures.
- **`physics/`** — high-level correctness tests: runs a full simulation in a limit where approximate analytic theory applies (Rabi oscillations, free mechanical oscillations, thermal equilibration, etc.) and checks numerical agreement against the analytic prediction. These tests verify *physical correctness*, not just code correctness.

#### 7.6.2 Test Suite for v3 / cppqed-multidiagonal

The `physics/` layer belongs entirely to CPPQEDcore and CPPQEDelements — it is out of scope for the standalone `cppqed-multidiagonal` library.

For `cppqed-multidiagonal`, the relevant layers are:

- **`boost/`** (unit tests): natural test cases:
  - `MultiArray`: construction, indexing, move semantics, JSON round-trip, `copyInit` extent mismatch
  - `SliceIterator`: `filterIn`/`filterOut` on known arrays, `calculateSlicesOffsets` against hand-computed values, `sliceRange` iteration count, `transpose` permutation
  - `MultiDiagonal`: `identity(N) | A == A`, `A | identity(N) == A`, `(A|B)|C == A|(B|C)`, `(A*B)† == A†*B†`, `hermitianConjugateOf(hermitianConjugateOf(A)) == A`, application of `aOp`/`aDagOp`/`nOp` to Fock states against analytic results

- **`compilefail/`**: a minimal Python script (~50 lines, extracted from v2's `CompileTarget`) for expected-failure-with-message-check tests. Natural targets:
  - Deleted copy constructor/assignment on `MultiArray` and `MultiDiagonal`
  - `consistent<retainedAxes, RANK>` violations in `sliceRange`
  - `scalar` concept violations in `operator*=`

- **`expected/`** (lightweight, JSON-based): apply known operators to known Fock states, serialize the result via `toStringJSON`, compare against stored JSON reference files. Plain `diff` suffices — no Python/numpy required. This is the primary justification for retaining JSON serialization in the standalone library.

#### 7.6.3 Test Suite for CPPQEDcore / CPPQEDelements (v3)

The full v2 architecture ports naturally to v3 for the higher layers:

- **`run/`** and **`continue/`**: verify the trajectory driver, checkpointing, and continuation
- **`expected/`** (full): trajectory output files compared via `numpy.allclose`; binary state files compared via a HDF5 reader once serialization is modernized
- **`physics/`**: same structure as v2 — `FunctionComparer` for analytic limits, specific targets:
  - Mode: coherent state evolution under drive, thermal equilibration under photon loss
  - Qbit: Rabi oscillations under coherent drive
  - Jaynes-Cummings: vacuum Rabi splitting, collapse-and-revival in photon number
  - Particle: wavepacket free propagation (exact vs Schrödinger picture agreement), mechanical oscillations in lattice potential

Benchmarks (for Paper 2): scaling of walltime with Hilbert space dimension, comparison against QuTiP for equivalent systems — run separately from the correctness suite, not part of CTest.

### 7.7 Planned: API Documentation

Doxygen-style documentation to be added to all headers. Style questions to resolve:
- Pure Doxygen (`/** */` with `\brief`, `\tparam`, `\param`) vs lighter `///` line comments (already used in places)
- Scope: public API only vs internal namespaces (e.g. `slice_iterator`)
- Convention for documenting C++20 concepts, which Doxygen doesn't natively render well

### 7.8 Planned: Build System and Documentation

#### Build system: CMake

Recommended modernization rather than replacement:
- Migrate to **modern target-based CMake** (no global `include_directories` / `add_definitions` — use `target_*` commands only)
- Add **CMake presets** (`CMakePresets.json`) for clean multi-configuration builds (Debug/Release/RelWithDebInfo) and CI
- Use **FetchContent** or **CPM.cmake** for dependency management where possible
- Ensure `CMAKE_EXPORT_COMPILE_COMMANDS ON` is set — required for `clangd` IntelliSense to work correctly through deep template instantiations

#### Documentation: Doxygen + Awesome + Markdown narrative

Recommended approach (two layers):

**Layer 1 — API reference:** Doxygen with the [Doxygen Awesome](https://jothepro.github.io/doxygen-awesome-css/) CSS theme. Minimal additional setup, modern appearance, integrates with existing `///` comments.

**Layer 2 — Narrative documentation:** Markdown files (hosted via GitHub Pages or ReadTheDocs) covering conceptual explanations that Doxygen cannot express automatically:
- What `hamiltonian<RANK>` means physically (not just syntactically)
- The $H/i$ convention and non-Hermitian picture
- The continuous interaction-picture reset method
- How to define a new SubSystem and compose it into a simulation

Alternative: **Sphinx + Breathe + Exhale** (ingests Doxygen XML, produces better output, used by LLVM/OpenCV) — more setup but substantially better results if documentation becomes a priority for the paper.

#### Development environment

The framework builds cleanly under WSL2 (Ubuntu or Arch Linux) with the standard gcc/clang toolchain, CMake, and VS Code Remote-WSL + CMake Tools + clangd extensions.

### 7.9 Planned: SciPost Physics Codebases Paper 1 — MultiDiagonal Library

A standalone C++ library for multi-diagonal sparse operators over tensor-product Hilbert spaces, extracted from C++QED and published independently.

**Key contributions:**
- Algebraic closure under composition, direct product, Hermitian conjugation, and addition — with exact offset arithmetic
- Composition rule: $\text{offset}_\text{result} = \text{offset}_a + \text{offset}_b$ — the numbers of quanta moved by the two factors simply add up
- Direct product rule: concatenation of offset arrays reflects tensor product structure exactly
- Products of truncated ladder operators are computed exactly in the truncated space (the result equals the truncated dense product, not the textbook commutation relations, which fail at the cutoff)
- Application is $O(N \cdot d)$ where $d$ is the number of diagonals; composition is $O(d_a \cdot d_b)$ — far cheaper than general sparse-sparse multiply
- In quantum optics, the generators of the dynamics (ladder, number, Kerr, cross-Kerr, two-photon drive, Jaynes–Cummings coupling, $e^{inkx}$) are multi-diagonal in the Fock or momentum basis — the algebra is closed over this physically relevant set. The unitaries they generate (displacement, squeezing) are dense and **not** multi-diagonal
- `MultiDiagonal` is deliberately free of physics-layer concerns — no frequency information, no interaction-picture assumptions (see §4.3)

**Prerequisites before submission:**
- Signed-offset refactor (section 7.4)
- Extract as standalone CMake library independent of C++QED
- Unit tests: algebraic identities
- API documentation

### 7.10 Planned: SciPost Physics Codebases Paper 2 — C++QED v3

The main framework paper, citing Paper 1 for the `MultiDiagonal` component.

**Key contributions / narrative arc:**
1. **Motivation** — v2 achieved compile-time performance for large Hilbert spaces but OO/TMP complexity made it hard to use and extend; C++QED targets a fundamentally different use case from Python/QuTiP — production HPC/HTC batch execution, not interactive prototyping (see section 1.1)
2. **Core contribution** — C++20 concepts + ranges + Hana enables a functional rewrite that preserves the compile-time performance advantage while being far more readable and composable
3. **Key design ideas:** `retainedAxes` as NTTPs and the slice machinery; `LDO` laziness and TAG dispatch; reuse of state-vector operations for density matrix evolution via $H/i$ convention; lazy ensemble averaging chain; unified parameter/JSON machinery; continuous interaction-picture reset as a numerically critical method; first-class MDoF support with exact kinetic propagation; `run()` capabilities (checkpointing, autostopping, adaptive output density); parallel-safe PRNG (`pcg64` with independent stream parameters)
4. **Benchmarks** — performance comparisons against QuTiP for large Hilbert spaces
5. **Examples** — Mode, Qbit, Jaynes-Cummings, Particle-in-cavity demonstrating the user-facing API

**Prerequisites before submission:**
- Paper 1 published (or at least submitted)
- Fix open implementation gaps (§7.2–7.3)
- Particle / MDoF subsystem migrated (§6.5, §7.3)
- API documentation
- Public repository in clean state with build instructions
- Reproducible benchmark suite

---

## 8. Key Conventions and Idioms

| Idiom | Meaning |
|---|---|
| `dcomp` | Complex double — the primary scalar type throughout |
| `RANK / R` | Number of subsystems / tensor indices (compile-time) |
| `retainedAxes<i...>` | NTTP: compile-time set of axes to retain in slicing |
| `extendedAxes<rA,RANK>` | Concatenation of `retainedAxes` with itself shifted by `RANK` (for density operator diagonal slicing) |
| `Extents<R>` | `std::array<size_t,R>` — shapes and strides |
| `LogTree` | JSON-like tree (`boost::json`) for structured logging |
| `FWD(x)` | Perfect-forwarding macro: `std::forward<decltype(x)>(x)` |
| `hana_sequence<T>` | Concept: T is a Boost.Hana sequence (tuple-like) |
| `passByValue_v<T>` | Trait: true for views (`MultiArrayView`, LDO), false for owning types |
| `multiArrayRank_v<T>` | Trait: the underlying `MultiArray` rank of a quantum type |
| `compileTimeOrdinals<RANK>` | Generates `retainedAxes<0,1,...,RANK-1>` — all axes retained, used e.g. in `superoperatorFromJump` and Master row iteration |
| `slv<"label">(value)` | Labelled value wrapper carrying a string label into `temporal_data_point` for output |
| $H/i$ convention | Functionals apply $H/i$ (not $H$) to state vectors; $H$ is potentially non-Hermitian (includes no-jump backaction in QJMC/Master contexts) |
| `_("opt","desc",default,binding)` | Parameter registration tuple — used in `Pars` constructors to register with popl and JSON simultaneously |
| `noJSON` | Sentinel to exclude a parameter from JSON serialization (e.g. initial conditions) |
| MDoF | Motional degree of freedom — particle moving in a periodic potential; represented in momentum space with exact kinetic propagation |

---

## 9. File Reference

### 9.1 CPPQEDutils
| File | Description |
|---|---|
| `SliceIterator.h` | Compile-time axis slicing and partial trace machinery |
| `ODE.h` | Adaptive ODE engine (Boost.Odeint wrapper) |
| `Trajectory.h` | Generic trajectory run/stream/checkpoint driver |
| `StochasticTrajectory.h` | `stochastic` concept, `Ensemble`, `AverageTrajectories` customization point |
| `Pars.h` / `Pars.cc` | Parameter machinery: popl extension with `mod` namespacing and JSON serialization |

### 9.2 CPPQEDcore / quantumdata
| File | Description |
|---|---|
| `LazyDensityOperator.h` | Lazy pure/mixed state abstraction via TAG-based overloading |

### 9.3 CPPQEDcore / quantumoperator
| File | Description |
|---|---|
| `SparseMatrix.h` / `SparseMatrix.cc` | General sparse operator as a list of `(row, column, value)` triplets, applicable as a Hamiltonian |

`MultiArray.h`, `MultiArrayComplex.h` and `MultiDiagonal.h` (the strided array and the multidiagonal operator algebra, deliberately free of frequency/picture information) come from cppqed-multidiagonal as `<multidiagonal/…>`; `CPPQEDutils/FrameworkTraits.h` adds the framework-only traits on top of `<multidiagonal/Traits.h>`.

### 9.4 CPPQEDcore / quantumtrajectory
| File | Description |
|---|---|
| `QuantumTrajectory.h` | Shared infrastructure: `initialTimeStep`, `TDP_DensityOperator` |
| `Master.h` | Lindblad master equation on `DensityOperator`, reusing state-vector broadcast operations |
| `QuantumJumpMonteCarlo.h` | QJMC integrating and stepwise algorithms, `AverageTrajectories` specialization |

### 9.5 CPPQEDcore / quantumdynamics
| File | Description |
|---|---|
| `Hamiltonian.h` | Recursive compile-time Hamiltonian ($H/i$, non-Hermitian) composition via Hana sequences |
| `ExactPropagator.h` | Per-subsystem diagonal propagator (not necessarily unitary) |
| `Liouvillian.h` | Runtime vector of Lindblad operators with jump/rate/superoperator variants |
| `ExpectationValues.h` | Recursive Hana-sequence expectation value functionals over `LazyDensityOperator` |
| `Composite.h` | `Broadcaster<RANK,ra...>`: lifts rank-`\|ra\|` SubSystem descriptors to rank-`RANK` via slice iteration |

### 9.6 CPPQEDelements / subsystems
| File | Description |
|---|---|
| `Mode.h` / `Mode.cc` | Harmonic oscillator: `MultiDiagonal` operators, Kerr nonlinearity, coherent drive, photon loss/gain |
| `Qbit.h` / `Qbit.cc` | Two-level system: lambda Hamiltonians, `UnaryDiagonalPropagator` (implemented but blocked), Pauli operators as `MultiDiagonal` |
| `Particle_.h` / `Particle.cc` | v2 — Motional degree of freedom: `Spatial` grid, `Exact` propagator, `expINKX`/`cosNKX`/`sinNKX`, wavepacket and HO initial conditions; not yet migrated to v3 |
| `ParsParticle.h` / `ParsParticle.cc` | v2 — Parameter structs for Particle and PumpedParticle |

### 9.7 CPPQEDelements / interactions
| File | Description |
|---|---|
| `JaynesCummings.h` | Cavity-atom coupling via `Placeholder` slicing, compile-time level indices `<i,j>` |
| `ParticleCavity_.h` / `ParticleCavity.cc` | v2 — Particle-cavity interactions (`ParticleOrthogonalToCavity`, `ParticleAlongCavity`); not yet migrated to v3 |

### 9.8 examples
| File | Description |
|---|---|
| `jc_qjmc.cc` | Example highest-level executable (target `jc_qjmc`): QJMC simulation of qubit + mode with JC interaction |

---

## 10. Compile-Time Subsystem Dimensions — Design Discussion

*Added April 2026*

### 10.1 Motivation

The framework's central principle is that everything known at compile time should be processed at compile time. `RANK` and `retainedAxes` already satisfy this. Hilbert space dimensions, however, are currently entirely runtime values (`Extents<R>` = `std::array<size_t,R>`), even for systems whose dimension is structurally fixed: a qubit is always dimension 2, a qutrit always 3, a spin-j system always 2j+1.

The concrete gains from encoding fixed dimensions in the type:

- Loop bounds in `calculateSlicesOffsets` over static axes become compile-time constants
- `MultiDiagonal` dimension compatibility becomes a type error rather than a runtime check
- `Broadcaster` over a fully-static subsystem needs no stored extent state — it becomes a pure type
- `StateVector` for small fully-static systems can be stack-allocated with a statically known size
- Operator bodies in `Qbit`, `Mode` etc. gain static bounds checking

The primary use case is structurally fixed systems: qubits, qutrits, natural atoms, spins. Oscillator modes, whose cutoff is a convergence parameter swept across runs, remain runtime.

### 10.2 The Leading-Dimensions Design

Rather than a per-axis static/dynamic tag (which would require combinatorial case-handling in all slice machinery), a cleaner invariant is imposed: **only leading axes may have compile-time extents**.

```cpp
template <size_t RANK, size_t... StaticDims>
struct StateVector { ... };
// sizeof...(StaticDims) <= RANK always
// axes 0 .. sizeof...(StaticDims)-1 : compile-time extents
// axes sizeof...(StaticDims) .. RANK-1 : runtime extents
```

The split point is a single number `N_static = sizeof...(StaticDims)` — no per-axis tag, no sentinel values, no combinatorial case analysis anywhere in the implementation.

The extent accessor unifies both categories:

```cpp
template <size_t RANK, size_t... StaticDims>
struct ExtentsFor {
    static constexpr size_t N_static  = sizeof...(StaticDims);
    static constexpr size_t N_dynamic = RANK - N_static;

    static constexpr std::array<size_t, N_static>  static_extents{StaticDims...};
    std::array<size_t, N_dynamic> dynamic_extents{};   // runtime remainder

    size_t operator[](size_t i) const {
        if constexpr (N_static > 0)
            if (i < N_static) return static_extents[i];
        return dynamic_extents[i - N_static];
    }
};
```

**Backward compatibility:** `StaticDims...` empty recovers exactly `Extents<RANK>` — the current design. Existing code requires no changes; the new parameterization is a pure extension.

### 10.3 Fit with Physical Composition Order

The leading-dimensions invariant maps directly onto natural composition order. Discrete systems with fixed dimensions (qubits, atoms) occupy low axes; oscillator modes with variable cutoffs occupy high axes. In the Jaynes-Cummings example:

```
axis 0 = qubit  (dim 2, always)   ← static prefix
axis 1 = mode   (dim cutoff, ?)   ← dynamic suffix
```

`qbit::init` returns `StateVector<1, 2>`; `mode::init` returns `StateVector<1>`. The tensor product `operator*` concatenates static dim lists:

```
StateVector<1, 2>  ⊗  StateVector<1>  →  StateVector<2, 2>
```

The rule: **compose subsystems in order of decreasing dimension-knowability** — fixed-dimension subsystems first, variable-cutoff modes last. The type system then enforces what was previously only a convention.

### 10.4 The Transposition Constraint

This is the central limitation of the design, arising from the freedom in `retainedAxes`.

`retainedAxes<r0, r1, ..., rk>` can specify axes in any order — the order encodes transposition of the resulting slice. A transposition that moves a dynamic axis into a leading (static-prefix) position of the result, or a static axis out of it, would violate the invariant.

**The precise rule:** the static prefix length of the result equals the length of the longest leading subsequence of the retained axes list that falls entirely within `[0, N_static)`. As soon as a dynamic axis appears in the retained list, the static prefix of the result ends — even if further static axes appear later.

```cpp
constexpr size_t resultStaticPrefix(size_t N_static, /* retainedAxes pack */) {
    size_t count = 0;
    for (size_t r : retainedAxes)
        if (r < N_static) ++count; else break;
    return count;
}
```

Consequences for the actual `retainedAxes` uses in the codebase:

| Use | Retained axes | Result static prefix |
|---|---|---|
| `Broadcaster<2,0>` (qubit) | `{0}` | 1 — static dim preserved ✓ |
| `Broadcaster<2,1>` (mode) | `{1}` | 0 — fully dynamic |
| `compileTimeOrdinals<RANK>` | `{0,1,...,R-1}` in order | `N_static` — fully preserved ✓ |
| `extendedAxes<rA,RANK>` (density op diagonal) | static prefix doubled | preserved if source is ✓ |
| `transpose` (all axes, reordered) | arbitrary permutation | depends on permutation |

A crossed permutation is not a type error — the fallback is well-defined (fully dynamic result) — but the compiler can detect the crossing and emit the correct result type.

### 10.5 Scope of the Benefit

The leading-dimensions design is **not** a global compile-time guarantee that propagates freely through all slice machinery. It is local to:

- **SubSystem-level operations**: operator application, expectation values, and propagators computed directly on rank-1 (or low-rank) subsystem slices
- **Order-preserving slicing** (`retainedAxes` a strictly increasing subsequence of `{0,...,RANK-1}`) — the natural case for subsystem extraction
- **Leading subsystems in `operator*` composition** — the static prefix accumulates correctly as long as static-dimension subsystems are composed first

The benefit does **not** propagate through:
- Any `Broadcaster` over a dynamic-axis subsystem — that branch becomes fully dynamic immediately
- Any `retainedAxes` permutation that crosses the static/dynamic boundary
- The `transpose` special case, unless the permutation happens to preserve the leading structure

In a typical Jaynes-Cummings simulation, the qubit's `Broadcaster<2,0>` preserves the static prefix; the mode's `Broadcaster<2,1>` does not. The composite `StateVector<2,2>` retains the qubit's static extent in axis 0, but the slice machinery operating on the mode side works in the fully dynamic path regardless.

### 10.6 Relationship to `MultiDiagonal` and `UnaryDiagonalPropagator`

`MultiDiagonal<RANK>` stores diagonal extents as runtime values and validates them in `calculateAndCheckDimensions`. With compile-time dimensions, a `MultiDiagonal<RANK, StaticDims...>` specialization could validate compatibility statically — a natural companion to the `StateVector` extension.

`UnaryDiagonalPropagator` stores the diagonal as `std::valarray<dcomp>` sized at runtime. For `Qbit`, the diagonal has exactly 2 elements — known at compile time. A `std::array<dcomp, 2>` replacement would be correct and stack-allocated, contingent on the operator bodies being implemented first.

### 10.7 Implementation Priority and Sequencing

This is a non-trivial refactor touching `MultiArray.h`, `SliceIterator.h`, `StateVector`, `DensityOperator`, `LDO`, `Broadcaster`, and all SubSystem `make()` / `init()` functions. It should not be attempted until:

1. The open correctness gaps are closed (§7.1 bugs, §7.2 `UnaryDiagonalPropagator`, QJMC stepwise, autostop fix)
2. The test suite is in place (§7.6)
3. The `MultiDiagonal` signed-offset refactor (§7.4) is done

The natural sequencing milestone: coordinate with the acronym-casing decision (§7.5) into a single "type system" release; the `SubSystem` and `quantumdynamics/` renames are already done.

For now, the design intent should be recorded in header comments: `Qbit.h` should note that `Dim=2` is a candidate for NTTP promotion, and `SliceIterator.h` should note the transposition constraint that any future static-dimension scheme must respect.

### 10.8 Summary

| Property | Assessment |
|---|---|
| No per-axis combinatorics | ✓ Single split point `N_static` |
| Backward compatible | ✓ Empty `StaticDims...` = current design |
| Matches physical composition order | ✓ Fixed-dim subsystems on leading axes |
| `operator*` concatenation rule | ✓ Well-defined, enforces convention as type rule |
| Propagates through arbitrary `retainedAxes` | ✗ Transpositions crossing static/dynamic boundary drop to fully dynamic |
| Global compile-time guarantee | ✗ Local to SubSystem level and order-preserving slices |
| Benefit for primary use case (large oscillator) | Indirect — qubit ops gain; mode hot loop unchanged |
| Ready to implement | Not yet — correctness and test prerequisites first |

---

*— end of document —*
