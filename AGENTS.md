# AGENTS.md

Canonical instructions for AI coding agents (Claude Code, Codex, Copilot,
Cursor, Gemini CLI, ...) working on PositiveIntegrators.jl. Most tools read this
file natively; Claude Code reads it via the `@AGENTS.md` import in `CLAUDE.md`.
Keep this file the single source of truth; put only tool-specific notes into
tool-specific files.

When editing this file, only add instructions and commands that have actually
been tested in this repository. Instructions that do not work are worse than no
instructions.

PositiveIntegrators.jl is a Julia package providing positivity-preserving time
integration methods (modified Patankar-Runge-Kutta methods and related schemes)
for production-destruction systems (PDS), built on the interface of
OrdinaryDiffEq.jl/SciMLBase.jl. Correctness of the numerics matters more than
anything else: many bugs are silent (e.g., a wrong index or a missing diagonal
entry that only reduces accuracy or only affects one matrix type), so be careful
and verify mathematical changes by running code.

Before editing, inspect `git status` and preserve all existing work. Do not
discard, reset, or stash unrelated changes. The working tree may contain
untracked directories of the developer (e.g., `debug_*/` scratch projects);
never stage them, so do not use `git add -A` or `git add .`.

## Repository layout

- `ext/JuMPExt.jl`: package extension implementing `SanduProjection` (loaded
  with JuMP.jl).
- `test/runtests.jl`: the complete test suite in a single file of nested
  `@testset`s.
- `utils/`: developer scripts used below.
- `run/` (if present, ignored by git): the human developer's own Julia project.
  Do not modify it; agents use `run_agents/` instead (see below).

## Architecture

- `src/PositiveIntegrators.jl` contains all imports and exports and includes the
  other files. New public names go into the export lists there and into
  `docs/src/api_reference.md`, with a docstring including references to
  publications with DOI.
- `src/proddest.jl`: the problem types `PDSProblem` (production matrix `P` and
  destruction vector `D`) and `ConservativePDSProblem` (only `P`).
  - Convention: `P[i, j]` (for `i != j`) is the rate at which `u[j]` is
    converted into `u[i]`, so it is a production term of `u[i]` and a
    destruction term of `u[j]`. The diagonal `P[i, i]` contains production
    terms without destruction counterpart, and `D[i]` destruction terms without
    production counterpart. The right-hand side is
    `du[i] = P[i, i] + sum_j (P[i, j] - P[j, i]) - D[i]`.
  - In-place problems store `P` in `p_prototype`, which can be dense,
    `Tridiagonal`, or sparse (`SparseMatrixCSC`); out-of-place problems in the
    library use StaticArrays.jl.
  - The keyword argument `linear_invariants` (rows are the linear invariants)
    is experimental: it is only stored and used in tests, not read by the
    solvers.
- `src/mprk.jl` contains the building blocks shared by all schemes;
  `build_mprk_matrix` and `build_mprk_matrix!` have separate methods for dense,
  `Tridiagonal`, sparse, and static matrices, so a fix in one method usually
  needs to be applied to the others.
- Each scheme follows the OrdinaryDiffEqCore.jl pattern: an algorithm type, a
  constant cache for out-of-place and a mutable cache for in-place problems
  (`alg_cache`), `initialize!`, and `perform_step!`. `src/interpolation.jl`
  provides dense output for all caches listed in the union `MPRKCaches`; add
  new cache types there.
- Internal functions (e.g., `evaluate_pds`, `lincomb`, `basic_patankar_step`)
  are useful for investigations but not part of the API; access them via
  `PositiveIntegrators.name`.
- To read the source of a dependency, use the version that is actually loaded,
  e.g., `pkgdir(PositiveIntegrators.OrdinaryDiffEqCore)` (OrdinaryDiffEqCore.jl
  is not a direct test dependency, so `using OrdinaryDiffEqCore` fails in the
  `run_agents` project). Never edit installed package sources.

## Running Julia code

### The `run_agents` Julia projects

Use the local Julia projects in `run_agents/` at the root of this repository for
running tests, scripts, and the documentation. They are ignored by git and must
never be committed.

- `run_agents/` develops the local PositiveIntegrators.jl and contains all test
  dependencies.
- `run_agents/docs/` develops the local PositiveIntegrators.jl and contains all
  documentation dependencies.
- `run_agents/formatter/` contains JuliaFormatter.jl pinned to the version used
  in CI.

Create or update all three projects **before running tests, the docs, or the
formatter** by running, from any directory,

```bash
julia utils/setup_run_agents.jl
```

It is idempotent: it takes a few seconds if nothing changed and recreates a
project only if it is outdated (about 1.5 minutes from scratch, plus
precompilation on first use).

Do **not** use

- `julia --project=.` or `TestEnv.activate()`: they operate on the package
  environment and its untracked root `Manifest.toml`, which may be stale (e.g.,
  missing a newly added dependency), and `Pkg.resolve()` there takes minutes;
- `julia --project=docs` together with `Pkg.develop`: the tracked
  `docs/Project.toml` intentionally does not list PositiveIntegrators.jl, and
  `Pkg.develop` adds it there;
- `julia --project=test`: the test environment lacks PositiveIntegrators.jl;
- stacked environments via `JULIA_LOAD_PATH` or `push!(LOAD_PATH, ...)`.

Add packages needed only for your own experiments (e.g., BenchmarkTools.jl) to
`run_agents/` (they are removed when the setup script recreates the project),
never to the tracked `Project.toml` files unless this is part of the task.
Keep temporary scripts and output outside the repository or in `run_agents/`.

### Persistent Julia sessions via MCP

Some developers configure MCP servers providing persistent Julia sessions, in
which packages are loaded and compiled only once. **If such tools are available
to you** (possibly only listed by name until they are loaded, see `CLAUDE.md`
for Claude Code), use them instead of starting `julia` in the shell for tests,
experiments, and the docs. Run the setup script above in the shell in any case,
even if the description of an MCP tool tells you to never run Julia via the
command line.

For [julia-mcp](https://github.com/aplavin/julia-mcp) (tools `julia_eval`,
`julia_restart`, and `julia_list_sessions`):

- Pass the **absolute path** of `run_agents` (or `run_agents/docs`,
  `run_agents/formatter`) as `env_path`. The working directory of the session
  is this directory, so use absolute paths.
- Pass a generous `timeout` (e.g., `1800`) for the first call in a session,
  which loads and compiles packages.
- Call `julia_restart` with the same `env_path` after re-creating `run_agents`.
- Do not override package methods with `@eval` in a session to try out a
  change; edit the source file instead (and revert it afterward). Such an
  override has crashed a session with a segmentation fault.
- Sessions may run with non-default Julia flags chosen by the developer. Check
  them with `Threads.nthreads()` and `Base.JLOptions().check_bounds` (`1` means
  bounds checking is enabled) if they matter.

## Testing

### Running tests

The complete test suite (`test/runtests.jl`, more than 120 000 tests) takes
about 25 minutes of CPU time locally and about one hour in CI (with coverage).
CI runs it on Linux, macOS, and Windows with Julia 1.10 and 1.13, plus a
`Downgrade` job with the lowest compatible versions of all dependencies.
**Do not run the complete test suite while developing**; run the relevant
testsets instead, and run the complete suite at most once at the end, and only
if the user agrees.

Run selected testsets by name with `utils/run_testsets.jl`, e.g., in a julia-mcp
session with `env_path` set to the absolute path of `run_agents`:

```julia
include("/path/to/PositiveIntegrators.jl/utils/run_testsets.jl")
run_testsets("Float32 input data", "HIRES reference solution");
```

or in the shell, from the root of the repository,

```bash
julia --project=run_agents -e '
include("utils/run_testsets.jl")
run_testsets("Float32 input data", "HIRES reference solution")'
```

- Copy names exactly from the `@testset "..."` lines (e.g., with
  `grep -n '@testset' test/runtests.jl`); they are matched wherever they are
  nested. Testsets with interpolated names (e.g., `"$alg"`) cannot be
  selected; select the enclosing testset. The same applies if a nested testset
  uses variables defined in its parent.
- Run the complete suite with
  `julia --project=run_agents test/runtests.jl` from the root of the
  repository (or `include` it in a julia-mcp session).

### Adding and changing tests

- New or changed behavior must be covered by focused tests; the overall test
  coverage must not decrease. For a bug fix, verify that the test detects the
  original behavior, e.g., by temporarily reverting only your own change.
- Add tests to the matching `@testset` in `test/runtests.jl` (search for
  neighboring tests of the same feature) and keep each new testset fast.
- Many bugs affect only some problem representations. For a new or changed
  algorithm, check (and extend the corresponding testsets for) in-place and
  out-of-place problems, conservative and nonconservative problems, dense,
  `Tridiagonal`, sparse, and static matrices, `Float32` and `Double64` data,
  and Unitful.jl units. Use enough time steps that multistage or start-up
  phases are actually passed (e.g., not only 4 steps).
- Compare with an analytic solution (`analytic` keyword) or a reference
  solution computed with a high-order standard solver at tight tolerances
  (e.g., `Rodas4P()` or `Vern9()` with `abstol = 1e-14, reltol = 1e-13`).
  Never loosen tolerances just to make a failing test pass; if a tolerance
  must change, understand and explain why.

### Adding a test problem

Adding a problem `prob_pds_xxx` to `src/PDSProblemLibrary.jl` touches:

1. `src/PDSProblemLibrary.jl`: `P`, `D` (if not conservative), a standard
   right-hand side `std_rhs`, `u0`, the problem with `linear_invariants` if
   available, and a docstring with the equations, the time span, and
   references (with DOI).
2. The export list in `src/PositiveIntegrators.jl` and the `@docs` block of
   example problems in `docs/src/api_reference.md`.
3. `NEWS.md` (section "Added" of the current lifecycle).
4. `test/runtests.jl`: the testsets
   `"Check that production-destruction form and standard ODE fit together in predefinded problems"`,
   `"PDS problem library (adaptive schemes)"`,
   `"PDS problem library (non-adaptive schemes)"`, and
   `"Check that the predefined linear invariant matrices fit in predefined problems"`,
   and a new testset comparing with a reference solution (e.g.,
   `"HIRES reference solution"`).

Transcribe the problem from primary sources (papers, books with page numbers,
or reference implementations such as the Fortran sources of the
[Test Set for IVP Solvers](https://archimede.uniba.it/~testset/)) and verify
it numerically: all entries of `P` and `D` are nonnegative for nonnegative
states, `std_rhs` agrees with the PDS form, the linear invariants are
preserved, and the solution agrees with published reference values.

## Documentation

- The complete documentation build runs all tutorials and benchmarks and takes
  about 40 minutes in CI. **Do not build the complete documentation**; build
  only the pages you changed with `utils/build_docs_pages.jl`, e.g., in a
  julia-mcp session with `env_path` set to the absolute path of
  `run_agents/docs`,
  ```julia
  include("/path/to/PositiveIntegrators.jl/utils/build_docs_pages.jl")
  build_docs_pages("scalar_pds.md")
  ```
  or in the shell, from the root of the repository,
  ```bash
  julia --project=run_agents/docs utils/build_docs_pages.jl scalar_pds.md
  ```
  The HTML pages and SVG figures are written to `run_agents/docs/build/`.
  Warnings about the missing `index.md` and about cross-references to other
  pages are expected there.
- In CI, Documenter.jl runs without `warnonly`, so broken cross-references
  (``[`name`](@ref)``), missing docstrings in `docs/src/api_reference.md`, and
  failing `@example` blocks are errors.
- Register new pages in the `pages` list of `docs/make.jl`.
- For a new benchmark page, copy the structure of an existing one (e.g.,
  `docs/src/hires_benchmark.md`). Standard solvers that may produce negative
  values are run with `isoutofdomain = isnegative`.
- Keep the generated figures small (SVG output of every docs build is stored in
  the `gh-pages` branch). The size is dominated by the number of plotted points
  and markers, not by the format: plot long reference solutions with
  `plotdensity = 1000` or a coarse `plotat`, and draw markers only if there are
  few steps (e.g., `markers = length(sol.t) <= 500 ? :circle : :none`). Check
  the size of the generated SVGs in `run_agents/docs/build/`.
- Plots.jl silently ignores unknown or misspelled keyword arguments (e.g.,
  `linewidht`, or `markevery`, which does not exist). Verify every plotting
  change by rendering the figure and checking that it changed as intended.

## Code style

- Format all changed Julia files with JuliaFormatter.jl v1.0.60 (configuration
  in `.JuliaFormatter.toml`, SciML style); CI checks the formatting of the
  whole repository. Other versions format differently, so use the pinned
  version in `run_agents/formatter/`, e.g., from the root of the repository,
  ```bash
  julia --project=run_agents/formatter -e 'using JuliaFormatter; println(format(["src", "test", "docs", "ext", "utils"]))'
  ```
  `format` returns `true` if all files were already formatted; otherwise, it
  formats them in place and returns `false`. After formatting, re-check
  numerically sensitive code you changed (the formatter rewrites, e.g.,
  the separators of literal `@SVector [...]`/`@SMatrix [...]` arrays).
- CI checks spelling with [typos](https://github.com/crate-ci/typos); if it is
  installed, run `typos` at the root of the repository.
- Use `@muladd` on numerical kernels as in the existing schemes.
- Keep code type-stable and generic in the number type (`Float32`, `Float64`,
  `Double64`, and Unitful.jl quantities are tested) and in the matrix type of
  `P`. Do not fix the matrix type of fields of mutable caches to `Matrix`;
  assigning a `Tridiagonal` or sparse matrix to such a field silently converts
  it to a dense matrix.
- Do not use `@inbounds` in problem definitions: after
  `remake(prob; u0 = ...)` with a different length, out-of-bounds accesses are
  silent. Elsewhere, use it only where profiling shows a benefit and the
  indices are provably in bounds.

## Investigations and reviews

- Verify claims by running code in the `run_agents` project instead of only
  reading code, and report which statements were verified and which are
  inferences. Give `file:line` locations.
- For a suspected bug, check whether it affects the shipped problems and
  default settings or only unusual input, and report this distinction.
- Behavior that looks like a bug can be inherent to the method. For example,
  MPRK schemes need many more steps than Rosenbrock methods for HIRES
  (`prob_pds_hires`) because the Patankar weights treat the stiff nonlinear
  reaction term inaccurately, not because of the step size controller.

## Pull requests

- Document user-visible changes in `NEWS.md`.
- Do not commit generated files (`docs/build/`, figures, `Manifest.toml`
  files, `run_agents/`) or temporary notes.
- Do not change the package version in `Project.toml` unless asked.
- Do not commit or push unless asked to do so by the user.
- `git fetch` is slow if it fetches the `gh-pages` branch (hundreds of MB of
  generated documentation).
