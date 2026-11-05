#!/usr/bin/env julia

# Create or update the local Julia projects used by AI coding agents (see `AGENTS.md`):
# - `run_agents/` develops the local version of PositiveIntegrators.jl and contains all
#   test dependencies (copied from `test/Project.toml`, including compat bounds and
#   preferences).
# - `run_agents/docs/` develops the local version of PositiveIntegrators.jl and contains
#   all dependencies of the documentation (copied from `docs/Project.toml`). The tracked
#   `docs/Project.toml` does not list PositiveIntegrators.jl and must not be modified.
# - `run_agents/formatter/` contains JuliaFormatter.jl pinned to the version used in CI.
# This script is idempotent: it only re-creates a project if it is outdated. Run it as
#   julia utils/setup_run_agents.jl
# from any directory; all paths are relative to the root of the repository.

using Pkg, TOML

const REPO_DIR = dirname(@__DIR__)
const RUN_AGENTS_DIR = joinpath(REPO_DIR, "run_agents")
const FORMATTER_DIR = joinpath(RUN_AGENTS_DIR, "formatter")
const FORMATTER_VERSION = "1.0.60"

parse_toml(file) = isfile(file) ? TOML.parsefile(file) : Dict{String, Any}()

# Check whether all entries of `required` are contained in `actual`, recursively for
# nested tables such as `[preferences.HostCPUFeatures]`. Additional entries are allowed,
# e.g., packages or compat bounds added by developers.
contains_entries(actual, required) = actual == required
function contains_entries(actual::AbstractDict, required::AbstractDict)
    return all(haskey(actual, key) && contains_entries(actual[key], value)
               for (key, value) in required)
end

# Create or update a project in `dir` with all entries of `source_project_file` that
# develops this checkout of PositiveIntegrators.jl
function setup_project(dir, source_project_file)
    mkpath(dir)
    required = TOML.parsefile(source_project_file)
    project_file = joinpath(dir, "Project.toml")
    project = parse_toml(project_file)
    manifest = parse_toml(joinpath(dir, "Manifest.toml"))

    deps = get(project, "deps", Dict{String, Any}())
    entries = get(get(manifest, "deps", Dict{String, Any}()), "PositiveIntegrators", [])
    sections_match = all(contains_entries(get(project, section, Dict{String, Any}()),
                                          get(required, section, Dict{String, Any}()))
                         for section in ("deps", "compat", "extras", "preferences"))
    develops_local = haskey(deps, "PositiveIntegrators") &&
                     any(haskey(entry, "path") &&
                         samefile(joinpath(dir, entry["path"]), REPO_DIR)
                         for entry in entries)

    Pkg.activate(dir)
    if !(sections_match && develops_local)
        rm(joinpath(dir, "Manifest.toml"); force = true)
        cp(source_project_file, project_file; force = true)
        # A relative path is stored relative to the project in the manifest
        Pkg.develop(path = relpath(REPO_DIR))
    end
    Pkg.instantiate()
end

function setup_formatter()
    manifest = parse_toml(joinpath(FORMATTER_DIR, "Manifest.toml"))
    entries = get(get(manifest, "deps", Dict{String, Any}()), "JuliaFormatter", [])
    pinned = any(get(entry, "version", nothing) == FORMATTER_VERSION &&
                 get(entry, "pinned", false) for entry in entries)

    pinned || rm(FORMATTER_DIR; force = true, recursive = true)
    Pkg.activate(FORMATTER_DIR)
    if !pinned
        Pkg.add(name = "JuliaFormatter", version = FORMATTER_VERSION)
        Pkg.pin("JuliaFormatter")
    end
    Pkg.instantiate()
end

setup_project(RUN_AGENTS_DIR, joinpath(REPO_DIR, "test", "Project.toml"))
setup_project(joinpath(RUN_AGENTS_DIR, "docs"), joinpath(REPO_DIR, "docs", "Project.toml"))
setup_formatter()
