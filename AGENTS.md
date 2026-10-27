# AGENTS.md

## Project Overview
A nextflow pipeline to process and analyze plant transcription factors.

### Environments
The main project environment is managed by pixi. There are subprojects:
- genesis_tools is a Python project managed by uv. It has software tools needed for some Processes in the main pipeline.
- environments/GENESIS_ALIGNMENT.yaml defines a conda environment used for alignment and indexing.
- environments/GENESIS_MACS3.yaml defines a conda environment used for peak calling and intersection.
- environments/GENESIS_SPP.yaml defines a conda environment used for QC.

### Compute resources
The project runs locally on macOS or Linux and on supported scheduler profiles.
- Use existing SSH authentication. If authentication fails, report the issue;
  do not create or modify credentials, keys, or SSH configuration.
- Use only accounts, queues and allocations authorized for the task. Keep local
  account names and filesystem configuration outside the repository.
- Access shared project storage only within the explicitly authorized scope;
  modifications require authorization. Do not infer access from mounted paths.

## Coding Standards
- Follow Python best practices and use type annotations.
- Follow bash script best practices.
- Prefer type narrowing or runtime validation. Use `typing.cast` when the type is already
  established but the checker cannot infer it. Avoid blanket ignores.

## Verification
- Add focused tests for changed behavior when practical, using existing test infrastructure.
  Ask before introducing a new test framework or tests requiring external services or substantial
  resources.
- Run checks for the affected component during development. Run the full suite before completing
  code changes when feasible. Reviews and documentation-only changes do not require the full suite.
- If the full suite is not run, state why and list the checks that were run.
- Some check tasks run `ruff check --fix` and `ruff format`, which modify files. Inspect the
  resulting diff and preserve unrelated user changes.
- Introduce no new failures or warnings. Fix those caused by the change. Report unrelated existing
  failures separately, and do not describe incomplete or failing checks as passing.
- Do not weaken assertions or suppress diagnostics merely to pass checks. Update test expectations
  when the requested behavior changes, while preserving meaningful coverage.

## Boundaries
Keep assistant settings, environment files and credentials local and untracked.
The public `conf/dap-images.lock` contains container image pins only; never put
credentials or machine-specific settings in it.
- Allowed without prompting:
  - Creation of temporary files as long as they are cleaned up.
  - Read-only bash commands that will not alter or create files, and will not read non-temporary files outside the project folder.
  - Linting, type checking, formatting, or running python tests.
- Agent-submitted Sherlock jobs must request at most five minutes of walltime, two CPUs, and 8 GB
  of total memory. Set explicit limits; do not launch jobs requiring more resources under these
  instructions.
- Do not run resource-intensive jobs on the login node.
- Avoid adding new dependencies to environments.
- Ask before pushing docker images.

Keep this root instruction file compact. Put detailed component-specific guidance in nested
`AGENTS.md` files when needed.
