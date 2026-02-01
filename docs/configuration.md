# Configuration (CRYSTAL23)

The workflow is designed to run with **CRYSTAL23** while keeping backwards
compatibility with older CRYSTAL14-era environment variables.

## Recommended setup

Source CRYSTAL23's environment in your shell (example):

```bash
source /path/to/CRYSTAL23/cry23.bashrc
```

This typically defines:

- `CRY23_EXEDIR`  (directory containing `crystal`, `Pcrystal`, etc.)
- `CRY23_UTILS`   (directory containing `runcry23`, `runPcry23`, `runprop23`, etc.)
- `CRY23_SCRDIR`  (scratch directory)

## Environment variables used by this repo

### CRYSTAL executables
- Primary: `CRY23_EXEDIR`
- Legacy aliases recognised:
  - `CRYSTAL_EXEDIR`
  - `CRY14_EXEDIR`

### Scratch directory
- Primary: `CRY23_SCRDIR`
- Legacy aliases recognised:
  - `CRYSTAL_TMPDIR`
  - `CRYSTAL_SCRDIR`
  - `CRY14_SCRDIR`

### Utilities (properties)
- Primary: `CRY23_UTILS` (optional)
- If not set, `runprop23` must be available on your `PATH`.

## Helper binaries

`billy` optionally uses auxiliary helpers:
- `helpbilly`, `helpbilly2` (compiled from `fortran/helpbilly*.f90`)
- `dfit` (compiled from `fortran/dfit/*.f90`)

Install these binaries into your `PATH` (recommended) and `billy` will find them automatically.




## One-time setup (line endings + executable bits)

If you copied this repo through Windows tooling (or downloaded as a ZIP), some scripts may
contain CRLF line endings (`\r\n`). Before first use on Linux/WSL/HPC, normalise line endings
*first*, then mark scripts as executable:

```bash
# from repo root
# 1) normalise CRLF -> LF for all executable scripts under ./scripts
sed -i 's/\r$//' scripts/billy scripts/runcrystal scripts/properties
find scripts -type f \( -name "*.sh" -o -name "*.csh" -o -name "*.py" \) -exec sed -i 's/\r$//' {} +

# 2) make scripts executable
chmod u+x scripts/billy scripts/runcrystal scripts/properties
find scripts -type f \( -name "*.sh" -o -name "*.csh" -o -name "*.py" \) -exec chmod u+x {} +

# 3) make helper executables executable (if you compiled / copied them here)
chmod u+x fortran/helpbilly fortran/helpbilly2 fortran/dfit/dfit
```

