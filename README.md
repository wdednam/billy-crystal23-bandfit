# billy-extended (CRYSTAL23)

This repository contains a modified version of the **billy** basis-set / geometry
optimisation script originally written by Mike Towler for CRYSTAL. The original author
granted explicit permission for unrestricted use, modification, and redistribution
(see `NOTICE`).

This fork focuses on:
- **CRYSTAL23-first** execution (while supporting legacy CRYSTAL14-era env vars).
- A cleaner, reproducible layout suitable for an academic software release.
- Additional band-structure fitting helpers used in the associated publication.

## Quick start

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



1) Ensure CRYSTAL23 is configured in your environment:

```bash
source /path/to/CRYSTAL23/cry23.bashrc
```

2) Put the repo scripts on your `PATH` (example):

```bash
export PATH="$PWD/scripts:$PATH"
```

3) Build helper binaries and add them to your `PATH`:

```bash
# Example (adapt to your compiler / environment)
cd fortran
gfortran -O2 -o helpbilly helpbilly.f90
gfortran -O2 -o helpbilly2 helpbilly2.f90
cd dfit
make
export PATH="$PWD:$PATH"
```

4) Run `billy` on an example:

```bash
cd examples/Al_POB_DZVP_rev2
billy -f -n 3 -np 1 -cost ref.band.npy -wE 0.15 -lambda 0.01 -cap 0.55 -erange -10.0 10.0 Al 15
```

## Configuration

See `docs/configuration.md`.

## Repository contents

- `scripts/` : executables and wrappers (`billy`, `runcrystal`, `properties`, helpers)
- `examples/`: worked examples and sweep scripts
- `fortran/` : helper utilities (optional but recommended)
- `upstream/`: original upstream README files (for reference)
