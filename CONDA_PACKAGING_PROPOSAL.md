# Proposal: Making SpecHLA Conda/Bioconda-Packageable

## Motivation

SpecHLA is a powerful HLA typing tool, but installing it requires cloning the repo, setting up a conda environment manually, running `index.sh`, and relying on vendored pre-compiled binaries. A bioconda package would let users install with a single `conda install spechla` command. However, several architectural changes are needed first.

This document outlines the required changes to the SpecHLA codebase to make it compatible with bioconda packaging.

---

## Current Blockers

| Issue | Why it blocks conda packaging |
|-------|-------------------------------|
| **Hardcoded relative paths** | Scripts locate binaries via `$dir/../../bin` and databases via `$dir/../../db`. In a conda install, binaries are on `$PATH` and data lives in `$PREFIX/share/`. |
| **Vendored pre-compiled binaries** | `bin/bcftools`, `bin/inchworm`, `bin/novoalign`, `bin/novoindex` are checked-in ELF binaries. Conda builds everything from source; pre-compiled binaries are not allowed. |
| **No install targets** | `ExtractHAIRs` CMakeLists.txt has no `install()` directive. The project has no mechanism to install scripts to a system location. |
| **Debug build flags** | Both CMakeLists.txt default to `CMAKE_BUILD_TYPE=DEBUG` with `-O0` and hardcode `g++` as the compiler. Conda requires Release builds and uses its own compiler wrappers. |
| **Fermikit bundled with full source** | The `bin/fermikit/` directory contains a full copy of fermikit source + pre-built binaries. The `fermikit` bioconda package already provides all needed tools. |
| **`index.sh` generates absolute-path config files** | The HLA config files (`db/HLA/HLA_*.config.txt`) contain hardcoded absolute paths set at install time. |

---

## Proposed Changes

### 1. Centralized Path Resolution

Create two small helper modules that all scripts source/import. They resolve paths using a priority chain:

1. **Environment variable** (`SPECHLA_DB`) — allows user override
2. **Conda location** (`$CONDA_PREFIX/share/spechla/db`) — for conda installs
3. **Relative path** (`../../db`) — backward compatible for development installs

**New file: `script/spechla_env.sh`** (sourced by all Bash scripts)
```bash
#!/bin/bash
# Resolve SpecHLA paths for conda and development installs
if [ -n "${SPECHLA_DB:-}" ] && [ -d "$SPECHLA_DB" ]; then
    SPECHLA_DB="$SPECHLA_DB"
elif [ -n "${CONDA_PREFIX:-}" ] && [ -d "$CONDA_PREFIX/share/spechla/db" ]; then
    SPECHLA_DB="$CONDA_PREFIX/share/spechla/db"
else
    _spechla_root=$(cd "$(dirname "$(realpath "${BASH_SOURCE[0]}")")/.." && pwd)
    SPECHLA_DB="$_spechla_root/db"
fi
export SPECHLA_DB
```

**New file: `script/spechla_paths.py`** (imported by all Python scripts)
```python
import os

def get_db_dir():
    if os.environ.get('SPECHLA_DB') and os.path.isdir(os.environ['SPECHLA_DB']):
        return os.environ['SPECHLA_DB']
    conda = os.environ.get('CONDA_PREFIX', '')
    conda_db = os.path.join(conda, 'share', 'spechla', 'db')
    if conda and os.path.isdir(conda_db):
        return conda_db
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'db')
```

### 2. Replace Vendored Binary Paths with PATH Lookups

Since conda puts all tool binaries on `$PATH`, every reference to a vendored binary should become a simple command name.

**Before:**
```bash
bin=$dir/../../bin
$bin/bcftools filter -R ...
$bin/novoalign -d $db/ref/...
```

**After:**
```bash
source "$(dirname $(realpath $0))/../spechla_env.sh"
bcftools filter -R ...
novoalign -d $SPECHLA_DB/ref/...   # still conditional on `which novoalign`
```

**Before (Python):**
```python
order = '%s/../bin/SpecHap/build/SpecHap --ncs ...' % sys.path[0]
order = '%s/../bin/bcftools consensus ...' % sys.path[0]
```

**After (Python):**
```python
order = 'SpecHap --ncs ...'
order = 'bcftools consensus ...'
```

**Before (Perl):**
```perl
my $db = "$Bin/../../db/HLA";
```

**After (Perl):**
```perl
my $db = $ENV{SPECHLA_DB} ? "$ENV{SPECHLA_DB}/HLA" : "$Bin/../../db/HLA";
```

#### Files requiring path changes

**Bash scripts:**
- `script/whole/SpecHLA.sh` — `$bin/bcftools`, `$bin/novoalign`, `$db` default
- `script/whole/test_SpecHLA.sh` — same pattern
- `script/run.assembly.realign.sh` — `$sdir/fermikit/fermi.kit/fermi2.pl` -> `fermi2.pl`; `$sdir/blast2sam.pl` -> `blast2sam.pl`
- `script/ExtractHLAread.sh` — `$db` default
- `script/ScanIndel/run_scanindel_sample.sh` — bin/db references

**Python scripts:**
- `script/phase_variants.py` — **~25 path references** (largest file): `SpecHap/build/SpecHap`, `extractHairs/build/ExtractHAIRs`, `bcftools`, db paths
- `script/long_read_typing.py` — `parameter.bin`, `sys.path[0]/../bin` patterns
- `script/typing_from_assembly.py` — `sys.path[0]/../db/`
- `script/refine_typing.py` — `--db` default
- `script/mask_low_depth_region.py` — `sys.path[0]/whole/exon_extent.bed`
- `script/whole/map_block2_database.py` — `$bin/bcftools`
- `script/whole/g_group_annotation.py` — db path
- `script/whole/top_allele_2_reads.py` — `sys.path[0]` db references

**Perl scripts:**
- `script/whole/annoHLA.pl` — `$Bin/../../db/HLA`
- `script/whole/select.combination.pl` — same
- `script/count.read.pl` — `$Bin/../db/ref`
- `script/cal.hla.copy.pl` — db references

### 3. Fix CMake Build System

**`bin/SpecHap/CMakeLists.txt`:**
```diff
- set(CMAKE_CXX_COMPILER g++)
- set(CMAKE_C_COMPILER gcc)
+ # Let the build system (conda or user) choose the compiler

- add_definitions(-D_GLIBCXX_USE_CXX11_ABI=0)
+ # Remove ABI override — causes linking issues with conda-provided libraries

- set(CMAKE_BUILD_TYPE DEBUG)
- set(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXX_FLAGS} -O0 -Wall -ggdb -fkeep-inline-functions")
+ if(NOT CMAKE_BUILD_TYPE)
+   set(CMAKE_BUILD_TYPE Release)
+ endif()

- if (NOT $ENV{CONDA_PREFIX} STREQUAL "")
-     include_directories($ENV{CONDA_PREFIX}/include)
-     link_directories($ENV{CONDA_PREFIX}/lib)
- endif ()
+ # Use CMAKE_PREFIX_PATH instead (passed by the build system)
```

**`bin/extractHairs/CMakeLists.txt`:**
Same fixes as above, plus add the missing install target:
```cmake
install(TARGETS ExtractHAIRs DESTINATION bin)
```

### 4. Remove Vendored Binaries

| File | Action | Replacement |
|------|--------|-------------|
| `bin/bcftools` | Delete | `bcftools` conda dependency |
| `bin/novoalign` | Delete | Users install separately; found on PATH |
| `bin/novoindex` | Delete | Users install separately; found on PATH |
| `bin/inchworm` | Delete | `trinity` conda dependency (provides inchworm) |
| `bin/fermikit/` | Delete entire directory | `fermikit` conda dependency |
| `bin/blast2sam.pl` | Keep, install to `$PREFIX/bin/` | Installed by build.sh |
| `bin/vcf-combine.py` | Keep, install to `$PREFIX/bin/` | Installed by build.sh |
| `bin/SpecHap/` | Keep source | Built from source by conda |
| `bin/extractHairs/` | Keep source | Built from source by conda |

### 5. Update `index.sh`

The current `index.sh` performs four tasks. In a conda install, most are handled by the build recipe:

| Task | Conda handling | Dev-mode handling |
|------|---------------|-------------------|
| Generate HLA config files | Done during `conda build` | Keep in `index.sh` |
| Build bowtie2 indexes | Done during `conda build` | Keep in `index.sh` |
| Symlink libncurses | Not needed (conda handles deps) | Remove |
| Build SpecHap/ExtractHAIRs | Done by CMake in build recipe | Keep in `index.sh` |

Update `index.sh` to use `spechla_env.sh` for path resolution, and change novoalign detection from checking a vendored binary to checking if `novoalign` is on PATH.

### 6. Database Handling

The `db/` directory (163MB) will be **bundled inside the conda package** at `$PREFIX/share/spechla/db/`. This is standard practice for bioconda tools with reference databases of this size.

The database includes:
- HLA allele sequences and frequency tables (`db/HLA/`, 106MB)
- Reference FASTA files with pre-built BWA/BLAST indexes (`db/ref/`, 57MB)
- Bowtie2 indexes are built at conda install time (not pre-existing in the repo)

### 7. Entry Point Wrappers

Conda will install wrapper scripts in `$PREFIX/bin/` that set environment variables and delegate to the real scripts:

| Command | Underlying script |
|---------|-------------------|
| `spechla` | `script/whole/SpecHLA.sh` |
| `spechla-extract-hla-reads` | `script/ExtractHLAread.sh` |
| `spechla-long-read` | `script/long_read_typing.py` |
| `spechla-assembly` | `script/typing_from_assembly.py` |
| `spechla-loh` | `script/cal.hla.copy.pl` |

Example wrapper:
```bash
#!/bin/bash
set -euo pipefail
export SPECHLA_DB="${SPECHLA_DB:-${CONDA_PREFIX}/share/spechla/db}"
exec bash "${CONDA_PREFIX}/share/spechla/script/whole/SpecHLA.sh" "$@"
```

---

## Bioconda Recipe (Overview)

### `meta.yaml` key sections

```yaml
requirements:
  build:
    - {{ compiler('cxx') }}
    - cmake >=3.6
    - make
  host:
    - htslib
    - arpack
    - zlib
  run:
    - python >=3.8
    - samtools
    - bcftools
    - bwa
    - bowtie2
    - freebayes
    - minimap2
    - blast
    - blat
    - bedtools
    - bamutil
    - fermikit
    - k8
    - parallel
    - trinity
    - pbmm2
    - pbsv
    - longshot
    - pysam
    - biopython
    - numpy
    - pandas
    - scipy
    - pyvcf
    - pulp
    - python-edlib
    - perl
```

### `build.sh` steps

1. Build SpecHap and ExtractHAIRs via CMake with `CMAKE_INSTALL_PREFIX=$PREFIX`
2. Copy `script/` to `$PREFIX/share/spechla/script/`
3. Copy `db/` to `$PREFIX/share/spechla/db/`
4. Install utility scripts to `$PREFIX/bin/`
5. Generate HLA config files with correct paths
6. Build bowtie2 indexes
7. Install wrapper scripts to `$PREFIX/bin/`

---

## Dependency Version Notes

The current `environment.yml` pins very old versions (e.g., `samtools=1.3.1`, `bedtools=2.26.0`). For most tools, the bioconda recipe should let the solver resolve compatible versions. However:

- **samtools >=1.10** is required because ExtractHAIRs uses htslib APIs (`hts_itr_multi_next`, `sam_hdr_destroy`) that were introduced in htslib 1.10. The old `samtools depth -d` flag was removed in samtools 1.13, but the `-d 1000000` option was only used to set an unlimited depth cap — which is the default in newer versions, so removing it is safe.
- For other tools, no exact version pins unless incompatibilities are found
- If incompatibilities are discovered, add minimum version constraints (e.g., `bedtools >=2.26`)

---

## Implementation Order

1. **Path resolution helpers** — `spechla_env.sh` + `spechla_paths.py` (non-breaking)
2. **CMake fixes** — build type, compiler, install target (non-breaking)
3. **Script path refactoring** — update all ~20 scripts (backward-compatible with fallback)
4. **Remove vendored binaries** — delete pre-compiled files from `bin/`
5. **Update `index.sh`** — use new path resolution
6. **Create bioconda recipe** — `meta.yaml` + `build.sh`
7. **Update tests** — adapt `example/test_*.sh` for conda installs
8. **Update README** — add conda installation instructions
9. **Tag release** — v1.0.9

---

## Backward Compatibility

All changes maintain backward compatibility for development installs:
- The path resolution chain falls back to relative paths when `SPECHLA_DB` and `CONDA_PREFIX` are not set
- `index.sh` continues to work for manual installations
- Existing users who clone the repo and follow current instructions will not be affected

---

## Risks

| Risk | Mitigation |
|------|------------|
| Newer samtools/bcftools break something | Test all example workflows with current versions before release |
| `trinity` package doesn't provide `inchworm` on PATH | Verify in a test conda env; if not, build inchworm separately |
| Config file absolute paths break on env relocation | Long-term: make ScanIndel resolve paths dynamically |
| `phase_variants.py` has deeply embedded shell commands with path interpolation | Test every code path (PE, TGS, Hi-C, 10X, Nanopore) after refactoring |
