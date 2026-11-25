# Licenses and notices of the shipped CleanSeams artifacts

CleanSeams is licensed **GPL-3.0-or-later**. The license text is `COPYING` (repository and
artifact root), byte-identical to the COPYING of the CleanSeams 0.1.0 add-on
(sha256 `230184f60bae2feaf244f10a8bac053c8ff33a183bcc365b4d8b876d2b7f4809`); a copy is in
`cleanseams/COPYING`.

Every runtime artifact (the Blender extension archive, the `cleanseams_native` wheel and the
CLI archive) carries `COPYING`, this directory, `deps.lock.json` (exact dependency
revisions and archive SHA-256) and `host-header-provenance.json` (the CPython header set
of the recorded macOS Blender binding build). `linux-host-header-provenance.json` records
the separate Linux Blender 5.2.1 header recipe, all header hashes and its pinned Blender
SDK revision. Linux packaging verifies the actual build inventory and every parsed
configuration macro against the target interpreter; its report identifies the selected
receipt. Other target header sets require their own build evidence. Third-party
texts here are copied unchanged from the
verified upstream archives listed in `deps.lock.json` (or, for CPython, from the Python
3.13.13 bundled with Blender 5.2.1); no third-party source file was edited.

| Component | Version | License | Compiled into | Notice files |
|---|---|---|---|---|
| CleanSeams | 0.3.0 extension / 0.2.0 engine | GPL-3.0-or-later | everything | `../COPYING`, `cleanseams/COPYING` |
| Eigen | 3.4.0 | MPL-2.0 (some files BSD-3-Clause, Apache-2.0, LGPL, MINPACK; see COPYING.README) | libcleanseams_core (headers) | `eigen/COPYING.*` |
| SuperLU (serial) | 7.0.0 | BSD-3-Clause (LBNL); its bundled CBLAS is the f2c translation of the netlib reference BLAS distributed inside SuperLU under the same file | libcleanseams_core (static) | `superlu/License.txt` |
| nlohmann/json | 3.12.0 | MIT | libcleanseams_core, cleanseams_cli, cleanseams_native | `nlohmann_json/LICENSE.MIT` |
| pybind11 | 3.0.1 | BSD-3-Clause | cleanseams_native | `pybind11/LICENSE` |
| NumPy (adapted introsort/heapsort, `src/quality/quantiles.cpp`) | 2.4.4 | BSD-3-Clause | libcleanseams_core | `numpy-2.4.4.txt` (provenance: `../third_party/quality_sort_provenance.json`) |
| CPython headers (inline functions) | 3.13.13 | PSF-2.0 | cleanseams_native | `cpython/LICENSE.txt` |

Blender supplies its own Python interpreter and NumPy installation.

Corresponding source and the build recipe: [SOURCES.md](SOURCES.md).
