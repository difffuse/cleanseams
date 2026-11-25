# File I/O and command-line contracts

The library reads and writes geometry and case data. Packing is the host application’s
responsibility; the standalone Blender integration uses Blender’s own packer.

## OBJ (`include/cleanseams/io.hpp`, `src/output/obj.cpp`)

`parse_obj`/`read_obj` follow the source reader: one-based indices, negative indices
relative to the elements defined so far, `#` comments, `v` (first three values), `vt`
(first two), `f v[/vt[/vn]]`, every other record ignored; per-corner UVs exist only when
every face carries texture indices. Malformed indices or numbers are `ContractError`s.
As the source reads the file in text mode, a line ends at `\n`, `\r\n` or a lone `\r`, and records split on Python
`str.split()` whitespace (ASCII `\t\n\v\f\r`, `\x1c`–`\x1f`, space, and the UTF-8
encoded Unicode separators). Indices follow Python `int()` and coordinates `float()`:
one sign only (`+-1` and `--3` are refused), PEP 515 underscores between digits
(`1_0`, `1_0.2_5`), and an index beyond int64 is reported as outside the mesh, as
Python does. Non-ASCII Unicode digits, which Python also accepts, are refused.
`format_obj`/`write_obj` write `%.17g` coordinates (so doubles round-trip exactly),
`f v/t` with running UV indices and one `l a b` line per seam edge. Polygon order and
corner order are preserved.

## Case files (`src/output/case_io.cpp`)

Every JSON text from outside (case files and CLI inputs) goes through
`parse_untrusted_json` (nesting at most 128 levels); descriptor paths with control
characters are refused and unresolvable paths are `io_error` (see [SAFETY.md](SAFETY.md)).

`read_case`/`case_from_json` read `cleanseams.case` version 1 documents with inline
arrays or binary-array descriptors (`file`, `dtype`, `shape`, `byte_offset`,
`byte_length`, `sha256`). Supported dtypes are float64-le, int64-le, uint64-le and uint8.
Before any bytes are used the reader checks the field's dtype, exact
`product(shape) * itemsize == byte_length` (overflow-checked), offset + length within the
file, a relative path without `.`/`..`, backslashes or drive letters, no symlink escape
from the case directory, the SHA-256 of the bytes, finiteness of float data and the
memory budget. Indices, offsets, pins and endpoint pairs are range-checked; unknown
top-level fields (fixture metadata, options, packing specs) are kept. `write_case` writes
arrays above a size threshold as descriptors with their hashes.

## External process utility (`include/cleanseams/process.hpp`)

* `run_process`: direct argument-vector launch (`posix_spawn` on POSIX, `CreateProcessW`
  on Windows), with no shell. Windows arguments and executable paths are UTF-8; quoting
  preserves spaces, quotes, backslashes and empty arguments. Embedded NUL arguments are
  rejected. Stdin, stdout and stderr progress independently; each output has its own
  byte cap and truncation flag, including a zero cap. `ProcessOptions.cancel` is a shared
  native `CancellationToken`; cancellation before launch does not start an executable,
  and in-flight cancellation returns `cancelled=true`, separately from `timed_out`.
  A zero timeout does not launch; negative timeouts are invalid. The caller owns the
  ordinary process tree: an isolated POSIX process group or a Windows kill-on-close job
  is terminated on timeout, cancellation, exceptional cleanup and direct-child exit.
  The direct child is reaped and pipe workers joined before returning. A descendant
  must not keep a completed command's pipes open indefinitely. POSIX descendants that
  deliberately detach into another process group are outside that ownership guarantee.
  Windows uses an inherited-handle allowlist and does not permit job breakaway.
  The Windows implementation uses Microsoft's [process handle/job attributes](https://learn.microsoft.com/en-us/windows/win32/api/processthreadsapi/nf-processthreadsapi-updateprocthreadattribute)
  and [job lifetime rules](https://learn.microsoft.com/en-us/windows/win32/procthread/job-objects).
  No process-wide SIGPIPE handler is replaced by the POSIX runner. `R_process` exercises
  argument/Unicode transfer, bounded duplex streams, closed pipes, early/nonzero exit,
  deadlines, cancellation, descendant cleanup and concurrent calls.
## Command line (`cli/main.cpp`)

`--version --json`, `inspect`, `run`, `validate`, `report`, `aggregate`, `convert`
and `export`, with exit codes (0 done,
2 invalid input, 3 missing/incompatible backend, 4 no accepted result or failed
validation, 5 external process failure, 6 cancelled, 7 timeout, 8 internal). `run` uses
the native pipeline. The behavior is `--behavior legacy_010|native_baseline|enhanced_v1`
when given, else the options' `policy`, whose default is `native_baseline`; `enhanced_v1`
runs the enhanced proposal pass (an accepted result that keeps a reported strict-only defect
exits 0). The options are the `--options` file or else the case file's `options`
object, applied whole through the options contract. `run` writes a versioned `RunResult` with independent status flags, the
finite report and a human summary; `--trace` stores the progress events and `--counters`
a measurement record (pipeline wall time, stage timings, repeated-work counters; used by
benchmark tooling, never part of the result). `validate`
recomputes topology, map validity and pins from the case and the saved result instead
of trusting `accepted`. Outputs may not overwrite inputs or each other.

`report --input case --result
result [--reference other-result] [--meta meta.json] [--directions K]` writes the
source's `seam_report` of the result (mirror plane detected natively, exposure over K view
directions, 64 by default as in the source's benchmark) and its `verdict` over the run
record the source's benchmark built (acceptance, report, seam report, quad share);
`--reference` adds the agreement block, the reference's own report and verdict, and
`compare`; `--meta` carries the verdict's asset facts (`expected_islands`, `symmetry`,
`tier`). `aggregate --comparisons list.json` applies the claim rule to `compare`
records.
