#pragma once
// Process-wide settings made at the start of every program (before any HDF5 call).

#include <cstdlib>

namespace ena {

// HDF5 (1.10 and later) locks every file it opens. On NFS a lock makes the client
// discard its cached copy of the file, so each open re-reads the whole file over the
// network (measured on an NFS workstation: 4.3 s instead of 0.3 s for a 620 MB
// snapshot). The programs never write a file another process reads at the same time
// (inputs are only read; outputs are written by one rank), so locking is switched
// off. An HDF5_USE_FILE_LOCKING value set by the user is kept.
inline void prepare_runtime() { setenv("HDF5_USE_FILE_LOCKING", "FALSE", 0); }

}  // namespace ena
