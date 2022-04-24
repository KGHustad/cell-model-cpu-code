# Rules to update the auto-generated command-line argument parser by running gengetopt
# Please use CMake to build the source code
cmdline_files_basename = dump_trace_cmdline

all: $(cmdline_files_basename).c

$(cmdline_files_basename).c $(cmdline_files_basename).h: $(cmdline_files_basename).ggo
	gengetopt -l -C -F $(cmdline_files_basename) < $<
