AM_CXXFLAGS = -Wall -Werror -g -O2 -I$(top_srcdir)/src

# Yaggo variable
noinst_HEADERS = $(YAGGO_SOURCES)
BUILT_SOURCES = $(YAGGO_SOURCES)
DISTCLEANFILES = $(BUILT_SOURCES)

# Yaggo automatic rules with silencing
V_YAGGO = $(V_YAGGO_$(V))
V_YAGGO_ = $(V_YAGGO_$(AM_DEFAULT_VERBOSITY))
V_YAGGO_0 = @echo "  YAGGO " $@;
.yaggo.hpp:
	$(V_YAGGO)$(YAGGO) -o $@ $<
YAGGO_SOURCES = # Append all file to be built by yaggo

dist_bin_SCRIPTS = src/eugene.sh src/my_maker_functional_fasta src/my_maker_functional_gff src/combine_gene_protein_gff.pl src/eviann.sh src/protein2genome.sh src/fix_suspect_introns.pl src/fix_splice_junctions.pl  src/add_cds_to_gff.pl
