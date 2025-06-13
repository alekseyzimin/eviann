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

dist_bin_SCRIPTS = src/ez_maker.sh \
	src/my_maker_functional_fasta \
	src/my_maker_functional_gff \
	src/combine_gene_protein_gff.pl \
	src/eviann.sh \
	src/eviprot.sh \
	src/fix_suspect_introns.pl \
	src/fix_splice_junctions.pl  \
	src/add_cds_to_gff.pl \
	src/check_cds.pl \
	src/filter_by_local_abundance.pl \
	src/filter_by_class_code.pl \
	src/filter_unused_proteins.pl \
	src/extract_utr_transcripts.pl \
	src/detect_readthroughs.pl \
	src/detect_readthrough_exons.pl \
	src/remove_readthrough_exons.pl \
	src/reassign_transcripts.pl \
	src/assign_class_code.pl \
	src/tblToGFF.pl \
	src/convert_trna.sh \
	src/add_features.pl \
	src/convert_miniprot_gff.pl \
	src/compute_junction_scores.pl \
	src/compute_junction_scores_bed.pl \
	src/compute_junction_counts.pl \
	src/score_transcripts_with_hmms.pl \
	src/fix_palign_splice_sites.pl \
	src/filter_junctions.pl \
	exonerate \
	makeblastdb blastp tblastn \
	gffread gffcompare trmap \
	stringtie \
	TransDecoder-v5.7.1a.tar.gz \
	samtools \
	miniprot.tgz
