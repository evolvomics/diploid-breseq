diploid-breseq
==============

This page describes extensions to |breseq| for calling single-base substitutions and small indels in non-haploid genomes (e.g., diploid or tetraploid genomes).

.. _invoking_multiploid_mode:   

Invoking multiploid mode
-------------------------

The simplest way to invoke multiploid mode is to supply the ``-y`` option before each reference sequence file argument. For example:

.. code-block:: bash

   $ diploid-breseq -y 2 -r diploid.gbk -y 3 -r triploid.gbk reads.fastq

There should be the same number of instances of ``-y <ploidy>`` in your command line as there are occurrences ``-r <reference_file>``.

Notes:
- You can only use CONSENSUS mode with multiploid mutation calling. Why? Since a heterozygous site in a genome has a 50% allele frequency, you can't tell the difference between a situation in which 100% of a mixed cell population is a heterozygote and 50% of the cell population is homozygous for each of the two alleles at that position, at least with the short read data used by |breseq|.
- All reference sequences in a given input file must have the same ploidy. If you would like to use a mixture of haploid and diploid reference sequences, for example, you would need to split them out so that all of the haploid ones are in one input file and all diploid ones are in a second input file.

Multiploid reference genomes
------------------------------

The above option is all you need if you started an experiment with a completely homozygous genome. Usually, however, you will have some existing allelic variation in your multiploid reference genome. The basic reference sequence formats used by |breseq|, such as FASTA, GenBank, and GFF3 have no way of providing information about heterozygous sites on their own. They provide only a consensus sequence.

There are two ways of providing this information about sites that have mixed ancestral alleles to |diploid-breseq|.

Diploid reference alleles using IUPAC degenerate base symbols
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The first option is for diploid genomes (and only diploid genomes). You can use one of the normal reference genome file types that includes `IUPAC degenerate base symbols <https://en.wikipedia.org/wiki/Nucleic_acid_notation#IUPAC_notation>`_ for two bases to indicate the heterozygous base possibilities. For example if a genome is heterozygous for A/G at a position, then you use the character R for that base. You must provide the ``-y 2`` option for this reference sequence. This option only provides limited support for mixed ancestral alleles, as you cannot encode variation between the two homologous chromosomes that involves an indel, such as a deletion of a base in one chromosome (e.g., A/.) using IUPAC codes.

Multiploid reference alleles using a VCF file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The second option is to include BOTH a normal consensus reference sequence file AND a VCF describing variant alleles as separate reference sequence option, e.g., ``-r consensus.vcf -r variants.vcf``. In this case, you don't need to use the ``-y`` option to define ploidy. |diploid-breseq| will automatically assign the correct ploidy based on your VCF file to that reference genome.

**Limitations:** Currently, |diploid-breseq| will only handle single nucleotide variants or single-base indel variants in VCF files. It does not support phased alleles in multiploid genomes or phased SNP calling.

An example of the format for a VCF file recognized by |diploid-breseq| is as follows. Note that columns are separated by **tabs** in VCF files::

   ##fileformat=VCFv4.0
   ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
   #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
   NC_001133	3	.	A	G	.	.	.	GT	2/1
   NC_001133	10000	.	A	G	.	.	.	GT	2/1
   NC_001133	20000	.	C	G	.	.	.	GT	1/2
   NC_001133	27129	.	G	C	.	.	.	GT	1/2

The last column lists the allelic states at that position in the genome. So for the first line ``2/1`` means that the first chromosome has an ``G`` and the second has an ``A``. Since |diploid-breseq| does not maintain phasing information, the order of ``1/2`` versus ``2/1`` does not matter.

A tetrapoloid file might look like this::

   ##fileformat=VCFv4.0
   ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
   #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
   NC_001133	3	.	A	G,T	.	.	.	GT	2/1/3/3
   NC_001133	10000	.	A	G	.	.	.	GT	2/1/1/1
   NC_001133	20000	.	C	G,A,T	.	.	.	GT	1/2/3/4
   NC_001133	27129	.	G	C	.	.	.	GT	1/2/1/2

In this example, the genotype is ``G/A/T/T`` for the first line.

For further information, see the `VCF format (version 4.0) definition <http://www.internationalgenome.org/wiki/Analysis/vcf4.0/>`_.

Multiploid Methods
------------------

As in normal (haploid) |breseq|, variant calling by |diploid-breseq| uses Bayesian methods. An empirical error model is fit from read alignments. Then a genotype score is calculated for each column in the resulting alignment to the reference genome. The genotype caller includes five  states in the model: the four bases 'A', 'T', 'C', 'G', plus a gap '.'. The genotype caller uses a prior probability for all mixed (heterozygous) states (e.g., A/T). The sum across all of these values can be set via the command line option ``--heterozygote-prior``. The default is 1E-6. Set this to be lower if you expect heterozygosity to be rare in a sample.

For RA entries the genotype score is -log:sub:`10`(P) - log:sub:`10`(L), where `P` is the probability that the genotype is **not** the one that is called and `L` is the length of the genome. (The latter term, involving `L`, is a crude correction for multiple testing.) Thus, a score of –5 theoretically means that there is only a 10:sup:`-5` chance of an incorrect genotype call. In reality, the error model does not capture all kinds of biases and this should be used as a relative indicator of confidence in a call rather than an absolute one. You can adjust the ``--base-quality-cutoff`` option to control this stringency.

**Limitation 1:** The current version of |diploid-breseq| does *not* attempt to call larger heterozygous mutations, such as deletions greater than a few bases in one copy of a chromosome or whole-chromosome aneuploidies. It only uses RA evidence (not JC or MC evidence) for calling heterozygous mutations.

**Limitation 2:** In genomes with long homopolymer repeats there are likely to be many false-positive indels predicted in these sequences due to higher error rates in this context that are not accounted for by the current error model. It is recommended to set ``--consensus-reject-indel-homopolymer-length 5`` and ``--consensus-reject-surrounding-homopolymer-length 8`` (or even lower values) to remove these from the output.

**Limitation 3:** Several of the normal |breseq| filters for consensus mutations do not apply in multiploid mode: ``--consensus-frequency-cutoff``, ``--consensus-minimum-variant-coverage``, ``--consensus-minimum-total-coverage``, ``--consensus-minimum-variant-coverage-each-strand``, ``--consensus-minimum-total-coverage-each-strand``.

Multiploid Output
------------------

The HTML and GenomeDiff output from |diploid-breseq| has a few extensions to communicate information about heterozygous sites in an evolved genome

In the HTML output, there will be multiple lines for point mutations describing how each allele affects any protein sequence that they overlap (amino acid change and whether it is synonymous or nonsynonymous, for example).

In the GenomeDiff output, these changes for each allele are separated by slashes: ``/``. For example, a RA and SNP lines may appear like this for a diploid genome::

   SNP	30	101	NC_001133	27126	C/C
   RA	101	.	NC_001133	27126	0	T/T	C/C	consensus_score=3.2	frequency=1	new_cov=17/18;17/18	ref_cov=0/1;0/1	total_cov=17/19

An annotated version of these lines contains additional information::

   SNP	30	101	NC_001133	27126	C/C	aa_new_seq=E/E	aa_position=281	aa_ref_seq=E/E	codon_new_seq=GAG/GAG	codon_number=281	codon_position=3	codon_ref_seq=GAA/GAA	gene_name=FLO9	gene_position=843	gene_product=flocculin FLO9	gene_strand=<	genes_overlapping=FLO9	locus_tag=YAL063C	locus_tags_overlapping=YAL063C	mutation_category=snp_synonymous/snp_synonymous	position_end=27126	position_start=27126	snp_type=synonymous/synonymous	transl_table=1
   RA	101	.	NC_001133	27126	0	T/T	C/C	aa_new_seq=E/E	aa_position=281	aa_ref_seq=E/E	codon_new_seq=GAG/GAG	codon_number=281	codon_position=3	codon_ref_seq=GAA/GAA	consensus_score=3.2	frequency=1	gene_name=FLO9	gene_position=843	gene_product=flocculin FLO9	gene_strand=<	locus_tag=YAL063C	new_cov=17/18;17/18	new_seq=C/C	ref_cov=0/1;0/1	ref_seq=T	snp_type=synonymous/synonymous	total_cov=17/19	transl_table=1

**Note:** The ``new_cov`` and ``ref_cov`` entries vary from the standard use of the ``/`` divider per allele. Because in normal |breseq| the ``/`` is used to separate the top versus bottom strand of read coverage, a semicolon ``;`` is used to separate genotypes for these entries.

Multiploid Tutorial
-------------------

For testing |diploid-breseq|, we used data from a yeast evolution experiment that began with haploid, diploid, and tetraploid ancestors:

Selmecki, A. M., Maruvka, Y. E., Richmond, P. a, Guillet, M., Shoresh, N., Sorenson, A. L., De, S., Kishony, R., Michor, F., Dowell, R., Pellman, D. (2015) Polyploidy can drive rapid adaptation in yeast. *Nature* **519**: 349–352. `Pubmed <https://www.ncbi.nlm.nih.gov/pubmed/25731168>`_



