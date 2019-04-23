#Configuration file for TraCeR#

[tool_locations]
#paths to tools used by TraCeR for alignment, quantitation, etc
bowtie2_path = /scg/apps/software/anaconda/3_5.0.1_20180125/envs/tracer_teichlab_2018-03-08/bin/bowtie2
bowtie2-build_path = /scg/apps/software/anaconda/3_5.0.1_20180125/envs/tracer_teichlab_2018-03-08/bin/bowtie2-build
igblastn_path = /scg/apps/software/anaconda/3_5.0.1_20180125/envs/tracer_teichlab_2018-03-08/bin/igblastn
makeblastdb_path = /scg/apps/software/anaconda/3_5.0.1_20180125/envs/tracer_teichlab_2018-03-08/bin/makeblastdb
kallisto_path = /scg/apps/software/anaconda/3_5.0.1_20180125/envs/tracer_teichlab_2018-03-08/bin/kallisto
salmon_path = /scg/apps/software/anaconda/3_5.0.1_20180125/envs/tracer_teichlab_2018-03-08/bin/salmon
trinity_path = /scg/apps/software/anaconda/3_5.0.1_20180125/envs/tracer_teichlab_2018-03-08/bin/Trinity
dot_path = /scg/apps/software/anaconda/3_5.0.1_20180125/envs/tracer_teichlab_2018-03-08/bin/dot
neato_path = /scg/apps/software/anaconda/3_5.0.1_20180125/envs/tracer_teichlab_2018-03-08/bin/neato


[trinity_options]
#line below specifies maximum memory for Trinity Jellyfish component. Set it appropriately for your environment.
max_jellyfish_memory = 1G

#uncomment the line below if you've got a configuration file for Trinity to use a computing grid
#trinity_grid_conf = /path/to/trinity/grid.conf

#uncomment the line below to explicitly specify Trinity version. Options are '1' or '2'
#trinity_version = 2

#### <---- beginning of trinity specialized options
# additional Trinity options in case you're dealing with very short reads (say 25 base reads)
# and want to achieve high sensitivity with questionable specificity, in other words
# trying to extract whatever you can from the data you have:

###  note, default Trinity kmer length is 25
#trinity_kmer_length = 17

### below stops trinity at the initial inchworm (greedy kmer extension) step.
#inchworm_only = True

#### end of trinity specialized options ---->


[IgBlast_options]
igblast_seqtype = TCR

[base_transcriptomes]
# reference transcriptomes for kallisto/salmon.  Just point to the raw transcriptome fasta files.
Mmus = /home/buckley7/Buckley/1.References/transcriptomes/gencode.vM16.transcripts.fa
Hsap = /path/to/kallisto/transcriptome_for_Hsap

[salmon_base_indices]
# salmon indices created from [base_transcriptomes] above; needed only when option --small_index is used
Mmus = /path/to/salmon/index_for_Mmus
Hsap = /path/to/salmon/index_for_Hsap

[kallisto_base_indices]
# kallisto indices created from [base_transcriptomes] above; needed only when option --small_index is used
Mmus = /home/buckley7/Buckley/1.References/kallisto/gencode.vM16.03152018
Hsap = /path/to/kallisto/index_for_Hsap

[salmon_options]
# line below specifies type of sequencing library for Salmon; if not specified, automatic detection (--libType A) is used
#libType = A

# line below specifies minimum acceptable length for valid match in salmon's quasi mapping; if not specified, default value of 31 is used
#kmerLen = 31

[tracer_location]
#Path to where TraCeR was originally downloaded
tracer_path = /home/buckley7/Buckley/5.BenNSCProject/tracer