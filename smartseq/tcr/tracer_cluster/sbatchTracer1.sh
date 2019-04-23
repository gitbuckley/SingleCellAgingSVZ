#!/bin/bash
#SBATCH --account=abrunet1
#SBATCH --job-name="tracerConda"
#SBATCH --nodes 1 # Number of tasks or threads
#SBATCH -t 1-00:01 # Runtime in D-HH:MM
#SBATCH --mem-per-cpu=16G # see also --mem 
#SBATCH -o out.%j # File to which STDOUT will be written
#SBATCH -e err.%j # File to which STDERR will be written
#SBATCH --mail-user=buckley7@stanford.edu
#SBATCH --mail-type=END,FAIL # What type of mail to send

# Purpose run "tracer assemble" on fastq files associated with a single cell. 

date

module load anaconda
source activate tracer_teichlab_2018-03-08

cd /srv/gsfs0/projects/brunet/Buckley/5.BenNSCProject/2.Tcell/

tracer assemble -c 5.Tracer/CONFIG \
				0.Raw/${F}_R1_001.fastq \
				0.Raw/${F}_R2_001.fastq \
				$F \
				5.Tracer/Output_Conda_Mar15 

date
echo "Finished"

## Environment packages:
# > conda list
# biopython                 1.70             py36h637b7d7_0
# blast                     2.6.0               boost1.64_2    bioconda
# boost                     1.64.0                   py36_4    conda-forge
# boost-cpp                 1.64.0                        1    conda-forge
# bowtie                    1.2.2            py36pl5.22.0_0    bioconda
# bowtie2                   2.3.4.1          py36pl5.22.0_0    bioconda
# bzip2                     1.0.6                h9a117a8_4
# ca-certificates           2017.08.26           h1d4fec5_0
# cairo                     1.14.12              h77bcde2_0
# certifi                   2018.1.18                py36_0
# collectl                  4.0.4                pl5.22.0_3    bioconda
# curl                      7.58.0               h84994c4_0
# cycler                    0.10.0           py36h93f1223_0
# cython                    0.27.3           py36h1860423_0
# dbus                      1.12.2               hc3f9b76_1
# decorator                 4.2.1                    py36_0
# expat                     2.2.5                he0dffb1_0
# fastool                   0.1.4                         2    bioconda
# fontconfig                2.12.6               h49f89f6_0
# freetype                  2.8                  hab7d2ae_1
# future                    0.16.0                   py36_1
# glib                      2.53.6               h5d9569c_2
# graphite2                 1.3.10               hf63cedd_1
# graphviz                  2.40.1               h25d223c_0
# gst-plugins-base          1.12.4               h33fb286_0
# gstreamer                 1.12.4               hb53b477_0
# harfbuzz                  1.7.4                hc5b324e_0
# hdf5                      1.8.17                        2
# icu                       58.2                 h9c2bf20_1
# igblast                   1.7.0                pl5.22.0_0    bioconda
# intel-openmp              2018.0.0             hc7b2577_8
# ipython                   6.2.1            py36h88c514a_1
# ipython_genutils          0.2.0            py36hb52b0d5_0
# jedi                      0.11.1                   py36_0
# jellyfish                 2.2.6                         0    bioconda
# jemalloc                  4.5.0                         0    bioconda
# jpeg                      9b                   h024ee3a_2
# kallisto                  0.44.0             hdf51.8.17_1    bioconda
# kiwisolver                1.0.1            py36h764f252_0
# libcurl                   7.58.0               h1ad7b7a_0
# libedit                   3.1                  heed3624_0
# libffi                    3.2.1                hd88cf55_4
# libgcc                    7.2.0                h69d50b8_2
# libgcc-ng                 7.2.0                hdf63c60_3
# libgfortran-ng            7.2.0                hdf63c60_3
# libpng                    1.6.34               hb9fc6fc_0
# libssh2                   1.8.0                h9cfc8f7_4
# libstdcxx-ng              7.2.0                hdf63c60_3
# libtiff                   4.0.9                h28f6b97_0
# libtool                   2.4.6                h544aabb_3
# libxcb                    1.12                 hcd93eb1_4
# libxml2                   2.9.7                h26e45fe_0
# matplotlib                2.2.0            py36hbc4b006_0
# mkl                       2018.0.1             h19d6760_4
# mock                      2.0.0            py36h3c5bf6c_0
# ncurses                   6.0                  h9df7e31_2
# networkx                  2.1                      py36_0
# nose                      1.3.7            py36hcdf7029_2
# numpy                     1.14.1           py36ha266831_2
# openjdk                   8.0.121                       1
# openssl                   1.0.2n               hb7f436b_0
# pandas                    0.22.0           py36hf484d3e_0
# pango                     1.41.0               hd475d92_0
# parafly                   r2013_01_21                   1    bioconda
# parso                     0.1.1            py36h35f843b_0
# patsy                     0.5.0                    py36_0
# pbr                       3.1.1            py36hb5f6b33_0
# pcre                      8.41                 hc27e229_1
# perl                      5.22.0.1                      0    conda-forge
# perl-app-cpanminus        1.7039               pl5.22.0_3    bioconda
# perl-module-build         0.4224               pl5.22.0_0    bioconda
# pexpect                   4.4.0                    py36_0
# pickleshare               0.7.4            py36h63277f8_0
# pip                       9.0.1                    py36_5
# pixman                    0.34.0               hceecf20_3
# prettytable               0.7.2                    py36_1    conda-forge
# prompt_toolkit            1.0.15           py36h17d85b1_0
# ptyprocess                0.5.2            py36h69acd42_0
# pydotplus                 2.0.2                    py36_0
# pygments                  2.2.0            py36h0d3125c_0
# pyparsing                 2.2.0            py36hee85983_1
# pyqt                      5.6.0            py36h0386399_5
# python                    3.6.4                hc3d631a_1
# python-dateutil           2.7.0                    py36_0
# python-levenshtein        0.12.0                   py36_1    bioconda
# pytz                      2018.3                   py36_0
# qt                        5.6.2               hd25b39d_14
# readline                  7.0                  ha6073c6_4
# salmon                    0.9.1                         1    bioconda
# samtools                  1.7                           1    bioconda
# scipy                     1.0.0            py36hbf646e7_0
# seaborn                   0.8.1            py36hfad7ec4_0
# setuptools                38.5.1                   py36_0
# simplegeneric             0.8.1                    py36_2
# sip                       4.18.1           py36h51ed4ed_2
# six                       1.11.0           py36h372c433_1
# slclust                   02022010                      2    bioconda
# sqlite                    3.22.0               h1bed415_0
# statsmodels               0.8.0            py36h8533d0b_0
# tbb                       4.4_20150728                  0    bioconda
# tk                        8.6.7                hc745277_3
# tornado                   5.0                      py36_0
# tracer                    0.5                       <pip>
# traitlets                 4.3.2            py36h674d592_0
# trimmomatic               0.36                          5    bioconda
# trinity                   2.5.1                         1    bioconda
# wcwidth                   0.1.7            py36hdf4376a_0
# wheel                     0.30.0           py36hfd4bba0_1
# xz                        5.2.3                h55aa19d_2
# zlib                      1.2.11               ha838bed_2
