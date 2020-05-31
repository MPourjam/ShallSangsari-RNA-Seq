# ShallSangsari-RNA-Seq
##### Analysis of ovary Transcriptome using DER Finder Package

The process is insired from the article "Incomplete annotation of OMIM genes .. ". The core fucntion is the "NewRegions" function aiming the goal of finding the best MCC and MRR for creating the ReginMatrix aligning the most on the previously annotated exons accorsing to the guiding article.\
The process is faciliated to be implemented easily on a low-performance computer by implementing the iteration outside R environment to duck the memory overloading problem in R. The goal has achieved through a perl document called "RPerl.pl".
The function is almost ready to use except for two main analysis.
I) Firstly, the filter for exons shorter than 3bp which is almost finished. functions.R line 81 / RegionMat_MCC_MRG
II) Secondly, filtering regions laying over multiple exons to avoid further complications might arise through downstream analysis.
