## Phylogenetic Tree Construction

![1200px-PhylogeneticTree,_Woese_1990](https://github.com/cmatth25/722_Group/assets/101157734/aac7c31f-dfc6-40c8-8360-9c3ec5728088)
<sub><sup>Maulucioni, CC BY-SA 3.0 <https://creativecommons.org/licenses/by-sa/3.0>, via Wikimedia Commons adapted from Woese CR, Kandler O, Wheelis ML (1990). *Proc Natl Acad Sci USA*. **87 (12)**: 4576–79. doi:10.1073/pnas.87.12.4576.</sup></sub>


### Sequences for alignment (or not)

#### please log into info2020


Remember using 16S in the metagenomics tutorial for microbial identification? It's a great ortholog for phylogenetic tree making for the same reason. Carl Woese famously used 16S to differentiate between Bacteria and Archea using 16S and used the small rRNA subunit more generally to propose the above Tree of Life.

Let's retrieve some 16S sequences from experimental isolates and some related type strains to get a better idea what we're working with and where they fall in the Actinobacterial tree. Oftentimes, when we're dealing with species and especially at the strain level, 16S just isn't enough (and I can confirm it simply isn't here either), so let's get 23S, part of the large ribosomal subunit, and concatenate the sequences for a multigene alignment.  Here's a link to some rRNA sequences.

``` 
mkdir trees ; cd trees
```
```
ln -s /2/scratch/CraigM/Group/test/*.fna . 
```
Lets check to make sure we've only got 1 sequence (barrnap may isolates multiple rRNA sequences) in each file. While we're at it, let's check they're about the right length, 16S is ~1.5 kb and 23S is ~3 kb with some variability.
```
for file in *S.fna; do echo "${file}" | tr '\n' '\t'; grep '^>' ${file} | wc -l | tr '\n' '\t';  grep -v '^>' ${file} | wc -c; done
```
Alright, let's concatenate... using paste
```
for file in `ls *_16S.fna | sed 's/_16S.fna@//g' `; do paste -d'\0' ${file}_16S.fna  ${file}_23S.fna > ${file}_combinedS.fna; done
```
```
mkdir 16S_23S/ ; mv *dS.fna 16S_23S/; cd 16S_23S/
```

Might be worth checking that you've got your sequences together and they're all about ~4.5 kb using the same script as before.

This could be done for many such selected marker genes, or some alignments take the entire CDS translated into AA sequences. 

An alternative aggregate method, building a concensus tree, would take many orthologous genes, construct trees from each gene, and the final tree would be the result of the plurality "vote" from each gene tree, but we would want many more than 2. Many tree building packages will pull out coding sequences from annotated genomes and allow for concatenated or concensus multigene trees. 


Back to the point, let's rename the fasta headers to make our final labels easier. This makes things easier for the sake of labelling but is not necessarily advisable since you may be losing some information, but since the paste function to add 16S and 23S together messed it up anyways, you might at least want to re add accurate metadata.

```
head -n1 SID10270_combinedS.fna
```
```
for file in *dS.fna; do sed -i "1s/.*/>${file%_combinedS.fna}/" $file; done
```
```
head -n1 SID10270_combinedS.fna
```
If you want a better idea of how I extracted these sequences from unannotated genome assemblies, check out the Appendix.

Let's pull them together into a multi-fasta file simply using cat, align them with mafft (various aligners such as muscle or clustal could work here but take a bit longer) and take a look.  

```
cat *dS.fna > 16S_23S_multi.fna
```
We can take a closer look and make sure we've got everything expected before align.
```
grep '^>' 16S_23S_multi.fna
```

```
/usr/local-centos6/clustal/clustalo_1.2.4 -i 16S_23S_multi.fna --outfmt=fa -o 16S_23S_al.afa
```

You'll want to note that whether you have phylip or fasta alignments, if your sequences didn't already have line breaks, they will now. 

Taking a closer look at the alignment, you'll see things get a bit messy at either end and where the genes have been concatenated together.

```
less 16S_23S_al.afa
```

There are plenty of good ways to handle gaps but since we concatenated 2 genes, and things can get a little messy at either end, I'd like to remove the heavily gapped areas them with trimAl. GBlocks is a more interactive alternative, if that's for you.

```
/usr/local-centos6/trimAl/trimAlv1.4/source/trimal -in 16S_23S_al.afa -out 16S_23S_trimal.afa -fasta -gt 0.25
```

#### SNP identification

Alignments are great and there are various packages to identify SNPs such as SNP-sites to pull SNPs out of the alignment, or snippy to pull out SNPs compared to a reference genome. We're going to move on from alignments, I hope you learned something about them (you can plug that alignment in a little later if there's time) but what if we want to compare the genomes of organisms where alignments aren't reliable? Microbial genomes are rife HGT, making alignments difficult or computationally expensive. How do we identify orthologs or SNPs?

Hopefully everyone remembers k-mers from the first metagenomics tutorial. Where kraken used exact k-mer matches (k-mer length 35) to classify sequences to a taxonomic rank based on a database, here kSNP uses split k-mers.

kSNP finds all kmers of k length for each sequence and identifies kmers that differ at exactly the mid-point. k=9 example below.

>seq1  
AAAA**T**AAAA  
>seq2  
AAAA**G**AAAA  

The appropriate k-mer length isn't necessarily straightforward, but kSNP provides a tool, Kchooser, to search for the optimal length for your sequences and to assess the FCK, fraction of core k-mers. Core k-mers are those that appear in all of your sequences and reflects the diversity of samples, and the fraction of core k-mers is the core k-mers/all k-mers. Low FCK scores means your samples are very diverse. The higher the FCK, the more reliable your kSNP tree will be.

Note, kSNP4.1pgk needs adding to path to use Kchooser because it is a stupid program.

It is our last class and Ben mentioned showing the class adding a vairable to the path so...

Before we mess around with your path in your .bash_profile, might as well take a look at what it's supposed to look like.

```
cat ~/.bash_profile
```
now we can add kSNP to your path
```
echo 'export PATH=/usr/local/kSNP/kSNP4.1pkg:$PATH' >> ~/.bash_profile
```
">>" to append is VERY important here. When in doubt, use your text editor of choice so you don't accidently mess **everything** up.

source bash profile to make sure it's added to your path now, but moving forward when you login you should be good to go.
```
source ~/.bash_profile
```
might as well take a quick look again
```
cat ~/.bash_profile
```

Before we can run Kchooser, we need an input file for kSNP (and Kchooser). kSNP comes with a build in untility for this(MakeKSNP4infile), but your unaligned files need to be the only thing in the directory. Instead you can run this in the directory with your files and whatever else you want:

```
for file in `ls *dS.fna `; do echo -e "$PWD/${file}\t${file::-14}" >> kSNP_input.txt;  done
```
now Kchooser can tell us what k-mer value to use and whether or not it's a good idea to run kSNP at all. It probably isn't because we have an alignment already and these sequences are short and well conserved but we'll continue on to familiarize you with some genome scale tools. 
```
/usr/local/kSNP/kSNP4.1pkg/Kchooser4 -in kSNP_input.txt
```

I seem to remember being told this is a bad idea before, but I think that was before I ran any code at all. Alas, we carry on. 

kSNP will build a tree for us, which is great to get an idea of what the tree looks like, but it doesn't provide much in the way of tree building options or any visualization.
```
mkdir kSNP_out
```
```
/usr/local/kSNP/kSNP4.1pkg/kSNP4 -k 11 -in kSNP_input.txt -min_frac 0.75 -core -vcf -NJ -outdir kSNP_out
```
-min_frac will provide outputs where only SNPs that appear in atleast 75% of the samples will be included and -core will provide outputs where only SNPs present in all samples are included.

-NJ will output a NJ tree and the distance matrix used to make it based on p-distance

```
cd kSNP_out/; less NJ.dist.matrix
```

kSNP lacks much of the functionality of packages designed specifically for tree building, like bootstrapping or changing substitution models. Luckily it does provide an output we can plug into other software packages. -vcf outputs a variant call format file, a file commonly used in SNP analysis that provides info about the position, alleles, etc, which you are welcome to inspect. I'm just going to stick with the fasta alignment for simplicity's sake.

You can see since I included -core and and -min_frac arguments, it's provided quite a number of outputs.

```  
less core_SNPs_matrix.fasta
```
```
cat COUNT_coreSNPs
```
```
less SNPs_in_majority0.75_matrix.fasta
```

There's a link to a kSNP tree built from these samples' whole genomes to try and get a better idea about some of the very closely related species in the Appendix.

For well characterized species, SNP databases are available, and SNP-sites can pull SNPs from an alignment, such as the one we prepared ourselves earlier.

#### Maximum liklihood tree with bootstrapping

We're going to keep it relatively simple and build a couple of trees from the alignment and the SNPs.

Let's start with SNPs

```
mkdir trees
```
First, let's test some of our tree building assumptions. I went through the MFP mode of iqtree that picks the best model for you based on log-liklihood scores and the Bayesian information criterion, which I won't get into. We can discuss the model it chose after.

While good, tree building assumptions are TRSH: treelikeness(bases descend from the same bifurcating tree), reversibility (mutations are just as likely to occur in reverse), stationarity (base frequencies remain constant) and homogeneity (subsitution rates remain constant). Here we use iqtree2's symtest to test for homogeneity and stationarity.
sym tests both, marginal symmetry tests stationarity and internal symmetry tests homogeneity.

```
 /usr/local-centos6/iqtree/version2.2/iqtree2 -s SNPs_in_majority0.75_matrix.fasta -m GTR+ASC+R2 -b 20 --symtest-only 
```
```
cat SNPs_in_majority0.75_matrix.fasta.symtest.csv
```
Now let's make a tree!

```
/usr/local-centos6/iqtree/version2.2/iqtree2 -s SNPs_in_majority0.75_matrix.fasta -m GTR+ASC+R2 -b 20 -redo
```
You'll notice the model includes +ASC, that is to account for ascertainment bias, because we've picked out all the variable sites and let behind invariable sites, we've biased the alignment and also increased the distances we would expect to find, since there will only be sites with variation. We've also included 20 bootstraps to support the tree, a number you would ideally have at least 100, if not thousands.

```
less SNPs_in_majority0.75_matrix.fasta.iqtree
```
Or for the just the tree, you can find it at the bottom in Newick format, which looks something like ((raccoon, bear),((sea_lion,seal),((monkey,cat), weasel)),dog) and was named after a Restaurant in New Hampshire of all things.

```
tail SNPs_in_majority0.75_matrix.fasta.iqtree
```

#### Visualize it

Great, now there are great software packages in the command line and R that can build publication ready trees such as FigTree or ggtree, but let's play around with our tree in the interactive Tree of Life (iTOL) GUI. It's simple, you can upload tree files (newick, plain text, phyloXML, .jplace from RaxML, .qza from QIIME), but even easier for our sake, copy and paste a newick tree.

https://itol.embl.de/

Now let's try annotating.


iTOL annotation file templates: 
https://itol.embl.de/help/example_data.zip 

Consensus tree in newick format:

(Actinacidaphila_yeochonensis:0.1666831938,(((((Bifidobacterium_bifidum:0.1428754323,Chloroflexus_aggregans:0.2932249821)98:0.5056190704,Mycobacterium_tuberculosis:0.2562105481)99:0.3816801535,(Micromonospora_antibiotica:0.1178612370,(Micromonospora_echinospora:0.0451309998,(SID11719:0.0048610328,SID8012:0.0070801205)100:0.0362111446)100:0.0818217338)100:0.1298468109)55:0.0339660549,(Pseudonocardia_abyssalis:0.1983869364,(Pseudonocardia_autotrophica:0.1378758308,Pseudonocardia_sp_EC080610-09:0.1332401010)100:0.1464868849)100:0.1041144798)98:0.0911355511,Kitasatospora_aureofaciens:0.1465795573)55:0.0214911296,(((((SID10270:0.0229261328,(SID7918:0.0000026526,SID9923:0.0000026496)100:0.0363756763)45:0.0056938377,Streptomyces_coelicolor:0.0273292280)100:0.0371432544,(SID11753:0.0896346430,(SID7959:0.0023950156,SID9905:0.0028849378)100:0.0914221756)63:0.0175993353)99:0.0358787923,((((SID11370:0.0017832824,SID7966:0.0018423293)100:0.0581998719,(SID7932:0.0035768388,SID8462:0.0017824813)100:0.0189778288)100:0.0700341367,SID7965:0.1749979952)50:0.0226990834,(Streptomyces_griseus:0.1255275073,Streptomyces_kanamyceticus:0.1362039465)94:0.0341240577)57:0.0167198815)93:0.0244699743,Streptomyces_albus:0.2018588730)90:0.0282543333)

## Appendix

### tree from alignment
I wanted to play around with an alignment, in part because if we didn't use just a gene or 2 we would have been waiting all afternoon. Heres a tree from the alignment with partitions for each gene, worth noting if you're planning on making multigene alignment trees where different genes are best modelled differently. 16S and 23S we know have some peculiarities, so it should now these alignments break stationarity and homogeneity assumptions even with partitioning.

If you have bad data that you feel you need to include, the suggestion is test the results with and without the bad data partitions to get an idea of how they're impacting the tree, ideally, they're not.

```
cd ..
```
Just so we didn't waste all that time aligning, we can look at partitioning the data alignment we made earlier. 
```
nano partition.txt
```
in partition.txt add the following, which defines the two genes. You can use partition to define base positions using a "\3" following the range of you CDS in frame (each position would look like. But I'm not sure about the reading frame here so don't do that. The program is smarter than me any ways and figures out what works best.
part1 = 1-100\3, 2-101\3
part2 = 3-102\3.
```
DNA, part1 = 1-1528
DNA, part2 = 1529-4654
```
We're going to use the UFBoot2 (-B) to approximate 1000 bootstraps, rather than actually performing them, its very fast, but again, an appproximation and not a true bootstrap. I didn't have time to look too deeply into it.
```
/usr/local-centos6/iqtree/version2.2/iqtree2 -s 16S_23S_trimal.afa -m GTR+F+I+I+R5 -B 1000 
```
```
less 16S_23S_trimal.afa.iqtree
```


#### rRNA retrieval by barrnap

This tool is available on compute canada and runs very quickly, extracting prokaryotic or eukaryotic rRNA sequences from unannotated genomes, outputting .gff or .fasta files. Here's a bash script that would work in compute canada. If your genome is annotated already, you can use faidx from samtools to extract your GOI from known ranges. Which yuou should be able to do in info.

```
#!/bin/bash
module load StdEnv/2020
module load barrnap/0.9

# variable 1 is complete path to directory containing sequences, including final "/"

for file in ${1}*.fna
do
        barrnap --kingdom bac --threads 15 --outseq ${file::-4}_rrna.fna < ${file} > ${file::-4}_rrna.gff
done

for file in ${1}*_rrna.fna
do
  cat ${file::-9}_rrna.fna | grep -A 1 '^>16' | head -n 2 > ${file::-9}_16S.fna
  cat ${file::-9}_rrna.fna | grep -A 1 '^>23' | head -n 2 > ${file::-9}_23S.fna
        paste -d'\0' ${file::-9}_16S.fna  ${file::-9}_23S.fna > ${file::-9}_combinedS.fna
done
mkdir 16_23S/
mv *S.fna 16_23S/

```


##### whole genome example using the tools covered here

here's a link to a tree built with the code that follows

```
ln-s /2/scratch/CraigM/Group/test/16S_23S/kSNP_out/WholeGenomekSNPIQtreeExample.iqtree .
```
the code that follows:
```
/usr/local/kSNP/kSNP4.1pkg/kSNP4 -k 17 -in kSNP_input.txt -min_frac 0.75 -core -ML -NJ -outdir kSNP_out
```
```
/usr/local-centos6/iqtree/version2.2/iqtree2 -s ../../../Genomes/Assembled/kSNP_out/SNPs_in_majority0.75_matrix.fasta -m MFP -b 100 
```
This only took about ~18 min


#### phylogenetic analaysis in R using APE

ape is a great R package for building trees from shorter sequences and alignments of smaller sample sizes and has great functionality (still much larger alignments and greater sequence numbers than this example, but will be much slower than many of the command line tools above).

Here's a little preliminary study that I initially used to look into the 16S sequences using ape to calculate a p-distance matrix, pull out the max p-value to get an idea of what tree-building method would be appropriate, and go from there. 

```
library(ape)
library(dplyr)
library(tidyverse)

setwd("")

#load alignment
seqsCombined <- read.dna("16S_23S_al.afa", format = "fasta")

#trim gaps (all)
seqsCombined < del.colgapsonly(seqsCombined)

#find p-distances
distmxCombined <- dist.dna(seqsCombined, model = "raw") 
#as matric and view
distmxtpCombined <- dist.dna(seqsCombined, model = "raw", as.matrix = TRUE) 
distmxtpCombined

#max p-distance?
maxpCombined <- max(as.numeric(unlist(distmxtpCombined)))
maxpCombined

#TN93 tree

distmxTNCombined <- dist.dna(seqsCombined, model = "TN93") 

#neighbor joining
njCombined <- nj(distmxTNCombined)

#root tree with outgroup
njCombinedRooted <- root(njCombined, outgroup = "Chloroflexus_aggregans_combinedS", resolve.root = TRUE)
#check it's rooted
is.rooted(njCombinedRooted)

#a tibble tree, another tree format used in ape and R packages
njtibble <- as_tibble(njCombinedRooted)
njtibble

#write Newick to file
write.tree(njCombinedRooted,file="NewickExample.txt")

#couple plot formats
plot(njCombinedRooted, "f", use.edge.length = FALSE, xpd = TRUE, no.margin = FALSE)

plot(njCombinedRooted,"p", use.edge.length = TRUE, xpd = TRUE, no.margin = FALSE)

```
