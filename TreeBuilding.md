## Phylogenetic Tree Construction

![1200px-PhylogeneticTree,_Woese_1990](https://github.com/cmatth25/722_Group/assets/101157734/aac7c31f-dfc6-40c8-8360-9c3ec5728088)


#### Sequences for alignment (or not)

Remember using 16S in the metagenomics tutorial for microbial identification? It's a great ortholog for phylogenetic tree making for the same reason.  
Let's retrieve some 16S sequences from experimental isolates and some related type strains to get a better idea what we're working with and where they fall in the Actinobacterial tree. Here's a link to some rRNA sequences 

```
{}
```

Oftentimes, when we're dealing with species and especially at the strain level, 16S just isn't enough (and I can confirm it simply isn't here either), so let's get 23S, part of the large ribosomal subunit, and concatenate the sequences for a multigene alignment. 

```
{}
```
This could be done for many such selected marker genes, or the entire CDS translated into AA sequences. 

An alternative aggregate method, building a concensus tree, would take many orthologous genes, construct trees from each gene, and the final tree would be the result of the plurality "vote" from each gene tree, but we would want many more than 2. Many tree building packages will pull out coding sequences from annotated genomes and allow for concatenated or concensus multigene trees. 

Back to the point, we'll rename the fasta headers to make our final labels easier.

```
{}
```

If you want a better idea of how I extracted these sequences from unannotated genome assemblies, check out the Appendix.

Let's align them and take a look.  

```
{}
```

There are plenty of good ways to handle gaps but since we concatenated 2 genes, and things can get a little messy at either end, I'd like to remove them with trimAl. GBlocks is a more interactive alternative, if that's for you.

```
#trimal
{}
```
OR
```
#gblocks
{}
```

#### SNP identification

Alignments are great and there are various packages to pull out SNPs from alignments or compared to a reference sequence (we'll plug that alignment in a little later) but what if we want to compare the genomes of organisms where alignments aren't reliable? Microbial genomes are rife HGT, making alignments difficult or computationally expensive. How do we identify orthologs or SNPs?

We can use a split k-mer approach offered by kSNP. kSNP finds all kmers of k length for each sequence and identifies kmers that differ at exactly the mid-point. k=9 example below.

>seq1  
AAAA**T**AAAA  
>seq2  
AAAA**G**AAAA  

The appropriate k-mer length isn't necessarily straightforward, but kSNP provides a tool, Kchooser, to search for the optimal length for your sequences and to assess the FCK, fraction of core k-mers. Core k-mers are those that appear in all of your sequences and reflects the diversity of samples, and the fraction of core k-mers is the core k-mers/all k-mers. Low FCK scores means your samples are very diverse. The higher the FCK, the more reliable your kSNP tree will be.

```
{}
```

kSNP will build a tree for us, which is great to get an idea of what the tree looks like, but it doesn't provide much in the way of tree building options or any visualization.

```
{}
```

kSNP lacks much of the functionality of packages designed specifically for tree building, like bootstrapping or changing substitution models. Luckily it does provide an output we can plug into other software packages.

```  
{}
```

There's a link to a kSNP tree built from these samples' whole genomes to try and get a better idea about some of the very closely related species in the Appendix.

For well characterized species, SNP databases are available

#### Maximum liklihood tree with bootstrapping

```
{}
```


#### Visualize it

Great, now there are great software packages in the command line and R that can build publication ready trees such as FigTree or ggtree, but let's play around with our tree in the interactive Tree of Life (iTOL) GUI. It's simple, you can upload tree files (newick, plain text, phyloXML, .jplace from RaxML, .qza from QIIME), but even easier for our sake, copy and paste a newick tree.

https://itol.embl.de/

```
{}
```

## Appendix

#### rRNA retrieval by barrnap

This tool is available on compute canada and runs very quickly, extracting prokaryotic or eukaryotic rRNA sequences from unannotated genomes, outputting .gff or .fasta files.

```
{}
```

#### GOI

If you're looking for something else, there may be a specific tool similar to barrnap. Alternatively, if your genome is annotated you can use gffread (in cufflinks on info).

```
{/usr/local-centos6/cufflinks1.1.0/gffread -}
```

#### kSNP, kChooser

What's the appropriate kmer length for kSNP? We'll it depends. Luckily, they made a tool for that too, kChooser.

```
{}
```

##### whole genome example using the tools covered here

here's a link to a tree built with the code that follows

```
{}
```

the code that follows:

```
{}
```

#### phylogenetic analaysis in R using APE

ape is a great package for building trees from shorter sequences and alignments of smaller sample sizes and has great functionality (still much larger alignments and greater sequence numbers than this example, but will be much slower than many of the command line tools above).

Here's a little preliminary study of the 16S sequences using ape to calculate a p-distance matrix, pull out the max p-value to get an idea of what tree-building method would be appropriate, and go from there. 

```
{}
```
