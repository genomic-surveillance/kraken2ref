# KRAKEN2REF  

This python package identifies suitable reference genomes as well as reads that should be associated with that reference genome in downstream analysis, using the kraken2 taxonomic report as input. Briefly, it first finds all graphs in the report that are rooted at species level ("S"); then analyses the leaf nodes of these graphs to identify one or more of them as suitable reference genomes; and finally outputs information about each selected reference to pass on to downstream processes.  

## Installation  

#### With pip  

```shell
pip install kraken2ref@git+ssh://git@gitlab.internal.sanger.ac.uk/malariagen1/misc_utils/kraken2ref.git
```  

#### From Source  

```shell
git clone https://gitlab.internal.sanger.ac.uk/malariagen1/misc_utils/kraken2ref.git
cd kraken2ref
pip install .
```  

Once installed, run as follows:  
```shell
## parse kraken2 report
kraken2ref -s sample_id -i path/to/kraken2/report.txt -o ./ -t min_read_threshold -m kmeans -s decomposed -q

## sort reads by reference (requires parse_report to have been run before)
sort_reads -s sample_id -fq1 path/to/fq1.fq -fq2 path/to/fq2.fq -k path/to/output.kraken -r path/to/kraken2ref.json -m tree
```  

#### From Singularity  
##TODO

#### From Docker  
##TODO

# List of Arguments  
## `kraken2ref`  
- `-v` [switch]: Print version  

- `-s` [str]: Sample ID [REQUIRED]  
- `-i` [path]: (Ideally the absolute) path to kraken2 taxonomy report file [REQUIRED]  
- `-t` [int]: Minimum number of reads assigned to a leaf node for it to be considered [OPTIONAL][Default = 100]  
- `-o` [path]: Path to output directory [OPTIONAL][Default = "./"]  
- `-m` [str]: Polling method to use [OPTIONAL][DEFAULT = "kmeans"]["kmeans", "tiles"]  
- `-x` [str]: Suffix to apply to `sample_id` when creating output JSON file [OPTIONAL][Default = "decomposed"]  
- `-q` [switch]: Whether to log to stderr or not [OPTIONAL][Default = True]  

## `sort_reads`  
- `-s` [str]: Sample ID [REQUIRED]  
- `-fq1` [path]: Path to R1 fastq file [REQUIRED]  
- `-fq2` [path]: Path to R2 fastq file [REQUIRED]  
- `-k` [path]: Path to kraken2 output.kraken file [REQUIRED]  
- `-r` [path]: Path to JSON file produced by `kraken2r parse_report` [OPTIONAL ONLY IF USING `-m unique`]  
- `-o` [path]: Path to output directory [OPTIONAL ONLY IF NOT USING `-r`][Default = "path/to/ref_json"]  
- `-m` [str]: Specify sorting mode [OPTIONAL][DEFAULT = "unique"]  
- `-u` [switch]: Whether to update the JSON file produced by `kraken2r parse_report` inplace or produce a new, updated copy [OPTIONAL][Default: produce new][ONLY USED IF `-r` SPECIFIED]  
- `-c` [switch]: Whether to dump all reads for a species into one file-pair (as opposed to producing a file-pair _per reference_)[ONLY USED IF USING `-m tree`]  

### Kraken2 Taxonomy Report  

Each line in the taxonomy report contains kraken2 output information for a single taxon; the information is presented with the following columns:  
1. **% of Reads Assigned**: The percentage of total reads assigned up to that taxon level  
2. **Cumulative #Reads Assigned**: The number of reads assigned up to that taxon level  
3. **#Reads Assigned**: The number of reads directly assigned to that taxon level  
4. **Taxon Level**: The short descriptor of taxonomic level (eg. "G" for Genus, "S" for Species, "S1" for sub-species or equivalent, etc)  
5. **Taxon ID**: The unique ID for that taxon  
6. **Descriptive Name**: The descriptive name of that taxon (eg. "Orthomyxoviridae", "Severe acute respiratory syndrome coronavirus 2")  

> If running kraken2 with `--report-minimizer-info`, two columns are added to between "#Reads Assigned" and "Taxon Level". There are: "Total #Minimizers Assigned to Taxon" and "#Unique Minimizers Assigned to Taxon".  

### Summary of the Kraken2 Read Assignment Algorithm  

Kraken2 depends on a database; the database contains a representation of the taxonomic tree/graph, and nucleotide sequences of the nodes in this tree/graph where available. Typically, sequence information is available for lower taxonomic levels, often only for what can be called leaf nodes in the taxonomic tree/graph; i.e. there is no single reference sequence for "Alphainfluenzavirus" but there do exist sequences for "Influenza A/H1N1/Isolate X".  

Briefly, once Kraken2 has a database set up, it attempts to assign each read in the input dataset to one of the taxonomic levels, starting at the leaf nodes. If a particular read could be assigned to more than one leaf node, then kraken2 instead assigns the read to the parent node of those leaf nodes. (Note that the parent node may or may not have actual sequence information associated with it in the kraken2 database.) This information is then summarised in the kraken2 taxonomy report as described above.  

For quick access to the contents of the kraken2 taxonomy report, we store salient information from it in a dictionary right at the start, where the keys are indexed nodes, and the value is a tuple containing that node's #Reads Directly Assigned and Taxonomic ID, before proceeding with the remainder of the process.  

##

In the context of viral data analysis pipelines, the root of the taxonomy tree/graph described in the kraken report will, of course, be the Domain "Viruses", with subtaxa following below. Since we expect kraken2 to virtually never encounter reads that can **only** be confidently assigned at very high taxonomic levels (any level higher than "Species", for example; see [here](#summary-of-the-kraken2-read-assignment-algorithm) and [here](https://github.com/DerrickWood/kraken2) for more), we start by "trimming" the taxonomy tree/graph to create one or more subtrees/subgraphs (in green, Fig. 1), each rooted at "Species" level and each extending to the leaf nodes in that branch of the graph.  

| ![Fig.1](assets/kraken_tax_example.png) |
|:--:|
| *Figure 1: Example kraken2 Taxonomy* |  

In the example, we retain a subtree/subgraph where the nodes can be represented as a list:  

```python
nodes = ["S", "S1", "S2", "S3", "S3", "S2", "S3", "S3"]
```

Next, we parse this list of nodes to obtain the two subtrees/subgraphs in this one:  
```python
subgraph1 = ["S", "S1", "S2", "S3", "S3"]
subgraph2 = ["S", "S1", "S2", "S3", "S3"]
```  

Noting that the two look the same, the simple fix is to work with a list of *indexed* nodes, then find paths in this graph:  

```python
nodes = ["S", "S1", "S2", "S3", "S3", "S2", "S3", "S3"]
indexed_nodes = nodes = [(0,"S"), (1,"S1"), (2,"S2"), (3,"S3"), (4,"S3"), (5,"S2"), (6,"S3"), (7,"S3")]

## 0 --> 1 --> 2 --> (3,4)
subgraph1 = [(0,"S"), (1,"S1"), (2,"S2"), (3,"S3"), (4,"S3")]
paths1 = [[(0,"S"), (1,"S1"), (2,"S2"), (3,"S3")],
            [(0,"S"), (1,"S1"), (2,"S2"), (4,"S3")]]

## 0 --> 1 --> 5 --> (6,7)
subgraph2 = [(0,"S"), (1,"S1"), (5,"S2"), (6,"S3"), (7,"S3")]
paths2 = [[(0,"S"), (1,"S1"), (5,"S2"), (6,"S3")],
            [(0,"S"), (1,"S1"), (5,"S2"), (7,"S3")]]
```

##
> Note that the indices used to construct the data dictionary and by extension the nodes/graphs, correspond to the line indices in the kraken2 report -- this lets us use the outputs of this package to query the kraken2 report directly if ever we need to
##

Now, we can evaluate the leaf nodes of each subgraph separately by checking the number of reads assigned to each of these leaf nodes by kraken2, using the data dictionary. Not all leaf nodes will pass the threshold number set by the user, and those paths through the graph will be discarded. If the entire graph contains insufficient reads, it is disregarged. However, we focus on checking firs the leaf nodes, and then their parent nodes. Any higher, and you usually risk evaluating at the species level or similar -- this is not bad, _per se_, but does mean that we then expect lots of data duplication when it is time to sort reads by the chosen reference.   

So, in case no leaf nodes in a subgraph pass the threshold, the algorithm jumps up one taxonomic level; for example, in the graph `[(0,"S"), (1,"S1"), (5,"S2"), (6,"S3"), (7,"S3")]`, if neither `(6,"S3")` nor `(7,"S3")` have more than the threshold number of reads directly assigned, the output will identify this and note the parent, in this case `(5,"S2")`, as the stopping point. But it will only include `(5,"S2")` in the output if the **cumulative** number of reads assigned at `(5,"S2")`, (i.e. all reads at `(5,"S2")` and below) passes the threshold.  

In most cases, there is expected to be at least one leaf node which passes the threshold; for example, in the subgraph `[(0,"S"), (1,"S1"), (2,"S2"), (3,"S3"), (4,"S3")]`, let us say leaf `(4,"S3")` passes, giving us a valid path `[(0,"S"), (1,"S1"), (2,"S2"), (4,"S3")]` through this subgraph. At this point, the output notes `(4,"S3")` as the chosen reference, records the path to that node, and also records the taxonomic IDs of **ALL** nodes in the entire parent graph `[(0,"S"), (1,"S1"), (2,"S2"), (3,"S3"), (4,"S3")]` -- this allows us to retain all read information associated with this parent graph, and potentially use all those reads to align to/call consensus on/analyse with the chose reference.  

## Polling  

### KMeans-Based Outlier Detection  

In this approach to outlier analysis, we use `sklearn.cluster`, specifically the `KMeans` module. Briefly, we conceptualise the list of references and their correspoding number of assigned reads as a frequency distribution. This frequency distribution is reshaped to a 2D `numpy` array so that `KMeans` can use it. Then, we "cluster" this distribution using `KMeans`, sort the frequenccies by their distance from the cluster centroid, and retain those outlier frequencies that are to the right of the median (i.e. those that are big numbers, rather than those that are outliers because they are small).  

| ![Fig.2](assets/polling_kmeans.png) |
|:--:|
| *Figure 2: KMeans-Based Outlier Analysis* |  

### Quantile-Based Outlier Detection  

In this approach, we use quartiles 1 and 3 (Q1 and Q3 respectively) to calculate the interquartile range (IQR) of the distribution, then we use these values to set up quantile "fences"; a "left fence": `Q1 - (1.5 x IQR)` and a "right fence": `Q3 + (1.5 x IQR)`. Finally, we retain frequences that lie beyond the "right fence".  

| ![Fig.3](assets/polling_tiles.png) |
|:--:|
| *Figure 3: Quantile-Based Outlier Analysis* |  




