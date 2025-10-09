# UV_damage_TFBS_analysis
Code and scripts for the paper "High-throughput characterization of TFs that modulate UV damage formation and repair at single-nucleotide resolution"

### Data used
* Bedfile containing genomic regions for the analysis. For e.g. all intergenic, open chromatin regions (ATAC-seq) for the cell-line of interest with blacklist regions substracted.
* Binding site calls of TF motif clusters from [Vierstra et al. (2020)](https://www.vierstra.org/resources/motif_clustering) intersected with previously described genomic regions.
* Manifest of TF motif clusters to use for the analysis. 
* CPD-seq v2.0 data from [Duan et al. (2024)](https://www.pnas.org/doi/10.1073/pnas.2310854121). Download from GEO accession number [GSE235483](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE235483).  
* hg19 reference genome. Download from [UCSC Genome Browser](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/).

## Running Pipelines
This section assumes all the necessary pre-processing described in `preprocessing.rmd` has been completed and this repo has been cloned to your desired directory.
We *strongly* reccomend that you run this pipeline on your HPC using slurm.
If you are interested in a step-by-step explanation of what these scripts are doing, see the provided jupyter notebook files.

**Note:** Input to the pipeline is defined and can be modified in `pipeline_config.bs`.

#### 1. Set up pipeline
Set up pipeline output folder structure:

`bash run_setup.bs`

#### 2. Prepare background regions
These regions that will be used background model development for both damage formation and repair:

`bash run_prep_background.bs config/pipeline_config.bs`

#### 3. Prepare TF binding site windows
Defines, collects kmer information, and aggregates damages in TF binding site windows:

`bash run_prep_TFBS.bs config/pipeline_config.bs <tf_cluster> <tf_len> <tf>`

Alternatively you can batch this for many TFs listed in a manifest file.
See `data/TF_clusters_motifs.csv` as an example. 
With this file included in your `pipeline_config.sh` file, simply run the following: 

`sbatch batch_prep_tfbs.slurm config/pipeline_config.sh`

### Damage formation analysis
Analyzes damage formation in TF binding site windows:

`bash run_analyze_TFBS_damage_formation.bs config/pipeline_config.bs <tf_len> <tf>`

To batch run this:

`sbatch batch_tfbs_damage_formation.slurm config/pipeline_config.sh`

### Damage repair analysis
Analyzes damage repair in TF binding site windows.
This analysis takes two CPD-seq (or equivalent damage mapping assay) timecourse datasets to perform a repair simulation.

#### 1. Define your timecourse datasets
Set the `TIMECOURSE_EXP1` and `TIMECOURSE_EXP2` variables in `pipeline_config.bs` to the desired timepoints for simulation.

#### 2. Prepare bootstrapping repair dictionary 
This is the tetranucleotide, damage timecourse dictionary the simulation will sample from: 

`bash run_make_repair_dictionary.bs config/pipeline_config.bs`

#### 3. Perform repair simulations
By default, this performs a 10,000 bootstrapping simulation per TF binding site position for both strands.
It takes a very long time (GPU implementation TBD). 
**Therefore, you should definitely run this pipeline only on your HPC with slurm:**

`sbatch batch_repair_simulations.slurm config/pipeline_config.bs`

This creates a job for each position in a TF binding site window.

#### 4. Analyze repair simulations
Scales and analyzes repair simulations:

`bash `


