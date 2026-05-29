# Australian-Structural-Biology-Computing/wisps: Usage

## Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 3 columns, and a header row as shown in the examples below.

```bash
--input '[path to samplesheet file]'
```

### Full samplesheet

An example of the sample sheet used for the pipeline is shown below:

```csv title="samplesheet.csv"
id,sequence,group,type
S1,S1.fasta,A,protein
S2,Sample2.fasta,B,protein
S3,S3.fasta,A,rna
```

| Column     | Description                                                                                                                                                                                                                                                                                                                                        |
| ---------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `id`       | The `id` identifiers have to be unique in the sample sheet and **must not** have an underscore or hyphen (`_`, `-`).                                                                                                                                                                                                                               |
| `sequence` | Full path to FASTA file for sample sequence. File has to have the extension ".fasta", or ".fa"                                                                                                                                                                                                                                                     |
| `type`     | The sequence type, can be either `protein`, `rna`, `smiles`. It is an optional column; by default, it will be considered as `protein`.                                                                                                                                                                                                             |
| `group`    | It is a way to group samples together to create custom interactions. It is an optional column and **must not** have a hyphen (`-`). If not provided, the interactions will be `all-all`, which means all samples against all samples, or you can use `group.a-group.b` to only consider the samples from `group.a` against samples from `group.b`. |

A sequence file can contain multiple entities but if the entities are different types, they must adhere to the colabfold standard described below.

An [example samplesheet](../assets/samplesheet.csv) has been provided with the pipeline.

## Creating Interactions from sample sheet

The interactions generated from the sample sheet are controlled with `--mode`. There are three supported modes:

1. `--mode all-all` (default): the workflow generates all unique sample pairs (order-independent) across the full sample sheet. This mode does not require the `group` column. You can limit pair generation with `--interaction_neighbours <n>`. The default is `0` (no limit). For any positive `n`, only pairs within `n` neighbours in the sample sheet order are kept. For example, with `--interaction_neighbours 2`, each sample is paired with itself and the next two samples in the sample sheet.

2. Group-based mode (for example, `--mode A-A,A-B`): interactions are generated only for the requested group combinations from the `group` column (comma-separated, no spaces). With groups `A` and `B`, `--mode A-A,A-B` will include all `A-A` pairs plus all `A-B` pairs.

3. `--mode manual`: no interaction pairing is generated from the sample sheet. Instead, each FASTA record is treated as one interaction input directly. In this mode, provide interaction FASTA entries in the required ColabFold format:
   - One FASTA record per interaction (`>interaction_id` then one sequence line).
   - Use `:` to separate entities/chains in a complex.
   - Use `molecule_type|sequence` for non-protein entities (for example `dna|ATCG`, `rna|AUGC`, `ccd|ATP`, `smiles|...`; optional copies like `ccd|ATP|2`).
   - If SMILES contains aromatic bonds, replace `:` with `;` as required by ColabFold parsing.

Reference: ColabFold input format examples in the official repository README: <https://github.com/sokrypton/ColabFold>.

## Pipeline parameters

### General Help Options

- `--help` `[boolean, string]`
  Show help for top-level parameters. If a parameter is provided, shows detailed help for that parameter.

- `--help_full` `[boolean]`
  Show help for all non-hidden parameters.

- `--show_hidden` `[boolean]`
  Show hidden parameters (must be used with `--help` or `--help_full`).

---

### Input / Output Options

- `--input` `[string]`
  Path to CSV file describing samples.

- `--outdir` `[string]`
  Output directory for results (must be an absolute path in cloud environments).

- `--mode` `[string]`
  Interaction mode:

  - `all-all`
  - `group-group` (uses `group` column in samplesheet)
    Multiple group-group combinations allowed (comma-separated, no spaces).
    **Default:** `all-all`

- `--tools` `[string]`
  Models to run (comma-separated, no spaces):
  `alphafold3`, `colabfold`, `boltz`
  **Default:** `boltz,colabfold`

- `--analysis_batch_size` `[integer]`
  Number of samples processed in one batch by boltz and Alphafold3.
  **Default:** `20`

- `--colabfold_batch_size` `[integer]`
  Number of samples processed in one batch by ColabFold.
  **Default:** `20`

- `--interaction_neighbours` `[integer]`
  Local neighbourhood size for all-by-all mode.
  **Default:** `0` which means no limit and all samples will be considered.

- `--pool` `[boolean]`
  Enable interaction pooling mode (uses `--mode` group pairs).
  **Default:** `false`

- `--pool_size` `[integer]`
  Max total sequence length per pool (required if pooling).
  **Default:** `2000`

- `--mmseqs_gpu` `[boolean]`
  Run mmseqs on GPU.

- `--use_spire_db` `[boolean]`
  Use Spire database.

---

## Database Options

- `--db` `[string]`
  Path to reference data and model parameters.

- `--colabfold_uniref30` `[string]`
  UniRef30 database.

- `--colabfold_envdb` `[string]`
  ColabFold environment database.

- `--colabfold_uniref30_prefix` `[string]`
  Default: `colabfold_uniref30/uniref30_2302_db`

- `--colabfold_envdb_prefix` `[string]`
  Default: `colabfold_envdb/colabfold_envdb_202108_db`

- `--spire_db` `[string]`
  Spire database.

- `--colabfold_alphafold2_params` `[string]`
  AlphaFold2 parameters for ColabFold.

- `--boltz2_aff` `[string]`
  Boltz affinity file.

- `--boltz2_conf` `[string]`
  Boltz-2 config file.

- `--boltz2_mols` `[string]`
  Boltz-2 molecule files.

- `--alphafold3_params` `[string]`
  AlphaFold3 parameters.

---

## tools Options

- `--colabfold_num_recycles` `[integer]`
  Number of recycles.
  **Default:** `3`

- `--ipsae_pae_cutoff` `[integer]`
  PAE cutoff.
  **Default:** `10`

- `--ipsae_dist_cutoff` `[integer]`
  Distance cutoff.
  **Default:** `10`

---

## Generic Options

- `--version` `[boolean]`
  Show version and exit.

- `--publish_dir_mode` `[string]`
  Output method:
  `symlink`, `rellink`, `link`, `copy`, `copyNoFollow`, `move`
  **Default:** `symlink`

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run Australian-Structural-Biology-Computing/wisps --input ./samplesheet.csv --outdir ./results -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

:::warning
Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).
:::

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run Australian-Structural-Biology-Computing/wisps -profile docker -params-file params.yaml
```

with:

```yaml title="params.yaml"
input: './samplesheet.csv'
outdir: './results/'
<...>
```

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull Australian-Structural-Biology-Computing/wisps
```

### Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [Australian-Structural-Biology-Computing/wisps releases page](https://github.com/Australian-Structural-Biology-Computing/wisps/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducbility, you can use share and re-use [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

:::tip
If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.
:::

## Core Nextflow arguments

:::note
These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).
:::

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

:::info
We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.
:::

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer enviroment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow ` 24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command).
Vii

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
