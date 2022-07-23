## 1. Prepping guides

### 1.1 zymoBIOMICS data on NCBI

The Sequence Read Archive (SRA) contains many sequencing databases from different projects.
Similar to their aims, we also hope to be able to reproduce the results and discover new finding from data analysis.
Here we are looking at research projects that sequenced the zymoBIOMIC Mock community.

The list of SRA files that are used for the analysis can be found in the `sra_files` directory.

#### 1.1.1 BGISEQ

Name, looks at different extraction methods for BGI sequencing.
Five different extraction methods were tested.


| SampleID   | SRA        | Protocol |
|------------|------------|----------|
| D6300-27   | ERR4097245 | MetaHIT  |
| D6300-28   | ERR4097111 | MetaHIT  |
| D6300-29   | ERR4097243 | MetaHIT  |
| D6300-30   | ERR4097237 | MetaHIT  |
| D6300-31   | ERR4097238 | MetaHIT  |
| D6300-32   | ERR4097276 | MetaHIT  |
| D6300-2-1  | ERR4097261 | MN       |
| D6300-2-2  | ERR4097241 | MN       |
| D6300-2-3  | ERR4097242 | MN       |
| D6300-3-1  | ERR4097244 | MN       |
| D6300-4-1  | ERR4097272 | MN       |
| D6300-4-2  | ERR4097208 | MN       |
| D6300-13   | ERR4097266 | MP       |
| D6300-14   | ERR4097269 | MP       |
| D6300-15   | ERR4097268 | MP       |
| D6300-16   | ERR4097271 | MP       |
| D6300-17   | ERR4097270 | MP       |
| D6300-18   | ERR4097262 | MP       |
| D6300-11   | ERR4097264 | PS       |
| D6300-12   | ERR4097232 | PS       |
| D6300-2-7  | ERR4097239 | PS       |
| D6300-2-8  | ERR4097240 | PS       |
| D6300-7    | ERR4097176 | PS       |
| D6300-9    | ERR4097177 | PS       |
| D6300-37   | ERR4097212 | Q        |
| D6300-38   | ERR4097211 | Q        |
| D6300-39   | ERR4097210 | Q        |
| D6300-40   | ERR4097172 | Q        |
| D6300-41   | ERR4097173 | Q        |
| D6300-42   | ERR4097174 | Q        |
| D6300-4-20 | ERR4097207 | ZYMO     |
| D6300-4-21 | ERR4097206 | ZYMO     |
| D6300-4-22 | ERR4097205 | ZYMO     |
| D6300-4-23 | ERR4097204 | ZYMO     |
| D6300-4-24 | ERR4097171 | ZYMO     |
| D6300-6-19 | ERR4097175 | ZYMO     |


### 1.2 Conda environments

Conda is great for reprduce... as you can install environment with specific tools/programs.
The basic command:

`conda create --name <environment_name> -c <channel> <tool1> <tool2>`

Create an environment called sra-tools from the channel bioconda with the program sra-tools


`conda create --name sra-tools -c bioconda sra-tools`


Activate the environment

`conda activate sra-tools`

To return to your previous environment, deactivate conda

`conda deactivate`

<details>
  <summary>Conda with yaml</summary>

    It is also possible to install conda environment from yaml files in the `envs` folder
    
    `conda env create -f envs/sra-tools.yaml`
    
</details>

