#!/bin/bash
#SBATCH -J repeatviz
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 6:00:00
#SBATCH -o repeatviz.%j.log
# Example SLURM wrapper for clusters that provide tools as apptainer/singularity images (e.g. BioContainers images in one directory).
# On a workstation you don't need this: install the conda env (environment.yml) and call
#   repeatviz run ASM.fa -o OUT -t 8 -m 32 && repeatviz plot OUT -o OUT/fig --edits edits.tsv
#
# usage: sbatch -p <partition> run_repeatviz.slurm.example.sh ASM.fa OUTDIR [plot options...]
#   e.g. sbatch -p <partition> run_repeatviz.slurm.example.sh curated.fasta rv --edits edits.tsv --depth depth.regions.bed.gz
set -euo pipefail
ASM=$1; OUT=$2; shift 2
VENV=${REPEATVIZ_VENV:-$HOME/tools/repeatviz-venv}
BT=${REPEATVIZ_SIF_DIR:?set REPEATVIZ_SIF_DIR to the directory that holds the images (layout <first letter>/<tool>:<version>)}
# make apptainer available first if your site needs it (e.g. module load apptainer)
command -v apptainer >/dev/null || { echo "apptainer not found" >&2; exit 1; }
# bind the file system that holds your data (e.g. REPEATVIZ_BIND=/data); defaults to $HOME
B="-B ${REPEATVIZ_BIND:-$HOME} -B ${TMPDIR:-/tmp}"
export REPEATVIZ_MERYL="apptainer exec $B $BT/m/meryl:1.4.2--hd19868c_0 meryl"
export REPEATVIZ_MERYL_LOOKUP="apptainer exec $B $BT/m/meryl:1.4.2--hd19868c_0 meryl-lookup"
export REPEATVIZ_ULTRA="apptainer exec $B $BT/u/ultra:1.2.1--h9948957_0 ultra"
T=${SLURM_CPUS_PER_TASK:-16}
MEM=$(( ${SLURM_MEM_PER_NODE:-65536} / 1024 ))
"$VENV/bin/repeatviz" run "$ASM" -o "$OUT" -t "$T" -m "$MEM"
"$VENV/bin/repeatviz" plot "$OUT" -o "$OUT/repeatviz" "$@"
