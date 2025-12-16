#!/bin/bash

#SBATCH --partition=rome 
#SBATCH --time=10:00 
#SBATCH --nodes=1 
#SBATCH --tasks-per-node=8
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=12
#SBATCH --mem=100G

INPUT_FILE=$1
BASENAME=$(basename "$INPUT_FILE" .mat)
DIR_NAME="$BASENAME"
ORIGINAL_PATH=$(pwd)

module load gcc/9.4.0-eewq4j6 

# 1. Δημιουργία προσωρινού Virtual Environment στο Scratch
# (Το ονομάζουμε 'my_py_env')
export PY_ENV_DIR="$SCRATCH/my_py_env"

if [ ! -d "$PY_ENV_DIR" ]; then
    echo "Creating virtual environment..."
    python3 -m venv "$PY_ENV_DIR"
    
    # Ενεργοποίηση
    source "$PY_ENV_DIR/bin/activate"
    
    # Εγκατάσταση βιβλιοθηκών (γίνεται μόνο την πρώτη φορά)
    echo "Installing numpy & matplotlib..."
    pip install --upgrade pip
    pip install numpy matplotlib scipy h5py
else
    # Αν υπάρχει ήδη, απλά το ενεργοποιούμε
    source "$PY_ENV_DIR/bin/activate"
fi

cd $SCRATCH
if [ ! -d "$DIR_NAME" ]; then
    echo "Ο φάκελος '$DIR_NAME' δεν υπάρχει. Εκτέλεση μετατροπής..."
    python3 "$ORIGINAL_PATH/convert_mat_to_bin.py" "$INPUT_FILE"
else
    echo "Ο φάκελος '$DIR_NAME' υπάρχει ήδη. Προσπέραση μετατροπής."
fi

cd "$ORIGINAL_PATH"
module purge
module load gcc/9.4.0-eewq4j6 openmpi/4.1.2-akxtzzl 

mpicc -O3 -fopenmp mpi_2.c -o mpi_2

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
srun -n $SLURM_NTASKS ./mpi_2 "$SCRATCH"/"$DIR_NAME"
