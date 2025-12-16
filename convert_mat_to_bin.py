import sys
import os
import numpy as np
import scipy.io
import h5py

def convert(mat_file):
    # 1. Βρες το όνομα του φακέλου (χωρίς το .mat)
    base_name = os.path.splitext(os.path.basename(mat_file))[0]
    
    # 2. Φτιάξε τον φάκελο (αν δεν υπάρχει)
    print(f"Creating directory: {base_name}/")
    os.makedirs(base_name, exist_ok=True)
    
    print(f"Checking file format for: {mat_file}...")
    
    # --- (Η λογική φόρτωσης μένει ίδια - Scipy vs H5PY) ---
    try:
        mat = scipy.io.loadmat(mat_file)
        print("Detected Version: < v7.3 (Standard)")
        if 'Problem' in mat:
            A = mat['Problem']['A'][0][0]
        else:
            keys = [k for k in mat.keys() if not k.startswith('_')]
            A = mat[keys[0]]
        indptr = A.indptr
        indices = A.indices
        n = A.shape[0]
        nnz = A.nnz
    except NotImplementedError:
        print("Detected Version: v7.3 (HDF5 Based)")
        with h5py.File(mat_file, 'r') as f:
            if 'Problem' in f:
                group = f['Problem']['A']
            else:
                key = list(f.keys())[0]
                group = f[key]
            indptr = np.array(group['jc']).flatten()
            indices = np.array(group['ir']).flatten()
            n = len(indptr) - 1
            nnz = len(indices)

    print(f"Graph loaded! Nodes: {n}, Edges: {nnz}")

    # --- ΑΠΟΘΗΚΕΥΣΗ ΣΤΟΝ ΦΑΚΕΛΟ ---
    
    # Συνθέτουμε τα paths: π.χ. "file1/row_ptr.bin"
    path_row = os.path.join(base_name, "row_ptr.bin")
    path_col = os.path.join(base_name, "col_ind.bin")
    path_meta = os.path.join(base_name, "meta.txt")
    
    print(f"Saving to {path_row} (int64)...")
    indptr.astype(np.int64).tofile(path_row)
    
    print(f"Saving to {path_col} (int32)...")
    indices.astype(np.int32).tofile(path_col)
    
    print(f"Saving to {path_meta}...")
    with open(path_meta, "w") as f:
        f.write(f"{n}\n{nnz}")
        
    print(f"Done! All files are inside the folder: {base_name}/")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python convert_organized.py file.mat")
    else:
        convert(sys.argv[1])