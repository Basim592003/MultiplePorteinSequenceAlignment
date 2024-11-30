import subprocess
import os
import tempfile
from Bio import AlignIO
import io
import re
from tkinter import messagebox
from visualization import format_newick_string

def run_alignment(algorithm, fasta_file):
    if algorithm == "ClustalW":
        # Run ClustalW alignment
        return run_clustalw(fasta_file)
    elif algorithm == "MUSCLE":
        # Run MUSCLE alignment
        return run_muscle(fasta_file)
    else:
        messagebox.showerror("Error", "Invalid algorithm selected.")
        return None, None, None


def run_clustalw(fasta_file):
    # Create temporary file for alignment output
    with tempfile.NamedTemporaryFile(suffix='.aln', delete=False) as temp_aln:
        output_file_path = temp_aln.name

    # ClustalW will create the guide tree file with .dnd extension in the same location as input
    guide_tree_file = os.path.splitext(fasta_file)[0] + '.dnd'

    print(f"ClustalW alignment file created: {output_file_path}")
    print(f"ClustalW guide tree file expected: {guide_tree_file}")

    script_dir = os.path.dirname(os.path.abspath(__file__))
    clustalw2 = os.path.join(script_dir, r"bin\clustalw2")
    command = f"{clustalw2} -INFILE={fasta_file} -OUTFILE={output_file_path}"
    print(f"Running command: {command}")
    
    try:
        subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        messagebox.showerror("Error", f"Error running ClustalW: {str(e)}")
        return None, None, None
    
    # Check if the output file was created and handle it
    if os.path.exists(output_file_path):
        with open(output_file_path, "r") as file:
            alignment_data = io.StringIO(file.read())
        
        alignment = AlignIO.read(alignment_data, "clustal")

        # Check if guide tree was created and is not empty
        if os.path.exists(guide_tree_file) and os.path.getsize(guide_tree_file) > 0:
            return alignment, output_file_path, guide_tree_file
        else:
            print("Warning: Guide tree file was not created or is empty")
            return alignment, output_file_path, None
    else:
        print("Error: Alignment output file was not created.")
        return None, None, None


def run_muscle(fasta_file):
    # Create temporary files for MUSCLE alignment output
    with tempfile.NamedTemporaryFile(suffix='.aln', delete=True) as temp_aln:
        output_file_path = temp_aln.name

    print(f"MUSCLE alignment file created: {output_file_path}")

    try:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        print("Contents of bin directory:", os.listdir(os.path.join(script_dir, "bin")))
        muscle_exe = os.path.join(script_dir, "bin", "muscle3.8.31_i86linux64")
        command = [muscle_exe, '-in', fasta_file, '-out', output_file_path]
        result = subprocess.run(command, capture_output=True, text=True)
        
        if result.returncode == 0 and os.path.isfile(output_file_path):
            with open(output_file_path, "r") as file:
                alignment_data = io.StringIO(file.read())
            alignment = AlignIO.read(alignment_data, "fasta")
            guide_tree_file = generate_tree_with_fasttree(output_file_path)
            if guide_tree_file:
                with open(guide_tree_file, 'r') as f:
                    tree_str = f.read()
                    formatted_tree = format_newick_string(tree_str)
                    # Write the formatted tree back to file
                    with open(guide_tree_file, 'w') as f:
                        f.write(formatted_tree)
            return alignment, output_file_path, guide_tree_file
        else:
            messagebox.showerror("Error", f"Error running MUSCLE: {result.stderr}")
            return None, None, None
    except PermissionError:
        messagebox.showerror("Permission Error", f"Permission denied when trying to write to {output_file_path}")
        return None, None, None
    except Exception as e:
        messagebox.showerror("Error", f"An unexpected error occurred: {str(e)}")
        return None, None, None


def generate_tree_with_fasttree(alignment_file):
    guide_tree_file = re.sub(r"\.aln$", ".dnd", alignment_file)
    print(f"FastTree guide tree file created: {guide_tree_file}")
    
    fasttree_exe = r"bin\FastTree"  # Replace with full path if FastTree is not in PATH

    command = [fasttree_exe, "-out", guide_tree_file, alignment_file]
    result = subprocess.run(command, capture_output=True, text=True)
    
    if result.returncode == 0 and os.path.isfile(guide_tree_file):
        return guide_tree_file
    else:
        print("Error generating guide tree:", result.stderr)
        return None
