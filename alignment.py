import subprocess
import os
import io
import re
import tempfile
import platform
from Bio import AlignIO
from visualization import format_newick_string

def run_clustalw(fasta_file):
    # Create temporary file for alignment output
    with tempfile.NamedTemporaryFile(suffix='.aln', delete=False) as temp_aln:
        output_file_path = temp_aln.name
    try:
        # ClustalW will create the guide tree file with .dnd extension in the same location as input
        guide_tree_file = os.path.splitext(fasta_file)[0] + '.dnd'

        print(f"ClustalW alignment file created: {output_file_path}")
        print(f"ClustalW guide tree file expected: {guide_tree_file}")

        script_dir = os.path.dirname(os.path.abspath(__file__))
        if platform.system() == 'Windows':
            # On Windows, use the Windows-specific executable
            clustalw2 = "clustalw2.exe"
        elif platform.system() == 'Linux':
            # On Linux, use the Linux-specific executable
            clustalw2 = "clustalw2"
        else:
            raise ValueError("Unsupported operating system")
        
        #clustalw2 = os.path.join(script_dir, r"bin\clustalw2")
        command = f"{clustalw2} -INFILE={fasta_file} -OUTFILE={output_file_path}"
        print(f"Running command: {command}")
    
        subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error: {e}")

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
        # Determine the platform (Windows or Linux)
        script_dir = os.path.dirname(os.path.abspath(__file__))
        if platform.system() == 'Windows':
            # On Windows, use the Windows-specific executable
            muscle_exe = "muscle3.8.31_i86win32.exe"
        elif platform.system() == 'Linux':
            # On Linux, use the Linux-specific executable
            muscle_exe = "muscle3.8.31_i86linux64"
        else:
            raise ValueError("Unsupported operating system")

        print(f"Using MUSCLE executable: {muscle_exe}")
        
        # Run the MUSCLE command
        command = [muscle_exe, '-in', fasta_file, '-out', output_file_path]
        result = subprocess.run(command, capture_output=True, text=True)
        
        if result.returncode == 0 and os.path.isfile(output_file_path):
            with open(output_file_path, "r") as file:
                alignment_data = io.StringIO(file.read())
            alignment = AlignIO.read(alignment_data, "fasta")
            return alignment, output_file_path
        else:
            print("Error running MUSCLE:", result.stderr)
            return None, None
            
    finally:
        # Files will be cleaned up after the app is done using them
        pass

def generate_tree_with_fasttree(alignment_file):
    guide_tree_file = re.sub(r"\.aln$", ".dnd", alignment_file)
    print(f"FastTree guide tree file created: {guide_tree_file}")
    
    # Set path to FastTree executable (adjust as needed)
    fasttree_exe = r"bin\FastTree"  # Replace with full path if FastTree is not in PATH

    command = [fasttree_exe, "-out", guide_tree_file, alignment_file]
    result = subprocess.run(command, capture_output=True, text=True)
    
    if result.returncode == 0 and os.path.isfile(guide_tree_file):
        return guide_tree_file
    else:
        print("Error generating guide tree:", result.stderr)
        return None

def run_alignment(algorithm, fasta_file):
    if algorithm == "ClustalW":
        alignment, output_file, guide_tree_file = run_clustalw(fasta_file)
        if guide_tree_file:
            with open(guide_tree_file, 'r') as f:
                tree_str = f.read()
                formatted_tree = format_newick_string(tree_str)
                # Write the formatted tree back to file
                with open(guide_tree_file, 'w') as f:
                    f.write(formatted_tree)
        return alignment, output_file, guide_tree_file
    elif algorithm == "MUSCLE":
        alignment, output_file = run_muscle(fasta_file)
        if alignment:
            guide_tree_file = generate_tree_with_fasttree(output_file)
            if guide_tree_file:
                with open(guide_tree_file, 'r') as f:
                    tree_str = f.read()
                    formatted_tree = format_newick_string(tree_str)
                    # Write the formatted tree back to file
                    with open(guide_tree_file, 'w') as f:
                        f.write(formatted_tree)
            return alignment, output_file, guide_tree_file
        return alignment, output_file, None
