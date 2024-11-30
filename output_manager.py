# output_manager.py
from Bio import AlignIO
from io import StringIO
from Bio.Align import MultipleSeqAlignment


# def save_aligned_output(output_file, alignment, format_type="clustal"):
#     """Save the aligned output to a file in the specified format."""
#     AlignIO.write(alignment, output_file, format_type)

# def load_aligned_sequences(file_path, algorithm):
#     """Load aligned sequences from a file."""
#     if algorithm == 'ClustalW':
#         alignment = AlignIO.read(file_path, "clustal") 
#     else:
#         alignment = AlignIO.read(file_path, "fasta") 
#     return alignment
def format_alignment_to_clustal_with_and_without_colors(alignment):
    """Format a MultipleSeqAlignment object to CLUSTAL format with and without color coding."""
    
    # Amino acid colors
    color_map = {
        'A': '#FF6347', 'C': '#32CD32', 'D': '#FFD700', 'E': '#4682B4',
        'F': '#FF69B4', 'G': '#ADFF2F', 'H': '#FF8C00', 'I': '#8A2BE2',
        'K': '#00FA9A', 'L': '#DA70D6', 'M': '#8B4513', 'N': '#0000CD',
        'P': '#FFD700', 'Q': '#6A5ACD', 'R': '#20B2AA', 'S': '#D2691E',
        'T': '#D3D3D3', 'V': '#FF4500', 'W': '#FF1493', 'Y': '#7FFF00',
        '-': '#B0C4DE',
    }

    # Define conservation groups
    strong_groups = [
        set("STA"), set("NEQK"), set("NHQK"), set("NDEQ"),
        set("QHRK"), set("MILV"), set("MILF"), set("HY"), set("FYW")
    ]
    weak_groups = [
        set("CSA"), set("ATV"), set("SAG"), set("STNK"), set("STPA"),
        set("SGND"), set("SNDEQK"), set("NDEQHK"), set("NEQHRK"),
        set("FVLIM"), set("HFY")
    ]

    # Updated CSS styles for better spacing
    css_styles = """
    <style>
        .seq-char {
            display: inline-block;
            width: 1.2ch;
            text-align: center;
            margin: 0;
            font-size: 14.5px;
        }
        .conservation-char {
            display: inline-block;
            width: 1.2ch;
            text-align: center;
            margin: 0;
            font-size: 14.5px;
        }
        .sequence-block {
            font-size: 12px;
            font-family: monospace;
        }
        .header {
            font-size: 14px;
            margin-bottom: 10px;
        }
    </style>
    """

    no_color_output = "\n"
    color_output = f"{css_styles}\n"
    
    seq_length = alignment.get_alignment_length()
    line_width = 60
    name_width = 15
    cumulative_positions = {record.id: 0 for record in alignment}

    for start in range(0, seq_length, line_width):
        chunk_columns = []
        for i in range(start, min(start + line_width, seq_length)):
            chunk_columns.append(set(alignment[:, i].replace("-", "").upper()))

        # Process sequences
        for record in alignment:
            stripped_id = record.id.split('|')[0][:name_width]
            sequence_chunk = str(record.seq[start:start + line_width])
            
            sequence_length = len(sequence_chunk.replace('-', ''))
            cumulative_positions[record.id] += sequence_length

            # No color output
            no_color_output += f"{stripped_id:<{name_width}} {sequence_chunk} {cumulative_positions[record.id]:>4}\n"

            # Color output with reduced line height
            color_output += f"<pre style='margin: 0; line-height: 1;'><span style='display: inline-block; width: {name_width}ch;'>{stripped_id}</span> "

            for pos, aa in enumerate(sequence_chunk):
                color = color_map.get(aa.upper(), 'black')
                display_char = '&#8209;' if aa == '-' else aa
                color_output += f"<span class='seq-char' style='color:{color};'>{display_char}</span>"

            color_output += f" {cumulative_positions[record.id]:>4}</pre>"

        # Add conservation symbols
        conservation_line = []
        for i in range(start, min(start + line_width, seq_length)):
            column = alignment[:, i]
            unique_residues = set(column.replace("-", ""))

            if len(unique_residues) == 1:
                conservation_line.append("*")
            elif any(unique_residues <= group for group in strong_groups):
                conservation_line.append(":")
            elif any(unique_residues <= group for group in weak_groups):
                conservation_line.append(".")
            else:
                conservation_line.append(" ")

        conservation_line_str = "".join(conservation_line)
        no_color_output += " " * (name_width + 1) + conservation_line_str + "\n\n"
        
        # Update conservation line style to match
        color_output += f"<pre style='margin: 0; line-height: 1;'><span style='display: inline-block; width: {name_width}ch;'></span> "
        for symbol in conservation_line_str:
            color_output += f"<span class='conservation-char'>{symbol}</span>"
        color_output += "</pre>\n\n"

    return no_color_output, color_output