import streamlit as st
import os
import tempfile
import time
from alignment import run_alignment
from benchmark import evaluate_alignment
from visualization import calculate_conservation_score, plot_guide_tree, plot_plotly_heatmap
from output_manager import format_alignment_to_clustal_with_and_without_colors

st.title("Multiple Protein Sequence Alignment App")

if 'example_input' not in st.session_state:
    st.session_state.example_input = ""
if 'input_file_name' not in st.session_state:
    st.session_state.input_file_name = ""
if 'sequence_input' not in st.session_state:
    st.session_state.sequence_input = ""

if st.button("Use Example"):
    example_file_path = r"temp_sequences.fasta" 
    with open(example_file_path, "r") as example_file:
        st.session_state.example_input = example_file.read()
    st.session_state.sequence_input = st.session_state.example_input 

with st.form("alignment_form"):
    uploaded_file = st.file_uploader("Upload your FASTA file", type=["fasta", "fa"])
    
    # Text area uses the value from session state for persistence
    sequence_input = st.text_area(
        "Or type your protein sequences here (in FASTA format):", 
        height=150, 
        value=st.session_state.sequence_input
    )
    
    algorithm = st.radio(
        "Select Alignment Algorithm", 
        ["ClustalW", "MUSCLE"],
        help='ClustalW: Reliable but slower for large datasets.\n\nMUSCLE: Faster and more accurate, ideal for large sequences.'
    )

    submit_button = st.form_submit_button("Run Alignment")

# After form submission, process alignment
if submit_button:
    st.session_state["overall_start_time"] = time.time()
    
    # Store text input back to session state to ensure persistence
    st.session_state.sequence_input = sequence_input

    # Save sequence data to a temporary file
    with tempfile.NamedTemporaryFile(delete=False, suffix=".fa", mode="w") as tmp_file:
        if uploaded_file is not None:
            content = uploaded_file.read().decode("utf-8")
            st.session_state.input_file_name = uploaded_file.name  
        elif sequence_input:  
            content = sequence_input
            st.session_state.input_file_name = "input_sequences.fasta"  
        else:
            st.warning("Please upload a file or enter sequences.")
            st.stop()
        
        # Validate FASTA format
        content = content.strip()  
        if not content.startswith('>'):
            st.error("Invalid FASTA format. Sequences must start with '>'")
            st.stop()
            
        tmp_file.write(content)
        input_file_path = tmp_file.name

    sequences = content.strip().split(">")
    sequence_count = len([seq for seq in sequences if seq.strip()])
    if sequence_count <= 1:
        st.error("Only one sequence found. Please provide multiple sequences for alignment.")
        st.stop()

    alignment, output_file, guide_tree_file = run_alignment(algorithm, input_file_path)
    
    if alignment:
        st.session_state["alignment"] = alignment
        st.session_state["output_file"] = output_file
        st.session_state["guide_tree_file"] = guide_tree_file
        st.session_state["formatted_alignment"] = format_alignment_to_clustal_with_and_without_colors(alignment)
        st.session_state["input_file_path"] = input_file_path 
        st.session_state["conservation_scores"] = calculate_conservation_score(alignment)

if "alignment" in st.session_state:
    alignment = st.session_state["alignment"]
    output_file = st.session_state["output_file"]
    guide_tree_file = st.session_state["guide_tree_file"]
    no_color, color = st.session_state["formatted_alignment"]
    input_file_path = st.session_state["input_file_path"]
    input_file_name = st.session_state.input_file_name  
    conservation_scores = st.session_state["conservation_scores"]

    tabs = st.tabs(["Tool Output", "Conserved Regions", "Dendrogram", "Evaluation Results", "Result Files"])

    with tabs[0]:
        st.header("Tool Output")
        st.markdown("### Plain Text (No Color)")
        st.markdown(f"```text\n{no_color}\n```", unsafe_allow_html=True)
        st.markdown("### Color-Coded Alignment")
        st.markdown(color, unsafe_allow_html=True)

    # Conserved Regions Tab
    with tabs[1]:
        st.header("Conserved Regions")
        evaluation_results = evaluate_alignment(algorithm, input_file_path)

        # Identity Matrix Section
        st.write("### Sequence Identity Matrix (%):")
        matrix = evaluation_results['identity_matrix']
        sequence_names = evaluation_results['sequence_names']
        

        
        for i, row in enumerate(matrix):
            row_str = f"{sequence_names[i][:10]:<10}"
            for val in row:
                row_str += f"{val:>10.1f}"
            st.text(row_str)
        
        st.write("### Conservation Score Heatmap:")
        plotly_fig = plot_plotly_heatmap(conservation_scores)
        st.plotly_chart(plotly_fig)

    with tabs[2]:
        st.header("Dendrogram")
        if guide_tree_file and os.path.exists(guide_tree_file):
            with open(guide_tree_file, "r") as file:
                guide_tree_content = file.read()
            st.text("Phylogenetic Tree:")
            st.markdown(f"```text\n{guide_tree_content}\n```", unsafe_allow_html=True)
            guide_tree_fig = plot_guide_tree(guide_tree_file)
            if guide_tree_fig:
                st.text("Guide Tree:")
                st.pyplot(guide_tree_fig)
        else:
            st.write("Guide tree file not found.")

    # Evaluation Results Tab
    with tabs[3]:
        overall_elapsed_time = (time.time() - st.session_state["overall_start_time"]) * 1000 
        st.write("### Evaluation Results:")
        st.write(f"**Algorithm Used:** {evaluation_results['algorithm']}")
        st.write(f"**Overall Time (ms):** {overall_elapsed_time:.2f} ms")
        st.write(f"**Total Sequences:** {evaluation_results['total_sequences']}")
        st.write(f"**Total Length:** {evaluation_results['total_length']} characters")
        # st.write(f"**Gap Count:** {evaluation_results['gap_count']}")

    # Result Files Tab
    with tabs[4]:  
        st.header("Result Files")

        if algorithm == "ClustalW":
            aligned_sequences_file = "aligned_clustal.aln"
            guide_tree_file_name = "guidetree_clustal.dnd"
        else: 
            aligned_sequences_file = "aligned_muscle.aln"
            guide_tree_file_name = "guidetree_muscle.dnd"

        matrix_content = "Sequence Identity Matrix (%)\n\n"
        # Add header
        matrix_content += "Sequence" + " " * 10
        for i in range(len(sequence_names)):
            matrix_content += f"{i+1:>10}"
        matrix_content += "\n"
        # Add rows
        for i, row in enumerate(matrix):
            row_str = f"{sequence_names[i][:10]:<10}"
            for val in row:
                row_str += f"{val:>10.1f}"
            matrix_content += row_str + "\n"

        output_files = [
            ("Input File", input_file_name),
            ("Aligned Sequences", aligned_sequences_file),
            ("Identity Matrix", "identity_matrix.txt"),
        ]

        if guide_tree_file and os.path.exists(guide_tree_file):
            output_files.append(("Guide Tree", guide_tree_file_name))

        for description, file_name in output_files:
            col1, col2, col3 = st.columns([2, 4, 2])  

            with col1:
                st.write(description)
            with col2:
                st.markdown(f"<div style='text-align: center'>{file_name}</div>", unsafe_allow_html=True)
            with col3:
                if description == "Input File":
                    data = open(input_file_path, "rb").read()
                    mime = "text/plain"
                elif description == "Aligned Sequences":
                    data = no_color  
                    mime = "text/plain"
                elif description == "Guide Tree":
                    data = open(guide_tree_file, "rb").read()
                    mime = "text/plain"
                elif description == "Identity Matrix":
                    data = matrix_content
                    mime = "text/plain"

                st.download_button(
                    label="Download",
                    data=data,
                    file_name=os.path.basename(file_name),
                    mime=mime
                )

        if guide_tree_file and not os.path.exists(guide_tree_file):
            st.write("Guide tree file not found.")
