import dash
from dash import dcc, html, Input, Output, State, dash_table
import dash_bootstrap_components as dbc
import dash_uploader as du
import os
import subprocess
import sys
import pandas as pd
import dash_cytoscape as cyto
from Bio import Phylo
import configparser
from flask import send_from_directory


# --- Load config.ini file settings ---
config = configparser.ConfigParser()
config.read('config.ini')
paths = config['PATHS']

# --- Define Path Variables from the Configuration File ---
JAVA_PATH = paths.get('JAVA', 'java')
TRIMMOMATIC_JAR = paths.get('TRIMMOMATIC_JAR')
TRIMMOMATIC_ADAPTERS = paths.get('TRIMMOMATIC_ADAPTERS_SE')
TRIMMOMATIC_ADAPTERS_PE = paths.get('TRIMMOMATIC_ADAPTERS_PE')
PEAR_PATH = paths.get('PEAR')
FASTQC_PATH = paths.get('FASTQC', 'fastqc')
VSEARCH_PATH = paths.get('VSEARCH', 'vsearch')
MAFFT_PATH = paths.get('MAFFT', 'mafft')
FASTTREE_PATH = paths.get('FASTTREE', 'fasttree')
RAXML_PATH = paths.get('RAXML_NG', 'raxml-ng')

# --- Folder Configuration ---
UPLOAD_DIRECTORY = os.path.join(os.getcwd(), "uploads")
os.makedirs(UPLOAD_DIRECTORY, exist_ok=True)

BOLDIGGER_RESULTS_DIR = os.path.join(UPLOAD_DIRECTORY, "boldigger3_results")
os.makedirs(BOLDIGGER_RESULTS_DIR, exist_ok=True)

app = dash.Dash(
    __name__,
    external_stylesheets=[dbc.themes.BOOTSTRAP],
    suppress_callback_exceptions=True
)
server = app.server
server.config['MAX_CONTENT_LENGTH'] = 10 * 1024 * 1024 * 1024  

du.configure_upload(app, UPLOAD_DIRECTORY)

def read_file_preview(path, n_lines=30):
    if not os.path.exists(path):
        return "File not found."
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as f:
            lines = []
            for i, line in enumerate(f):
                if i >= n_lines:
                    lines.append("... (truncated output) ...")
                    break
                lines.append(line.rstrip())
        return "\n".join(lines)
    except Exception as e:
        return f"Error reading file: {str(e)}"

def fastq_to_fasta(fastq_path, fasta_path):
    if not os.path.exists(fastq_path) or os.path.getsize(fastq_path) == 0:
        return False
    with open(fastq_path, "r") as fq, open(fasta_path, "w") as fa:
        while True:
            header = fq.readline()
            if not header:
                break
            seq = fq.readline()
            plus = fq.readline()
            qual = fq.readline()
            if not (header and seq and plus and qual):
                break
            if header.startswith('@'):
                fa.write(">" + header[1:])
                fa.write(seq)
    return os.path.exists(fasta_path) and os.path.getsize(fasta_path) > 0

def run_vsearch_otus(input_fasta, output_dir, vsearch_path=VSEARCH_PATH, identity=0.97):
    os.makedirs(output_dir, exist_ok=True)
    derep_fasta = os.path.join(output_dir, "derep.fasta")
    otus_fasta = os.path.join(output_dir, "otus.fasta")
    uc_file = os.path.join(output_dir, "clusters.uc")
    derep_cmd = [
        vsearch_path, "--derep_fulllength", input_fasta,
        "--output", derep_fasta,
        "--sizeout"
    ]
    derep_res = subprocess.run(derep_cmd, capture_output=True, text=True)
    if not os.path.exists(derep_fasta) or os.path.getsize(derep_fasta) == 0:
        return derep_res.stdout + derep_res.stderr, "", "", ""
    cluster_cmd = [
        vsearch_path, "--cluster_fast", derep_fasta,
        "--id", str(identity),
        "--centroids", otus_fasta,
        "--uc", uc_file
    ]
    cluster_res = subprocess.run(cluster_cmd, capture_output=True, text=True)
    return derep_res.stdout + derep_res.stderr, cluster_res.stdout + cluster_res.stderr, otus_fasta, uc_file

def run_mafft(input_fasta, output_fasta):
    if not os.path.exists(input_fasta) or os.path.getsize(input_fasta) == 0:
        return False, "Input FASTA for MAFFT not found or empty."
    cmd = [MAFFT_PATH, "--auto", input_fasta]
    with open(output_fasta, "w") as out_f:
        result = subprocess.run(cmd, stdout=out_f, stderr=subprocess.PIPE, text=True)
    if not os.path.exists(output_fasta) or os.path.getsize(output_fasta) == 0:
        return False, "MAFFT failed: output file not created."
    return True, result.stderr

def run_fasttree(input_fasta, output_tree):
    from Bio import SeqIO
    seqs = [rec for rec in SeqIO.parse(input_fasta, "fasta")]
    if len(seqs) < 2:
        return False, "FastTree requires at least 2 valid sequences."
    if not os.path.exists(input_fasta) or os.path.getsize(input_fasta) == 0:
        return False, "Input FASTA for FastTree not found or empty."
    cmd = [FASTTREE_PATH, "-gtr", "-nt", input_fasta]
    with open(output_tree, "w") as out_f:
        result = subprocess.run(cmd, stdout=out_f, stderr=subprocess.PIPE, text=True)
    if not os.path.exists(output_tree) or os.path.getsize(output_tree) == 0:
        return False, "FastTree failed: output file not created."
    return True, result.stderr

def run_raxml(input_fasta, output_dir, n_bootstrap=100):
    if not os.path.exists(input_fasta) or os.path.getsize(input_fasta) == 0:
        return False, "Input FASTA for RAxML-NG not found or empty."
    os.makedirs(output_dir, exist_ok=True)
    prefix = os.path.join(output_dir, "raxml") 
    cmd = [
        RAXML_PATH, "--all",
        "--msa", input_fasta,
        "--model", "GTR+G",
        "--prefix", prefix,
        "--bs-trees", str(n_bootstrap),
        "--redo"
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    best_tree_file_path = prefix + ".raxml.bestTree" 
    if not os.path.exists(best_tree_file_path) or os.path.getsize(best_tree_file_path) == 0:
        return False, result.stdout + "\n" + result.stderr
    return True, result.stdout + "\n" + result.stderr

def uniquify_fasta(in_fasta, out_fasta):
    from Bio import SeqIO
    seen = {}
    with open(out_fasta, "w") as out_handle:
        for i, rec in enumerate(SeqIO.parse(in_fasta, "fasta")):
            base_id = rec.id
            if base_id in seen:
                seen[base_id] += 1
                rec.id = f"{base_id}_{seen[base_id]}"
                rec.description = ""
            else:
                seen[base_id] = 1
            SeqIO.write(rec, out_handle, "fasta")

def filter_empty_aligned_fasta(in_fasta, out_fasta, min_real=10):
    from Bio import SeqIO
    records = []
    for rec in SeqIO.parse(in_fasta, "fasta"):
        real_bases = sum(1 for b in str(rec.seq) if b not in "-.")
        if real_bases >= min_real:
            records.append(rec)
    SeqIO.write(records, out_fasta, "fasta")

def sanitize_ids(in_fasta, out_fasta):
    from Bio import SeqIO
    import re
    with open(out_fasta, "w") as out_handle:
        for rec in SeqIO.parse(in_fasta, "fasta"):
            safe_id = re.sub(r'[^A-Za-z0-9_.-]', '_', rec.id)
            rec.id = safe_id
            rec.name = safe_id
            rec.description = ""
            SeqIO.write(rec, out_handle, "fasta")

def get_uploaded_fastqs(upload_dir):
    fastqs = []
    for root, dirs, files in os.walk(upload_dir):
        for f in files:
            if (f.endswith('.fastq') or f.endswith('.fastq.gz')) and f != "trimmed.fastq":
                file_path = os.path.join(root, f)
                if os.path.getsize(file_path) > 0:
                    fastqs.append(file_path)
    return fastqs

def run_trimmomatic(input_file, output_file, trimmomatic_jar, adapters_file):
    """Executes Trimmomatic for single-end reads."""
    log_message = f"Running Trimmomatic on {os.path.basename(input_file)}..."
    print(log_message)
    cmd = [
        JAVA_PATH, "-jar", trimmomatic_jar, "SE", "-phred33",
        input_file, output_file,
        f"ILLUMINACLIP:{adapters_file}:2:30:10",
        "LEADING:3", "TRAILING:3", "SLIDINGWINDOW:4:15", "MINLEN:36"
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Trimmomatic Error: {result.stderr}")
    return result.stdout, result.stderr

def run_trimmomatic_pe(input_r1, input_r2, output_r1_paired, output_r1_unpaired, output_r2_paired, output_r2_unpaired, trimmomatic_jar, adapters_file):
    """Executes Trimmomatic for paired-end reads."""
    log_message = f"Running Trimmomatic on {os.path.basename(input_r1)} and {os.path.basename(input_r2)}..."
    print(log_message)
    cmd = [
        JAVA_PATH, "-jar", trimmomatic_jar, "PE", "-phred33",
        input_r1, input_r2,
        output_r1_paired, output_r1_unpaired,
        output_r2_paired, output_r2_unpaired,
        f"ILLUMINACLIP:{adapters_file}:2:30:10",
        "LEADING:3", "TRAILING:3", "SLIDINGWINDOW:4:15", "MINLEN:36"
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Trimmomatic PE Error: {result.stderr}")
    return result.stdout, result.stderr

def run_fastqc(input_file, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    cmd = [FASTQC_PATH, input_file, "-o", output_dir]
    result = subprocess.run(cmd, capture_output=True, text=True)
    return result.stdout, result.stderr

def run_pear(forward_file, reverse_file, output_dir, pear_path):
    cmd = [
        pear_path,
        "-f", forward_file,
        "-r", reverse_file,
        "-o", os.path.join(output_dir, "merged")
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    return result.stdout, result.stderr

def biopython_to_cytoscape_elements(tree):
    elements = []
    node_counter = 0

    def _clade_to_elements_recursive(clade, parent_id=None):
        nonlocal node_counter, elements
        node_counter += 1
        
        clade_name_safe = str(clade.name).replace('.', '_').replace(':', '_') if clade.name else f'internal_{node_counter}'
        current_id = f"node_{node_counter}_{clade_name_safe}"

        label = ""
        if clade.name:
            label = clade.name
        elif hasattr(clade, 'confidence') and clade.confidence is not None:
            try:
                label = str(round(float(clade.confidence), 2))
            except (ValueError, TypeError):
                label = str(clade.confidence) 
        
        elements.append({'data': {'id': current_id, 'label': label}})

        if parent_id:
            elements.append({'data': {'source': parent_id, 'target': current_id}})

        if clade.clades: 
            for child_clade in clade.clades:
                _clade_to_elements_recursive(child_clade, current_id)

    if tree and tree.root:
        _clade_to_elements_recursive(tree.root)
    return elements

app.layout = dbc.Container([
    dcc.Store(id='uploaded-files', data={}),
    dcc.Store(id='temp-r1-path-store', data=None), 
    dcc.Store(id='temp-r2-path-store', data=None), 
    dcc.Store(id='temp-se-status-store', data=None), 
    dcc.Store(id='workflow-progress', data={
        'qc_done': False,
        'preprocess_done': False,
        'analysis_done': False,
        'db_done': False,
        'phylo_done': False,
        'results_done': False
    }),
    html.Div([
        html.Img(src="/uploads/logo.png", style={"height": "70px", "marginRight": "8px", "verticalAlign": "middle", "display": "inline-block"}),
        html.H1("PARAFISH", style={"display": "inline-block", "verticalAlign": "middle", "margin": "0"}),
    ], style={"textAlign": "left", "marginLeft": "0", "marginTop": "10px", "marginBottom": "10px"}),
    html.Div([
        html.Label("Choose the step where you want to start:"),
        dcc.Dropdown(
            id="start-step",
            options=[
                {"label": "Quality Control", "value": "tab-qc"},
                {"label": "Preprocessing", "value": "tab-preprocess"},
                {"label": "Sequence Analysis", "value": "tab-analysis"},
                {"label": "Reference Databases", "value": "tab-databases"},
                {"label": "Phylogenetic Analysis", "value": "tab-phylo"},
            ],
            value=None,
            clearable=True,
            style={"width": "350px", "marginBottom": "20px"}
        )
    ]),
    dbc.Row([
        dbc.Col([
            dcc.Tabs(
                id="workflow-tabs",
                value='tab-qc',  
                vertical=True,
                children=[
                    dcc.Tab(label='Quality Control', value='tab-qc', id='tab-qc'),
                    dcc.Tab(label='Preprocessing', value='tab-preprocess', id='tab-preprocess'),
                    dcc.Tab(label='Sequence Analysis', value='tab-analysis', id='tab-analysis'),
                    dcc.Tab(label='Reference Databases', value='tab-databases', id='tab-databases'),
                    dcc.Tab(label='Phylogenetic Analysis', value='tab-phylo', id='tab-phylo'),
                ],
                style={"height": "90vh"}
            ),
        ], width=2, style={"borderRight": "1px solid #ddd", "paddingRight": "0"}),
        dbc.Col([
            html.Div(id='tab-content', className="p-4"),
        ], width=10)
    ], style={"height": "90vh"}),
    html.Div(id='upload-status-qc', style={'display': 'none'}),
    html.Div(id='qc-results', style={'display': 'none'}),
    html.Div(id='preprocess-uploader-container', style={'display': 'none'}),
    html.Div(id='upload-status-preprocess', style={'display': 'none'}),
    html.Div(id='preprocess-results', style={'display': 'none'}),
    html.Div(id='upload-status-analysis', style={'display': 'none'}),
    html.Div(id='analysis-results', style={'display': 'none'}),
    html.Div(id='upload-status-databases', style={'display': 'none'}),
    html.Div(id='db-info', style={'display': 'none'}),
    html.Div(id='boldigger-results', style={'display': 'none'}),
    html.Div(id='upload-status-phylo', style={'display': 'none'}),
    html.Div(id='phylo-results', style={'display': 'none'}),
    html.Div(id='final-results', style={'display': 'none'}), 

], fluid=True)

@app.callback(
    Output('run-pear-btn', 'color'),
    Input('uploaded-files', 'data'),
    prevent_initial_call=False
)
def update_pear_btn_color(uploaded_files):
    mode = (uploaded_files or {}).get('preprocess_mode', 'se')
    if mode == 'pe':
        return "primary"  
    return "secondary"   

@app.callback(
    Output('workflow-tabs', 'value'),
    Input('start-step', 'value'),
    prevent_initial_call=True
)
def go_to_selected_tab(selected_tab):
    return selected_tab

@app.callback(
    [Output('tab-qc', 'disabled'),
     Output('tab-preprocess', 'disabled'),
     Output('tab-analysis', 'disabled'),
     Output('tab-databases', 'disabled'),
     Output('tab-phylo', 'disabled')],
    Input('workflow-progress', 'data')
)
def update_tab_states(progress):
    return [
        False,
        not progress['qc_done'],
        not progress['preprocess_done'],
        not progress['analysis_done'],
        not progress['db_done']
    ]

@app.callback(
    Output('tab-content', 'children'),
    Input('workflow-tabs', 'value'),
    State('start-step', 'value')
)
def render_tab(tab, start_step):
    if not start_step:
        return html.Div([
            html.H3("Welcome to PARAFISH"),
            html.P("PARAFISH is an interactive platform for streamlined analysis of sequencing data from quality control to phylogenetic inference."),
            html.P("Start by selecting the workflow step where your data is ready, and upload your files to proceed."),
        ], style={"textAlign": "center", "marginTop": "10em"})
    elif tab == 'tab-qc':
        upload_component = []
        if start_step == 'tab-qc':
            upload_component = [
                du.Upload(
                    id='du-upload-qc',
                    text='Upload FASTQ para QC',
                    max_files=1,
                    max_file_size=100000, 
                    filetypes=['fastq.gz', 'gz', 'fastq'],
                ),
                html.Div(
                    "Only one file can be analyzed at a time in Quality Control. "
                    "To analyze both, repeat the upload and analysis for the second file.",
                    style={"color": "orange", "marginTop": "1em"}
                ),
                html.Div(id='upload-status-qc', className="mt-3"),
            ]
        return dbc.Card([
            dbc.CardHeader("Quality Control (FastQC)"),
            dbc.CardBody(
                upload_component + [
                    html.P("Visualize quality metrics: read abundance, adapter content, GC distribution, per-base quality."),
                    dbc.Button("Run Quality Control", id="run-qc-btn", color="primary"),
                    dcc.Loading(
                        id="qc-loading",
                        type="default",
                        children=html.Div(id="qc-results", className="mt-3")
                    ),
                ]
            )
        ])
    elif tab == 'tab-preprocess':
        upload_component = [
            html.Label("Modo de Preprocessing:"),
            dbc.RadioItems(
                id="preprocess-mode",
                options=[
                    {"label": "Single-end (1 file)", "value": "se"},
                    {"label": "Paired-end (2 files)", "value": "pe"},
                ],
                value="se", 
                inline=True,
                className="mb-2"
            ),
            dbc.Button("Confirm mode", id="preprocess-confirm-btn", color="info", className="mb-4"),
            html.Div(id="preprocess-uploader-container"), 
            html.Div(id='upload-status-preprocess', className="mt-2"), 
        ]
        return dbc.Card([
            dbc.CardHeader("Preprocessing (Trimmomatic & PEAR)"),
            dbc.CardBody(
                upload_component + [
                    html.P("Trim low-quality regions and adapters (Trimmomatic). Merge paired-end reads (PEAR)."),
                    dbc.Button("Run Trimmomatic", id="run-trim-btn", color="primary", className="me-2"),
                    dbc.Button("Run PEAR", id="run-pear-btn", className="ms-2"),
                    dcc.Loading(
                        id="preprocess-loading",
                        type="default",
                        children=html.Div(id="preprocess-results", className="mt-3"),
                    ),
                ]
            )
        ])
    elif tab == 'tab-analysis':
        upload_component = []
        if start_step == 'tab-analysis': 
            upload_component = [
                du.Upload(
                    id='du-upload-analysis',
                    text='Upload FASTA for Analysis (if starting here)', 
                    max_files=1,                                      
                    max_file_size=100000, 
                    filetypes=['fasta', 'fa'],                        
                ),
                html.Div(id='upload-status-analysis', className="mt-3"),
            ]
        return dbc.Card([
            dbc.CardHeader("Sequence Analysis"),
            dbc.CardBody(
                upload_component + [ 
                    html.P("Choose analysis method:"),
                    dbc.RadioItems(
                        options=[
                            {"label": "OTU Clustering (VSEARCH)", "value": "vsearch"},
                            {"label": "ASV Inference (DADA2)", "value": "dada2"},
                            {"label": "Direct Taxonomic Classification (Kraken2)", "value": "kraken2"},
                        ],
                        value="vsearch",
                        id="analysis-method",
                        inline=True,
                    ),
                    html.Div([
                        html.Label("Identity threshold (OTU clustering):"),
                        dcc.Slider(
                            id="vsearch-identity",
                            min=0.8, max=1.0, step=0.01, value=0.97,
                            marks={0.8: "80%", 0.85: "85%", 0.9: "90%", 0.95: "95%", 0.97: "97%", 1.0: "100%"},
                            tooltip={"placement": "bottom", "always_visible": True}
                        ),
                    ], style={"marginTop": "1em", "marginBottom": "1em"}),
                    dbc.Button("Run Analysis", id="run-analysis-btn", color="primary", className="mt-2"),
                    dcc.Loading(
                        id="analysis-loading",
                        type="default",
                        children=html.Div(id="analysis-results", className="mt-3"),
                    ),
                ]
            )
        ])
    elif tab == 'tab-databases':
        upload_component = []
        if start_step == 'tab-databases':
            upload_component = [
                du.Upload(
                    id='du-upload-databases',
                    text='Upload FASTA to BOLDigger (if starting here)', 
                    max_files=1,                                         
                    max_file_size=100000, 
                    filetypes=['fasta', 'fa'],                           
                ),
                html.Div(id='upload-status-databases', className="mt-3"),
            ]
        return dbc.Card([
            dbc.CardHeader("Reference Databases"),
            dbc.CardBody(
                upload_component + [ 
                    html.P("Select reference database(s) for classification:"),
                    dbc.Checklist(
                        options=[
                            {"label": "NCBI GenBank", "value": "ncbi"},
                            {"label": "SILVA", "value": "silva"},
                            {"label": "PR2", "value": "pr2"},
                            {"label": "BOLD (via BOLDigger)", "value": "bold"},
                        ],
                        value=["bold"], 
                        id="db-selection",
                        inline=True,
                    ),
                    html.Div(id="db-info", className="mt-3"),
                    html.Hr(),
                    html.H5("BOLDigger3 Customization"),
                    html.Div([
                        html.Label("BOLD Database (1-8):"),
                        dcc.Dropdown(
                            id="bold-db-dropdown",
                            options=[
                                {"label": f"{i}: {name}", "value": i} for i, name in [
                                    (1, "ANIMAL LIBRARY (PUBLIC)"),
                                    (2, "ANIMAL SPECIES-LEVEL LIBRARY (PUBLIC + PRIVATE)"),
                                    (3, "ANIMAL LIBRARY (PUBLIC+PRIVATE)"),
                                    (4, "VALIDATED CANADIAN ARTHROPOD LIBRARY"),
                                    (5, "PLANT LIBRARY (PUBLIC)"),
                                    (6, "FUNGI LIBRARY (PUBLIC)"),
                                    (7, "ANIMAL SECONDARY MARKERS (PUBLIC)"),
                                    (8, "VALIDATED ANIMAL RED LIST LIBRARY"),
                                ]
                            ],
                            value=1,
                            clearable=False,
                            style={"width": "350px"}
                        ),
                        html.Label("BOLD Operating Mode (1-3):"),
                        dcc.Dropdown(
                            id="bold-mode-dropdown",
                            options=[
                                {"label": "1: Rapid Species Search", "value": 1},
                                {"label": "2: Genus and Species Search", "value": 2},
                                {"label": "3: Exhaustive Search", "value": 3},
                            ],
                            value=1,
                            clearable=False,
                            style={"width": "350px"}
                        ),
                        html.Label("Custom Thresholds (Species, Genus, Family, Order, Class):"),
                        dcc.Input(id="bold-thresholds", type="text", placeholder="Ex: 99 97 90 85 80", style={"width": "350px", "marginBottom": "5px"}),
                        html.Div(
                            [
                                html.Small(
                                    "Optional. If left blank, BOLDigger's default thresholds will be used. ",
                                    style={"fontStyle": "italic"}
                                ),
                                html.Small(
                                    "Enter values separated by spaces (e.g., 99 97 90 85 80).",
                                    style={"fontStyle": "italic"}
                                )
                            ],
                            style={"color": "grey", "marginBottom": "1em"}
                        ),
                        dbc.Button("Run BOLDigger3 Identification", id="run-boldigger-btn", color="primary", className="mt-2"),
                        dcc.Loading(
                            id="boldigger-loading",
                            type="default",
                            children=html.Div(id="boldigger-results", className="mt-3"),
                        ),
                    ], style={"marginTop": "2em"}),
                ]
            )
        ])
    elif tab == 'tab-phylo':
        upload_component = []
        if start_step == 'tab-phylo': 
            upload_component = [
                du.Upload(
                    id='du-upload-phylo',
                    text='Upload FASTA for Phylogenetic Analysis (if starting here)', 
                    max_files=1,                                                    
                    max_file_size=100000, 
                    filetypes=['fasta', 'fa'], 
                ),
                html.Div(id='upload-status-phylo', className="mt-3"),
            ]
        return dbc.Card([
            dbc.CardHeader("Phylogenetic Analysis"),
            dbc.CardBody(
                upload_component + [ 
                    html.P([
                        "Through the use of MAFFT:",
                        html.Ul([
                            html.Li("Recognised for its precision, speed and high-quality alignments."),
                            html.Li("Essential for robust phylogenetic analyses.")
                        ]),
                        "Construction of Phylogenetic Trees:",
                        html.Ul([
                            html.Li("FastTree 2: Initial topology of the phylogenetic tree (GTR model)."),
                            html.Li("RAxML-NG: Optimises branch length for more accurate and informative trees.")
                        ])
                    ]),
                    html.Label("Number of bootstraps for RAxML-NG:"),
                    dcc.Input(
                        id="raxml-bootstrap-input",
                        type="number",
                        min=10,
                        max=1000,
                        step=10,
                        value=100,
                        style={"width": "120px", "marginBottom": "1em"}
                    ),
                    dbc.Button("Run Phylogenetic Analysis", id="run-phylo-btn", color="primary"),
                    dcc.Loading(
                        id="phylo-loading",
                        type="default",
                        children=html.Div(id="phylo-results", className="mt-3"),
                    ),
                ]
            )
        ])
    elif tab == 'tab-results':
        return dbc.Card([
            dbc.CardHeader("Results & Visualization"),
            dbc.CardBody([
                html.P("Interactive results: tables, plots, phylogenetic trees."),
                dcc.Loading([ 
                    html.Div(id="final-results"),
                ]),
            ])
        ])
    else:
        return html.Div("Invalid tab selected.")

@du.callback(
    output=[
        Output('upload-status-qc', 'children'),
        Output('workflow-progress', 'data', allow_duplicate=True),
        Output('uploaded-files', 'data', allow_duplicate=True)
    ],
    id='du-upload-qc'
)
def get_a_list_qc(filenames):
    if filenames is None:
        return html.Div(), dash.no_update, dash.no_update
        
    if not filenames:
        status_children = html.Div("No files uploaded.")
        return status_children, dash.no_update, dash.no_update
    
    status_children = html.Ul([html.Li(f"Successfully uploaded: {os.path.basename(x)}") for x in filenames])
    progress_update = {
        'qc_done': True,
        'preprocess_done': False,
        'analysis_done': False,
        'db_done': False,
        'phylo_done': False,
        'results_done': False
    }
    uploaded_files_update = {'qc': filenames}
    return status_children, progress_update, uploaded_files_update

@du.callback(
    output=[
        Output('upload-status-analysis', 'children'),
        Output('workflow-progress', 'data', allow_duplicate=True),
        Output('uploaded-files', 'data', allow_duplicate=True)
    ],
    id='du-upload-analysis'
)
def get_a_list_analysis(filenames):
    if filenames is None:
        return html.Div(), dash.no_update, dash.no_update
    
    uploaded_files_update = {}
    status_message = html.Div("No files uploaded for analysis.")
    
    progress_update = {
        'qc_done': True, 
        'preprocess_done': True, 
        'analysis_done': False,
        'db_done': False, 
        'phylo_done': False, 
        'results_done': False
    }

    if filenames and isinstance(filenames, list) and len(filenames) > 0 and filenames[0]:
        actual_filename = filenames[0]
        uploaded_files_update['analysis'] = [actual_filename]
        status_message = html.Ul([html.Li(f"File for Analysis: {os.path.basename(actual_filename)}")])
        progress_update['analysis_done'] = True
        return status_message, progress_update, uploaded_files_update
    else:
        uploaded_files_update['analysis'] = []
        return status_message, dash.no_update, uploaded_files_update

@du.callback(
    output=[
        Output('upload-status-databases', 'children'),
        Output('workflow-progress', 'data', allow_duplicate=True),
        Output('uploaded-files', 'data', allow_duplicate=True)
    ],
    id='du-upload-databases'
)
def get_a_list_databases(filenames):
    if filenames is None:
        return html.Div(), dash.no_update, dash.no_update
        
    uploaded_files_update = {}
    status_message = html.Div("No FASTA files uploaded to BOLDigger.")

    progress_update = {
        'qc_done': True, 
        'preprocess_done': True, 
        'analysis_done': True, 
        'db_done': False,
        'phylo_done': False, 
        'results_done': False
    }

    if filenames and isinstance(filenames, list) and len(filenames) > 0 and filenames[0]:
        actual_filename = filenames[0]
        uploaded_files_update['databases'] = [actual_filename]
        status_message = html.Ul([html.Li(f"FASTA file for BOLDigger: {os.path.basename(actual_filename)}")])
        progress_update['db_done'] = True
        return status_message, progress_update, uploaded_files_update
    else:
        uploaded_files_update['databases'] = []
        return status_message, dash.no_update, uploaded_files_update

@du.callback(
    output=[
        Output('upload-status-phylo', 'children'),
        Output('workflow-progress', 'data', allow_duplicate=True),
        Output('uploaded-files', 'data', allow_duplicate=True)
    ],
    id='du-upload-phylo'
)
def get_a_list_phylo(filenames):
    if filenames is None:
        return html.Div(), dash.no_update, dash.no_update
        
    uploaded_files_update = {}
    status_message = html.Div("No FASTA files uploaded for phylogenetic analysis.")

    progress_update = {
        'qc_done': True, 
        'preprocess_done': True, 
        'analysis_done': True, 
        'db_done': True, 
        'phylo_done': False,
        'results_done': False
    }
    if filenames and isinstance(filenames, list) and len(filenames) > 0 and filenames[0]:
        actual_filename = filenames[0]
        uploaded_files_update['phylo'] = [actual_filename]
        status_message = html.Ul([html.Li(f"FASTA file for phylogeny: {os.path.basename(actual_filename)}")])
        progress_update['phylo_done'] = True
        return status_message, progress_update, uploaded_files_update
    else:
        uploaded_files_update['phylo'] = []
        return status_message, dash.no_update, uploaded_files_update

@du.callback(
    output=[
        Output('temp-se-status-store', 'data'),
        Output('uploaded-files', 'data', allow_duplicate=True)
    ],
    id='du-upload-preprocess-se' 
)
def upload_preprocess_se_to_temp_store(filenames):
    if filenames is None:
        return None, dash.no_update
        
    store_data_main = dash.callback_context.states.get('uploaded-files.data') or {}
    if not isinstance(store_data_main, dict):
        store_data_main = {} 
    
    updated_uploaded_files = store_data_main.copy()
    status_content_for_store = None

    try:
        updated_uploaded_files['preprocess_mode'] = 'se'
        
        current_filenames = []
        if filenames and isinstance(filenames, list):
            current_filenames = [f for f in filenames if f] 

        updated_uploaded_files['preprocess'] = current_filenames
        
        if 'preprocess_r1' in updated_uploaded_files: del updated_uploaded_files['preprocess_r1']
        if 'preprocess_r2' in updated_uploaded_files: del updated_uploaded_files['preprocess_r2']
        
        if current_filenames:
            status_content_for_store = html.Ul([html.Li(f"File Single-End: {os.path.basename(x)}") for x in current_filenames])
        else:
            status_content_for_store = html.Div("No single-end files uploaded or uploads cancelled.")
            
        return status_content_for_store, updated_uploaded_files
    except Exception as e:
        print(f"Error in upload_preprocess_se_to_temp_store: {e}")
        error_status_content = html.Div(f"Error uploading SE: {str(e)}", style={'color': 'red'})
        return error_status_content, store_data_main

@du.callback(
    output=Output('temp-r1-path-store', 'data'), 
    id='du-upload-preprocess-r1', 
)
def upload_preprocess_r1_to_temp_store(filenames):
    if filenames and isinstance(filenames, list) and len(filenames) > 0 and filenames[0]:
        return filenames[0] 
    return None 


@du.callback(
    output=Output('temp-r2-path-store', 'data'), 
    id='du-upload-preprocess-r2', 
)
def upload_preprocess_r2_to_temp_store(filenames):
    if filenames and isinstance(filenames, list) and len(filenames) > 0 and filenames[0]:
        return filenames[0] 
    return None 


def get_first_uploaded_fastq(upload_dir):
    for root, dirs, files in os.walk(upload_dir):
        for f in files:
            if (f.endswith('.fastq') or f.endswith('.fastq.gz')) and f != "trimmed.fastq":
                file_path = os.path.join(root, f)
                if os.path.getsize(file_path) > 0:
                    return file_path
    return None

@app.callback(
    Output('qc-results', 'children'),
    Input('run-qc-btn', 'n_clicks'),
    State('uploaded-files', 'data'),
    prevent_initial_call=True
)
def qc_callback(n_clicks, uploaded_files):
    if not n_clicks:
        return ""
    fastqs = uploaded_files.get('qc', []) if uploaded_files else []
    if not fastqs:
        return html.Div("No valid FASTQ files found for QC. Please upload a non-empty FASTQ file.")
    results = []
    for input_file in fastqs[:2]:  
        
        base_for_subdir = os.path.splitext(os.path.basename(input_file))[0] 
        output_dir_fastqc = os.path.join(UPLOAD_DIRECTORY, "fastqc_results", base_for_subdir)
        os.makedirs(output_dir_fastqc, exist_ok=True)
        
        stdout, stderr = run_fastqc(input_file, output_dir_fastqc)

       
        true_input_basename = os.path.basename(input_file) 
        if true_input_basename.endswith(".fastq.gz"):
            fastqc_report_stem = true_input_basename[:-9] 
        elif true_input_basename.endswith(".fastq"):
            fastqc_report_stem = true_input_basename[:-6]
        elif true_input_basename.endswith(".fq.gz"):
            fastqc_report_stem = true_input_basename[:-6]
        elif true_input_basename.endswith(".fq"):
            fastqc_report_stem = true_input_basename[:-3]
        else:
            fastqc_report_stem = os.path.splitext(true_input_basename)[0]
            if fastqc_report_stem.endswith((".fastq", ".fq")): 
                 fastqc_report_stem = os.path.splitext(fastqc_report_stem)[0]


        actual_report_filename = fastqc_report_stem + "_fastqc.html" 

        
        report_url_path = os.path.join("uploads", "fastqc_results", base_for_subdir, actual_report_filename).replace(os.sep, '/')
        
        download_link = html.A(
            f"Download FastQC Report for {os.path.basename(input_file)}",
            href=f"/{report_url_path}", 
            target="_blank",
            download=actual_report_filename 
        )
        results.append(html.Div([
            html.H5("Downloads"),
            download_link,
            dbc.Accordion([
                dbc.AccordionItem([
                    html.H6("View FastQC HTML Report"),
                    html.Iframe(
                        src=f"/{report_url_path}", 
                        style={"width": "100%", "height": "600px", "border": "1px solid #ccc"}
                    )
                ], title="View FastQC Report"),
                dbc.AccordionItem([
                    html.H6("FastQC stdout"),
                    html.Pre(stdout),
                    html.H6("FastQC stderr"),
                    html.Pre(stderr),
                ], title="Log FastQC"),
            ], start_collapsed=True, always_open=False, style={"marginTop": "1em"})
        ]))
    return html.Div(results)

@app.callback(
    [Output("preprocess-uploader-container", "children"),
     Output('uploaded-files', 'data', allow_duplicate=True),      
     Output('upload-status-preprocess', 'children', allow_duplicate=True), 
     Output('temp-r1-path-store', 'data', allow_duplicate=True),  
     Output('temp-r2-path-store', 'data', allow_duplicate=True),
     Output('temp-se-status-store', 'data', allow_duplicate=True)], 
    Input("preprocess-confirm-btn", "n_clicks"),
    [State("preprocess-mode", "value"), 
     State('uploaded-files', 'data')],
    prevent_initial_call=True
)
def show_preprocess_uploaders_and_reset_state(n_clicks, new_mode, current_uploaded_files_data):
    if not n_clicks: 
        return dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update

    updated_uploaded_files = (current_uploaded_files_data or {}).copy()
    updated_uploaded_files['preprocess_mode'] = new_mode 

    uploader_children_to_render = None
    initial_status_message = ""
    temp_r1_val_to_set = None 
    temp_r2_val_to_set = None 
    temp_se_status_val_to_set = None 

    if new_mode == "se":
        uploader_children_to_render = du.Upload(
            id='du-upload-preprocess-se', 
            text='Upload FASTQ (single-end)', max_files=1, filetypes=['fastq.gz', 'fastq']
        )
        initial_status_message = "Waiting for single-end file..."
        if 'preprocess_r1' in updated_uploaded_files: del updated_uploaded_files['preprocess_r1']
        if 'preprocess_r2' in updated_uploaded_files: del updated_uploaded_files['preprocess_r2']
    
    elif new_mode == "pe":
        uploader_children_to_render = html.Div([
            html.Div(
                "Important: For Paired-End mode, ensure your FASTQ files have the mandatory suffixes '_R1' and '_R2' (e.g., 'my_sample_R1.fastq.gz', 'my_sample_R2.fastq.gz').",
                style={"color": "orange", "marginBottom": "1em", "padding": "10px", "border": "1px solid orange", "borderRadius": "5px"}
            ),
            html.Div([
                html.Div([html.H6("File _R1"), du.Upload(id='du-upload-preprocess-r1', text='Upload _R1.fastq', max_files=1, filetypes=['fastq.gz', 'fastq'])]),
                html.Div([html.H6("File _R2"), du.Upload(id='du-upload-preprocess-r2', text='Upload _R2.fastq', max_files=1, filetypes=['fastq.gz', 'fastq'])])
            ], style={"display": "flex", "gap": "2rem"})
        ])
        initial_status_message = html.Ul([html.Li("File R1: Waiting..."), html.Li("File R2: Waiting...")])
        if 'preprocess' in updated_uploaded_files: del updated_uploaded_files['preprocess']
        if 'preprocess_r1' in updated_uploaded_files: del updated_uploaded_files['preprocess_r1']
        if 'preprocess_r2' in updated_uploaded_files: del updated_uploaded_files['preprocess_r2']
        temp_se_status_val_to_set = None 
    else: 
        return dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update

    return uploader_children_to_render, updated_uploaded_files, initial_status_message, temp_r1_val_to_set, temp_r2_val_to_set, temp_se_status_val_to_set


@app.callback(
    Output('upload-status-preprocess', 'children', allow_duplicate=True),
    Input('temp-se-status-store', 'data'),
    [State('workflow-tabs', 'value'), 
     State('uploaded-files', 'data')], 
    prevent_initial_call=True
)
def update_se_status_from_temp_store(status_data_from_temp_store, active_tab, uploaded_files_data):
    if not status_data_from_temp_store: 
        return dash.no_update

    current_mode = (uploaded_files_data or {}).get('preprocess_mode')
    
    if active_tab == 'tab-preprocess' and current_mode == 'se':
        return status_data_from_temp_store 
    
    return dash.no_update


@app.callback(
    [Output('upload-status-preprocess', 'children', allow_duplicate=True),
     Output('uploaded-files', 'data', allow_duplicate=True)],
    [Input('temp-r1-path-store', 'data'),
     Input('temp-r2-path-store', 'data')],
    [State('uploaded-files', 'data'), 
     State('workflow-tabs', 'value')], 
    prevent_initial_call=True 
)
def consolidate_pe_uploads_and_update_status(r1_path_from_temp, r2_path_from_temp, current_uploaded_files_data, active_tab):
    ctx = dash.callback_context
    if not ctx.triggered_id: 
        return dash.no_update, dash.no_update

    current_mode_from_store = (current_uploaded_files_data or {}).get('preprocess_mode')

    if current_mode_from_store != 'pe':
        return dash.no_update, dash.no_update

    updated_data = (current_uploaded_files_data or {}).copy()

    if ctx.triggered_id == 'temp-r1-path-store':
        if r1_path_from_temp is not None:
            updated_data['preprocess_r1'] = r1_path_from_temp
        else: 
            if 'preprocess_r1' in updated_data: del updated_data['preprocess_r1']

    if ctx.triggered_id == 'temp-r2-path-store':
        if r2_path_from_temp is not None:
            updated_data['preprocess_r2'] = r2_path_from_temp
        else: 
            if 'preprocess_r2' in updated_data: del updated_data['preprocess_r2']

    if 'preprocess' in updated_data:
        del updated_data['preprocess']

    actual_r1_for_status = updated_data.get('preprocess_r1')
    actual_r2_for_status = updated_data.get('preprocess_r2')
    
    status_items = []
    if actual_r1_for_status:
        status_items.append(html.Li(f"File R1: {os.path.basename(actual_r1_for_status)}"))
    else:
        status_items.append(html.Li("File R1: Waiting..."))
    
    if actual_r2_for_status:
        status_items.append(html.Li(f"File R2: {os.path.basename(actual_r2_for_status)}"))
    else:
        status_items.append(html.Li("File R2: Waiting..."))
        
    status_children_content = html.Ul(status_items)
    
    if active_tab == 'tab-preprocess' and updated_data.get('preprocess_mode') == 'pe':
        return status_children_content, updated_data
    elif updated_data.get('preprocess_mode') == 'pe':
        return dash.no_update, updated_data
    else:
        return dash.no_update, dash.no_update

@app.callback(
    Output('preprocess-results', 'children'),
    [Input('run-trim-btn', 'n_clicks'), Input('run-pear-btn', 'n_clicks')],
    State('uploaded-files', 'data'),
    prevent_initial_call=True
)
def preprocess_tools_callback(n_trim, n_pear, uploaded_files):
    ctx = dash.callback_context
    if not ctx.triggered:
        return ""
    button_id = ctx.triggered[0]['prop_id'].split('.')[0]

    if not uploaded_files:
        return html.Div("No files uploaded in the application status.")

    mode = uploaded_files.get('preprocess_mode', 'se')
    r1, r2 = None, None
    
    current_fastqs_for_this_step = []

    if mode == 'se':
        se_files = uploaded_files.get('preprocess', []) 
        if se_files and isinstance(se_files, list) and len(se_files) > 0:
            r1 = se_files[0]
            if r1:
                current_fastqs_for_this_step.append(r1)
    else: 
        r1 = uploaded_files.get('preprocess_r1')
        r2 = uploaded_files.get('preprocess_r2')
        if r1:
            current_fastqs_for_this_step.append(r1)
        if r2:
            current_fastqs_for_this_step.append(r2)
    
    valid_fastqs = [f for f in current_fastqs_for_this_step if f and (f.endswith('.fastq') or f.endswith('.fastq.gz'))]

    print("Uploaded files state:", uploaded_files) 
    print(f"Mode: {mode}, R1: {r1}, R2: {r2}") 
    print("Valid FASTQs for processing:", valid_fastqs) 

    if button_id == "run-trim-btn":
        if mode == 'se':
            if not r1:
                return html.Div("Single-End Mode: No R1 file found. Please upload the file.")
            if r1 not in valid_fastqs:
                 return html.Div(f"Single-End Mode: File R1 ({os.path.basename(r1) if r1 else 'N/A'}) it is not a valid FASTQ.")
        elif mode == 'pe':
            if not r1 or not r2:
                return html.Div("Paired-End Mode: R1 and R2 files are required. Please upload both.")
            if r1 == r2:
                return html.Div("Paired-End Mode: R1 and R2 files cannot be the same.")
            if not (r1 in valid_fastqs and r2 in valid_fastqs):
                return html.Div("Paired-End Mode: One or both files are not valid FASTQ files.")
    elif button_id == "run-pear-btn":
        if mode != 'pe':
            return html.Div("PEAR can only be run in Paired-End mode.")
        if not r1 or not r2:
            return html.Div("PEAR: Files R1 and R2 are required. Please upload both.")
        if r1 == r2:
            return html.Div("PEAR: Files R1 and R2 cannot be the same.")

    if button_id == "run-trim-btn":
        if mode == 'pe' and r1 and r2:
            base1 = os.path.splitext(os.path.basename(r1))[0]
            base2 = os.path.splitext(os.path.basename(r2))[0]
            out_r1_paired = os.path.join(UPLOAD_DIRECTORY, f"{base1}_trimmed_R1_paired.fastq")
            out_r1_unpaired = os.path.join(UPLOAD_DIRECTORY, f"{base1}_trimmed_R1_unpaired.fastq")
            out_r2_paired = os.path.join(UPLOAD_DIRECTORY, f"{base2}_trimmed_R2_paired.fastq")
            out_r2_unpaired = os.path.join(UPLOAD_DIRECTORY, f"{base2}_trimmed_R2_unpaired.fastq")
            
            stdout_trim, stderr_trim = run_trimmomatic_pe(
                r1, r2, out_r1_paired, out_r1_unpaired, out_r2_paired, out_r2_unpaired,
                TRIMMOMATIC_JAR, TRIMMOMATIC_ADAPTERS_PE
            )
            return html.Div([
                html.H5("Downloads Trimmomatic (Paired-End)"),
                html.A(f"Download trimmed R1 (paired)", href=f"/uploads/{os.path.basename(out_r1_paired)}", target="_blank", download=os.path.basename(out_r1_paired)), 
                html.Br(),
                html.A(f"Download trimmed R2 (paired)", href=f"/uploads/{os.path.basename(out_r2_paired)}", target="_blank", download=os.path.basename(out_r2_paired)), 
                dbc.Accordion([
                    dbc.AccordionItem(html.Pre(read_file_preview(out_r1_paired)), title="View R1 Paired Trimmed"),
                    dbc.AccordionItem(html.Pre(read_file_preview(out_r2_paired)), title="View R2 Paired Trimmed"),
                    dbc.AccordionItem([html.H6("Trimmomatic stdout"), html.Pre(stdout_trim), html.H6("Trimmomatic stderr"), html.Pre(stderr_trim)], title="Log Trimmomatic"),
                ], start_collapsed=True, always_open=False, style={"marginTop": "1em"})
            ])
        elif mode == 'se' and r1:
            input_file = r1
            base = os.path.splitext(os.path.basename(input_file))[0]
            output_file_name = f"{base}_trimmed.fastq"
            output_file_path = os.path.join(UPLOAD_DIRECTORY, output_file_name)
            
            stdout_trim, stderr_trim = run_trimmomatic(input_file, output_file_path, TRIMMOMATIC_JAR, TRIMMOMATIC_ADAPTERS)
            return html.Div([
                html.H5("Downloads Trimmomatic (Single-End)"),
                html.A("Download trimmed (single-end)", href=f"/uploads/{output_file_name}", target="_blank", download=output_file_name), 
                dbc.Accordion([
                    dbc.AccordionItem(html.Pre(read_file_preview(output_file_path)), title="View Trimmed FASTQ"),
                    dbc.AccordionItem([html.H6("Trimmomatic stdout"), html.Pre(stdout_trim), html.H6("Trimmomatic stderr"), html.Pre(stderr_trim)], title="Log Trimmomatic"),
                ], start_collapsed=True, always_open=False, style={"marginTop": "1em"})
            ])
        else:
            return html.Div("Unexpected error: Invalid file configuration for Trimmomatic.")

    elif button_id == "run-pear-btn":
        if mode != 'pe': 
            return html.Div("PEAR can only be run in Paired-End mode.")
        
        if not (r1 and r2): 
             return html.Div("PEAR: Original R1 and R2 files not found in the state. Run the upload and Trimmomatic first.")

        base1_orig = os.path.splitext(os.path.basename(r1))[0]
        base2_orig = os.path.splitext(os.path.basename(r2))[0]
        
        trimmed_r1_for_pear = os.path.join(UPLOAD_DIRECTORY, f"{base1_orig}_trimmed_R1_paired.fastq")
        trimmed_r2_for_pear = os.path.join(UPLOAD_DIRECTORY, f"{base2_orig}_trimmed_R2_paired.fastq")

        if not (os.path.exists(trimmed_r1_for_pear) and os.path.exists(trimmed_r2_for_pear)):
            return html.Div(f"PEAR: Trimmed files R1 ({os.path.basename(trimmed_r1_for_pear)}) and/or R2 ({os.path.basename(trimmed_r2_for_pear)}) not found. Run Trimmomatic in paired-end mode first.")
        
        output_dir_pear = os.path.join(UPLOAD_DIRECTORY, "pear_results", f"{base1_orig}_{base2_orig}")
        os.makedirs(output_dir_pear, exist_ok=True)
        
        stdout_pear, stderr_pear = run_pear(trimmed_r1_for_pear, trimmed_r2_for_pear, output_dir_pear, PEAR_PATH)
        
        merged_file_name = "merged.assembled.fastq"
        merged_file_path = os.path.join(output_dir_pear, merged_file_name) 

        if not os.path.exists(merged_file_path):
             stderr_pear += "\n[ERROR] File merged.assembled.fastq was not created by PEAR."

        return html.Div([
            html.H5("Downloads PEAR"),
            html.A("Download merged PEAR output", href=f"/uploads/pear_results/{base1_orig}_{base2_orig}/{merged_file_name}", target="_blank", download=merged_file_name), 
            dbc.Accordion([
                dbc.AccordionItem(html.Pre(read_file_preview(merged_file_path)), title="View Merged FASTQ"),
                dbc.AccordionItem([html.H6("PEAR stdout"), html.Pre(stdout_pear), html.H6("PEAR stderr"), html.Pre(stderr_pear)], title="Log PEAR"),
            ], start_collapsed=True, always_open=False, style={"marginTop": "1em"})
        ])
    return ""

@app.callback(
    Output('analysis-results', 'children'),
    Input('run-analysis-btn', 'n_clicks'),
    State('analysis-method', 'value'),
    State('vsearch-identity', 'value'),
    State('uploaded-files', 'data'),
    prevent_initial_call=True
)
def run_analysis(n_clicks, method, identity, uploaded_files):
    if not n_clicks or not method:
        return ""
    if method == "vsearch":
        fastas = uploaded_files.get('analysis', []) if uploaded_files else []
        fasta_file = None
        
        if fastas and isinstance(fastas, list) and len(fastas) > 0:
            fasta_file = fastas[0]
        else:
            r1_orig = uploaded_files.get('preprocess_r1')
            r2_orig = uploaded_files.get('preprocess_r2')
            
            fastq_file_pear = None
            if r1_orig and r2_orig:
                base1_orig = os.path.splitext(os.path.basename(r1_orig))[0]
                base2_orig = os.path.splitext(os.path.basename(r2_orig))[0]
                pear_output_subdir = f"{base1_orig}_{base2_orig}"
                fastq_file_pear = os.path.join(UPLOAD_DIRECTORY, "pear_results", pear_output_subdir, "merged.assembled.fastq")

            se_files_trimmed = uploaded_files.get('preprocess', []) 
            fastq_file_single_trimmed = None
            if se_files_trimmed and isinstance(se_files_trimmed, list) and len(se_files_trimmed) > 0:
                r1_se_orig = se_files_trimmed[0]
                base_se_orig = os.path.splitext(os.path.basename(r1_se_orig))[0]
                fastq_file_single_trimmed = os.path.join(UPLOAD_DIRECTORY, f"{base_se_orig}_trimmed.fastq")


            if fastq_file_pear and os.path.exists(fastq_file_pear) and os.path.getsize(fastq_file_pear) > 0:
                fastq_file = fastq_file_pear
                fasta_file = os.path.join(UPLOAD_DIRECTORY, "vsearch_otus", "input_from_pear.fasta")
                os.makedirs(os.path.join(UPLOAD_DIRECTORY, "vsearch_otus"), exist_ok=True)
            elif fastq_file_single_trimmed and os.path.exists(fastq_file_single_trimmed) and os.path.getsize(fastq_file_single_trimmed) > 0:
                fastq_file = fastq_file_single_trimmed
                fasta_file = os.path.join(UPLOAD_DIRECTORY, "vsearch_otus", "input_from_trimmed_se.fasta")
                os.makedirs(os.path.join(UPLOAD_DIRECTORY, "vsearch_otus"), exist_ok=True)
            else:
                return html.Div("No FASTQ files (PEAR merged or SE trimmed) found for OTU clustering. Run Trimmomatic and/or PEAR first, or upload a FASTA file directly.")

            if not fasta_file or not fastq_to_fasta(fastq_file, fasta_file):
                return html.Div(f"Error converting FASTQ ({os.path.basename(fastq_file)}) to FASTA.")

        if not fasta_file or not os.path.exists(fasta_file) or os.path.getsize(fasta_file) == 0:
            return html.Div("Input FASTA file not found or empty for OTU clustering.")

        output_dir_vsearch = os.path.join(UPLOAD_DIRECTORY, "vsearch_otus")
        derep_out, cluster_out, otus_fasta, uc_file = run_vsearch_otus(fasta_file, output_dir_vsearch, VSEARCH_PATH, identity)
        
        if not otus_fasta or not os.path.exists(otus_fasta):
            return html.Div(f"Error in VSEARCH clustering. Check the input and logs.\nDerep log:\n{derep_out}\nCluster log:\n{cluster_out}")

        rel_otus_fasta_path = os.path.relpath(otus_fasta, UPLOAD_DIRECTORY) 
        rel_uc_file_path = os.path.relpath(uc_file, UPLOAD_DIRECTORY)

        return html.Div([
            html.H5("Downloads VSEARCH"),
            html.A("Download OTUs FASTA", href=f"/uploads/{rel_otus_fasta_path.replace(os.sep, '/')}", download=os.path.basename(otus_fasta), target="_blank"),
            html.Br(),
            html.A("Download Cluster Assignments (.uc)", href=f"/uploads/{rel_uc_file_path.replace(os.sep, '/')}", download=os.path.basename(uc_file), target="_blank"),
            dbc.Accordion([
                dbc.AccordionItem([
                    html.H6("View FASTA OTUs"),
                    html.Pre(read_file_preview(otus_fasta))
                ], title="View FASTA OTUs"),
                dbc.AccordionItem([
                    html.H6("View Clusters (.uc)"),
                    html.Pre(read_file_preview(uc_file))
                ], title="View Clusters (.uc)"),
                dbc.AccordionItem([
                    html.H6("VSEARCH Dereplication Log"),
                    html.Pre(derep_out),
                    html.H6("VSEARCH Clustering Log"),
                    html.Pre(cluster_out),
                ], title="Log VSEARCH"),
            ], start_collapsed=True, always_open=False, style={"marginTop": "1em"})
        ])
    elif method == "dada2":
        return html.Div("DADA2 analysis not implemented yet.")
    elif method == "kraken2":
        return html.Div("Kraken2 analysis not implemented yet.")
    return ""

@app.callback(
    Output('boldigger-results', 'children'),
    Input('run-boldigger-btn', 'n_clicks'),
    State('bold-db-dropdown', 'value'),
    State('bold-mode-dropdown', 'value'),
    State('bold-thresholds', 'value'),
    State('uploaded-files', 'data'),
    prevent_initial_call=True
)
def run_boldigger_callback(n_clicks, db_nr, mode, thresholds, uploaded_files):
    if not n_clicks:
        return ""
    
    fasta_path_input_step = None
    db_uploads = uploaded_files.get('databases', []) if uploaded_files else []
    if db_uploads and isinstance(db_uploads, list) and len(db_uploads) > 0:
        for uploaded_db_file in db_uploads:
            if uploaded_db_file and (uploaded_db_file.lower().endswith('.fasta') or uploaded_db_file.lower().endswith('.fa')):
                fasta_path_input_step = uploaded_db_file
                break
    
    if not fasta_path_input_step:
        fasta_path_input_step = os.path.join(UPLOAD_DIRECTORY, "vsearch_otus", "otus.fasta")

    if not fasta_path_input_step or not os.path.exists(fasta_path_input_step) or os.path.getsize(fasta_path_input_step) == 0:
        return html.Div("No FASTA file found for BOLDigger. Run VSEARCH first or upload a FASTA file in the 'Reference Databases' tab.")

    boldigger_input_fasta_dir = os.path.join(UPLOAD_DIRECTORY, "boldigger_input")
    os.makedirs(boldigger_input_fasta_dir, exist_ok=True)
    
    sanitized_input_basename = "boldigger_input_sanitized" 
    fasta_path_unique_sanitized = os.path.join(boldigger_input_fasta_dir, f"{sanitized_input_basename}.fasta")
    temp_unique_path = os.path.join(boldigger_input_fasta_dir, "temp_unique_for_boldigger.fasta")

    try:
        uniquify_fasta(fasta_path_input_step, temp_unique_path)
        sanitize_ids(temp_unique_path, fasta_path_unique_sanitized)
        if os.path.exists(temp_unique_path):
            os.remove(temp_unique_path)
    except Exception as e:
        return html.Div(f"Error preparing FASTA file for BOLDigger (uniquify/sanitize): {e}")

    if not os.path.exists(fasta_path_unique_sanitized) or os.path.getsize(fasta_path_unique_sanitized) == 0:
        return html.Div("Error creating sanitized FASTA file for BOLDigger. The file is empty or was not created.")

    cmd = [
        sys.executable, "-m", "boldigger3", "identify",
        "--db", str(db_nr),
        "--mode", str(mode),
        fasta_path_unique_sanitized 
    ]

    if thresholds and str(thresholds).strip():
        try:
            threshold_values = [str(int(val)) for val in str(thresholds).strip().split()]
            if threshold_values:
                cmd += ["--thresholds"] + threshold_values
        except ValueError:
            return html.Div([
                html.P("Invalid thresholds. Use numbers separated by spaces, e.g. 99 97 90 85 80", style={"color": "red"}),
                dbc.Accordion([
                    dbc.AccordionItem([
                        html.H6("BOLDigger3 stdout"),
                        html.Pre("N/A - Input error"),
                        html.H6("BOLDigger3 stderr"),
                        html.Pre("User error in configuring thresholds."),
                    ], title="Log BOLDigger"),
                ], start_collapsed=True, always_open=False, style={"marginTop": "1em"})
            ])

    print("Running BOLDigger command:", " ".join(cmd))
    

    original_cwd = os.getcwd()
    os.makedirs(BOLDIGGER_RESULTS_DIR, exist_ok=True) 
    
    result = None
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=False, cwd=BOLDIGGER_RESULTS_DIR)
    except Exception as e:
        os.chdir(original_cwd) 
        return html.Div(f"Failed to execute BOLDigger: {str(e)}")
    finally:
        os.chdir(original_cwd) 

    log_accordion = dbc.Accordion([
        dbc.AccordionItem([
            html.H6("BOLDigger3 stdout"),
            html.Pre(result.stdout if result and result.stdout else "No stdout."),
            html.H6("BOLDigger3 stderr"),
            html.Pre(result.stderr if result and result.stderr else "No stderr."),
        ], title="Log BOLDigger"),
    ], start_collapsed=True, always_open=False, style={"marginTop": "1em"})

    if result and result.returncode != 0:
        if "No object named additional_data" in (result.stderr or "") or \
           "No matches found for your query sequences" in (result.stderr or ""):
            return html.Div([
                html.H5("BOLDigger3 Output"),
                html.P("No results were found for your sequences in BOLD. "
                       "This can happen if the sequences are too divergent, short, "
                       "or do not represent typical regions of the BOLD. "
                       "Try using less restrictive thresholds, another search mode, "
                       "or check that your sequences are suitable for identification in BOLD.", style={"color": "orange"}),
                log_accordion
            ])
        return html.Div([
            html.H5("BOLDigger3 Output"),
            html.P(f"BOLDigger3 error (return code {result.returncode}). Check the logs below.", style={"color": "red"}),
            log_accordion
        ])
    elif not result: 
        return html.Div([html.H5("BOLDigger3 Output"), html.P("BOLDigger execution failed to produce a result object.", style={"color": "red"}), log_accordion])


    output_file_basename = f"{sanitized_input_basename}_identification_result"
    
    output_xlsx_name = f"{output_file_basename}.xlsx"
    output_parquet_name = f"{output_file_basename}.parquet"
    
    output_xlsx_path = os.path.join(boldigger_input_fasta_dir, output_xlsx_name)
    output_parquet_path = os.path.join(boldigger_input_fasta_dir, output_parquet_name)

    xlsx_exists = os.path.exists(output_xlsx_path) and os.path.getsize(output_xlsx_path) > 0
    parquet_exists = os.path.exists(output_parquet_path) and os.path.getsize(output_parquet_path) > 0
    
    main_content = [html.H5("BOLDigger3 Output")]

    if not xlsx_exists and not parquet_exists:
         main_content.append(html.P(f"BOLDigger3 ran, but no main result file ({output_xlsx_name} or {output_parquet_name}) was found in the expected location: {boldigger_input_fasta_dir}. Check the logs.", style={"color": "orange"}))
    else:
        if xlsx_exists:
            relative_xlsx_download_path = os.path.relpath(output_xlsx_path, UPLOAD_DIRECTORY).replace(os.sep, '/')
            main_content.append(html.A(f"Download BOLDigger3 Results ({output_xlsx_name})", href=f"/uploads/{relative_xlsx_download_path}", download=output_xlsx_name, target="_blank"))
            main_content.append(html.Br())
            
            try:
                df = pd.read_excel(output_xlsx_path)
                table_preview = dash_table.DataTable(
                    data=df.to_dict('records'),
                    columns=[{'name': i, 'id': i} for i in df.columns],
                    page_size=10,
                    style_table={'overflowX': 'auto', 'marginTop': '1em'},
                    style_header={'backgroundColor': 'rgb(230, 230, 230)', 'fontWeight': 'bold'},
                    style_cell={
                        'height': 'auto',
                        'minWidth': '100px', 'width': '150px', 'maxWidth': '250px',
                        'whiteSpace': 'normal',
                        'textAlign': 'left',
                        'padding': '5px'
                    },
                    filter_action="native",
                    sort_action="native",
                    sort_mode="multi",
                )
                main_content.extend([html.Hr(), html.H6("Results Table Preview:"), table_preview])
            except Exception as e:
                main_content.append(html.P(f"Error reading or displaying Excel file: {str(e)}", style={"color": "orange"}))
        else:
            main_content.append(html.P(f"File {output_xlsx_name} not found or is empty.", style={"color": "red"}))
        
        if parquet_exists:
            relative_parquet_download_path = os.path.relpath(output_parquet_path, UPLOAD_DIRECTORY).replace(os.sep, '/')
            main_content.append(html.A(f"Download BOLDigger3 Results ({output_parquet_name})", href=f"/uploads/{relative_parquet_download_path}", download=output_parquet_name, target="_blank"))
            main_content.append(html.Br())
        elif not xlsx_exists: 
            main_content.append(html.P(f"File {output_parquet_name} not found or is empty.", style={"color": "red"}))

    return html.Div(main_content + [log_accordion])


@app.callback(
    Output('phylo-results', 'children'),
    Input('run-phylo-btn', 'n_clicks'),
    State('raxml-bootstrap-input', 'value'),
    State('uploaded-files', 'data'),
    prevent_initial_call=True
)
def run_phylo(n_clicks, n_bootstrap, uploaded_files):
    if not n_clicks:
        return ""
    if n_bootstrap is None or n_bootstrap < 10:
        n_bootstrap = 100 

    phylo_uploads = uploaded_files.get('phylo', []) if uploaded_files else []
    otus_fasta_input = None
    if phylo_uploads and isinstance(phylo_uploads, list) and len(phylo_uploads) > 0:
        for uploaded_phy_file in phylo_uploads:
            if uploaded_phy_file and (uploaded_phy_file.lower().endswith('.fasta') or uploaded_phy_file.lower().endswith('.fa')):
                otus_fasta_input = uploaded_phy_file
                break
    
    if not otus_fasta_input:
        otus_fasta_input = os.path.join(UPLOAD_DIRECTORY, "vsearch_otus", "otus.fasta")

    if not otus_fasta_input or not os.path.exists(otus_fasta_input) or os.path.getsize(otus_fasta_input) == 0:
        return html.Div("FASTA file of OTUs not found. Run VSEARCH first or upload a FASTA file in the “Phylogenetic Analysis” tab.")

    phylo_intermediate_dir = os.path.join(UPLOAD_DIRECTORY, "phylo_intermediate")
    os.makedirs(phylo_intermediate_dir, exist_ok=True)
    
    otus_fasta_unique = os.path.join(phylo_intermediate_dir, "otus_unique_for_phylo.fasta")
    uniquify_fasta(otus_fasta_input, otus_fasta_unique)

    aligned_fasta = os.path.join(phylo_intermediate_dir, "aligned.fasta")
    filtered_aligned_fasta = os.path.join(phylo_intermediate_dir, "aligned_filtered.fasta")
    aligned_safe_ids = os.path.join(phylo_intermediate_dir, "aligned_safe_ids.fasta")
    
    tree_file = os.path.join(phylo_intermediate_dir, "tree.nwk") 
    raxml_output_dir = os.path.join(UPLOAD_DIRECTORY, "raxmlng_results")

    mafft_ok, mafft_log = run_mafft(otus_fasta_unique, aligned_fasta)
    if not mafft_ok:
        return html.Div(f"MAFFT alignment failed: {mafft_log}")

    filter_empty_aligned_fasta(aligned_fasta, filtered_aligned_fasta, min_real=10)
    sanitize_ids(filtered_aligned_fasta, aligned_safe_ids)

    if not os.path.exists(aligned_safe_ids) or os.path.getsize(aligned_safe_ids) == 0:
        return html.Div("No valid sequences remained after alignment and filtering for phylogenetic analysis.")

    fasttree_ok, fasttree_log = run_fasttree(aligned_safe_ids, tree_file)
    raxml_ok, raxml_log = run_raxml(aligned_safe_ids, raxml_output_dir, n_bootstrap)
    
    raxml_best_tree_leaf_filename = "raxml.raxml.bestTree"
    raxml_best_tree_file = os.path.join(raxml_output_dir, raxml_best_tree_leaf_filename) 

    rel_aligned_safe_ids_path = os.path.relpath(aligned_safe_ids, UPLOAD_DIRECTORY)
    rel_tree_file_path = os.path.relpath(tree_file, UPLOAD_DIRECTORY)
    rel_raxml_best_tree_path = None
    if raxml_ok and os.path.exists(raxml_best_tree_file) and os.path.getsize(raxml_best_tree_file) > 0:
        rel_raxml_best_tree_path = os.path.relpath(raxml_best_tree_file, UPLOAD_DIRECTORY)

    fasttree_viz_content = html.P("FastTree failed or tree file is empty/invalid.", style={"color": "red"})
    if fasttree_ok and os.path.exists(tree_file) and os.path.getsize(tree_file) > 0:
        try:
            fasttree_obj = Phylo.read(tree_file, "newick")
            fasttree_elements = biopython_to_cytoscape_elements(fasttree_obj)
            
            root_node_id_ft = None
            if fasttree_elements and fasttree_elements[0] and 'data' in fasttree_elements[0] and 'id' in fasttree_elements[0]['data']:
                 root_node_id_ft = f"#{fasttree_elements[0]['data']['id']}"

            if fasttree_elements:
                fasttree_viz_content = cyto.Cytoscape(
                    id='cytoscape-fasttree',
                    elements=fasttree_elements,
                    layout={
                        'name': 'breadthfirst',
                        'directed': True,
                        'spacingFactor': 1.25, 
                        'grid': False,
                        'roots': root_node_id_ft
                    },
                    style={'width': '100%', 'height': '800px', 'border': '1px solid #ccc'}, 
                    stylesheet=[
                        {'selector': 'node', 'style': {'label': 'data(label)', 'font-size': '9px', 'text-valign': 'center', 'text-halign': 'center', 'background-color': '#ADD8E6', 'width': 'mapData(label.length, 0, 20, 15, 80)', 'height': 'mapData(label.length, 0, 20, 15, 80)', 'text-wrap': 'wrap', 'text-max-width': '70px'}},
                        {'selector': 'edge', 'style': {'width': 1, 'line-color': '#9dbaea', 'curve-style': 'bezier'}}, 
                        {'selector': '[label = ""]', 'style':{'width': '5px', 'height': '5px', 'background-color': 'grey', 'font-size': '7px'}}, 
                        {'selector': 'node[!label]', 'style':{'width': '5px', 'height': '5px', 'background-color': 'grey', 'font-size': '7px'}}, 
                        {'selector': 'node[label^="."]', 'style':{'font-size': '7px', 'width': '5px', 'height': '5px'}}, 
                    ]
                )
            else:
                fasttree_viz_content = html.P("Could not generate FastTree visualization elements.", style={"color": "orange"})
        except Exception as e:
            fasttree_viz_content = html.P(f"Error parsing FastTree Newick for visualization: {e}", style={"color": "red"})
    elif not fasttree_ok:
         fasttree_viz_content = html.P(f"FastTree failed: {fasttree_log}", style={"color": "red"})

    raxml_viz_content = html.P("RAxML-NG failed or tree file is empty/invalid.", style={"color": "red"})
    if raxml_ok and rel_raxml_best_tree_path and os.path.exists(raxml_best_tree_file) and os.path.getsize(raxml_best_tree_file) > 0:
        try:
            raxml_tree_obj = Phylo.read(raxml_best_tree_file, "newick")
            raxml_elements = biopython_to_cytoscape_elements(raxml_tree_obj)

            root_node_id_raxml = None
            if raxml_elements and raxml_elements[0] and 'data' in raxml_elements[0] and 'id' in raxml_elements[0]['data']:
                 root_node_id_raxml = f"#{raxml_elements[0]['data']['id']}"

            if raxml_elements:
                raxml_viz_content = cyto.Cytoscape(
                    id='cytoscape-raxml',
                    elements=raxml_elements,
                    layout={
                        'name': 'breadthfirst',
                        'directed': True,
                        'spacingFactor': 1.25, 
                        'grid': False,
                        'roots': root_node_id_raxml
                    },
                    style={'width': '100%', 'height': '800px', 'border': '1px solid #ccc'}, 
                    stylesheet=[ 
                        {'selector': 'node', 'style': {'label': 'data(label)', 'font-size': '9px', 'text-valign': 'center', 'text-halign': 'center', 'background-color': '#90EE90', 'width': 'mapData(label.length, 0, 20, 15, 80)', 'height': 'mapData(label.length, 0, 20, 15, 80)', 'text-wrap': 'wrap', 'text-max-width': '70px'}},
                        {'selector': 'edge', 'style': {'width': 1, 'line-color': '#9dbaea', 'curve-style': 'bezier'}},
                        {'selector': '[label = ""]', 'style':{'width': '5px', 'height': '5px', 'background-color': 'grey', 'font-size': '7px'}},
                        {'selector': 'node[!label]', 'style':{'width': '5px', 'height': '5px', 'background-color': 'grey', 'font-size': '7px'}},
                        {'selector': 'node[label^="."]', 'style':{'font-size': '7px', 'width': '5px', 'height': '5px'}},
                    ]
                )
            else:
                raxml_viz_content = html.P("Could not generate RAxML-NG visualization elements.", style={"color": "orange"})
        except Exception as e:
            raxml_viz_content = html.P(f"Error parsing RAxML-NG Newick for visualization: {e}", style={"color": "red"})
    elif not raxml_ok:
        raxml_viz_content = html.P(f"RAxML-NG failed. Check log.", style={"color": "red"})

    return html.Div([
        html.H5("Downloads Phylogenetics"),
        html.A("Download Aligned FASTA (Sanitized)", href=f"/uploads/{rel_aligned_safe_ids_path.replace(os.sep, '/')}", download=os.path.basename(aligned_safe_ids), target="_blank"),
        html.Br(),
        html.A("Download Newick Tree (FastTree)", href=f"/uploads/{rel_tree_file_path.replace(os.sep, '/')}", download=os.path.basename(tree_file), target="_blank") if fasttree_ok and os.path.exists(tree_file) else html.P("FastTree Newick file not available.", style={"color":"orange"}),
        html.Br(),
        (html.A("Download RAxML Best Tree", href=f"/uploads/{rel_raxml_best_tree_path.replace(os.sep, '/')}", download=os.path.basename(raxml_best_tree_file), target="_blank")
         if rel_raxml_best_tree_path else html.P("RAxML best tree not found or RAxML failed.", style={"color": "red" if not raxml_ok else "orange"})),
        html.Br(),
        html.Hr(),
        dbc.Accordion([
            dbc.AccordionItem(html.Pre(read_file_preview(aligned_safe_ids)), title="View FASTA Alignment (Sanitized)"),
            dbc.AccordionItem(fasttree_viz_content, title="View Phylogenetic Tree (FastTree)"),
            dbc.AccordionItem(raxml_viz_content, title="View Phylogenetic Tree (RAxML-NG)"),
            dbc.AccordionItem([html.H6("MAFFT Alignment Log"), html.Pre(mafft_log if mafft_log else "No MAFFT log output.")], title="Log MAFFT"),
            dbc.AccordionItem([html.H6("FastTree Log"), html.Pre(fasttree_log if fasttree_log else "No FastTree log output.")], title="Log FastTree"),
            dbc.AccordionItem([html.H6("RAxML-NG Output"), html.Pre(raxml_log if raxml_log else "No RAxML-NG log output.")], title="Log RAxML-NG"),
        ], start_collapsed=True, always_open=False, style={"marginTop": "2em"})
    ])

@app.server.route('/uploads/<path:path>')
def serve_file(path):
    return send_from_directory(UPLOAD_DIRECTORY, path)

@app.callback(
    Output('workflow-progress', 'data', allow_duplicate=True),
    Input('run-qc-btn', 'n_clicks'),
    State('workflow-progress', 'data'),
    prevent_initial_call=True
)
def mark_qc_done(n, progress):
    if n:
        progress['qc_done'] = True
    return progress

@app.callback(
    Output('workflow-progress', 'data', allow_duplicate=True),
    [Input('run-trim-btn', 'n_clicks'), Input('run-pear-btn', 'n_clicks')],
    State('workflow-progress', 'data'),
    prevent_initial_call=True
)
def mark_preprocess_done(n_trim, n_pear, progress):
    if n_trim or n_pear:
        progress['preprocess_done'] = True
    return progress

@app.callback(
    Output('workflow-progress', 'data', allow_duplicate=True),
    Input('run-analysis-btn', 'n_clicks'),
    State('workflow-progress', 'data'),
    prevent_initial_call=True
)
def mark_analysis_done(n, progress):
    if n:
        progress['analysis_done'] = True
    return progress

@app.callback(
    Output('workflow-progress', 'data', allow_duplicate=True),
    Input('db-selection', 'value'),
    State('workflow-progress', 'data'),
    prevent_initial_call=True
)
def mark_db_done(selected, progress):
    if selected:
        progress['db_done'] = True
    return progress

@app.callback(
    Output('workflow-progress', 'data', allow_duplicate=True),
    Input('run-phylo-btn', 'n_clicks'),
    State('workflow-progress', 'data'),
    prevent_initial_call=True
)
def mark_phylo_done(n, progress):
    if n:
        progress['phylo_done'] = True
    return progress

@app.callback(
    Output('workflow-progress', 'data', allow_duplicate=True),
    Input('final-results', 'children'),
    State('workflow-progress', 'data'),
    prevent_initial_call=True
)
def mark_results_done(children, progress):
    if children:
        progress['results_done'] = True
    return progress

@app.callback(
    Output('db-info', 'children'),
    Input('db-selection', 'value')
)
def show_db_info(selected):
    if not selected:
        return html.Div("No databases selected.")

    db_descriptions = {
        "ncbi": "NCBI GenBank: Largest nucleotide database, broad taxonomic coverage.",
        "silva": "SILVA: Curated rRNA database (16S, 18S), high quality, regular updates.",
        "pr2": "PR2: Protist Ribosomal Reference database, curated 18S rRNA sequences.",
        "bold": "BOLD: Barcode of Life Data System, COI gene focus, accessed via BOLDigger.",
    }

    return html.Ul([html.Li(db_descriptions[db]) for db in selected if db in db_descriptions])

if __name__ == '__main__':
    app.run(debug=True)