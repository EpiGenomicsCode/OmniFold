import argparse
import json
import re
import shutil
import subprocess
import zipfile
from pathlib import Path
from typing import Dict, List, Any

import numpy as np
import plotly.graph_objects as go
from jinja2 import Environment, FileSystemLoader

# --- Data Parsing Functions ---

def _get_plddt_and_chains_from_cif(cif_file: Path):
    """
    Extracts pLDDT for C-alpha atoms of standard amino acid residues and detailed chain and residue information from a CIF file.
    Returns a tuple of (plddt_list, chain_info_dict, res_indices_list).
    """
    plddts = []
    chain_info = {}
    res_indices = []

    if not cif_file.is_file():
        print(f"Warning: CIF file not found for pLDDT extraction: {cif_file}")
        return [], {}, []
    try:
        # Use MMCIFParser to get the structure
        from Bio.PDB import MMCIFParser
        from Bio.PDB.Polypeptide import is_aa
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure(cif_file.stem, str(cif_file))
        model = structure[0] # Assuming only one model for simplicity
        
        for chain in model:
            res_list = [res for res in chain if is_aa(res.get_resname(), standard=True)]
            num_residues = len(res_list)
            
            if num_residues > 0:
                chain_info[chain.id] = {'count': num_residues, 'type': 'protein'}
                total_residues = sum(chain_info[c_id]['count'] for c_id in chain_info.keys())
                res_indices.extend([res.id[1] for res in res_list])
                
                # Extract pLDDT for C-alpha atoms
                for residue in res_list:
                    if "CA" in residue:
                        plddts.append(residue["CA"].get_bfactor())

    except Exception as e:
        print(f"Warning: Could not parse {cif_file} for pLDDT or chain info: {e}")
        return [], {}, []

    return plddts, chain_info, res_indices

def parse_all_af3_outputs(directory: Path) -> List[Dict[str, Any]]:
    """Parses all AlphaFold 3 model outputs in a directory and sorts them by the master ranking scores CSV."""
    all_models = []
    
    # --- Find and parse the master ranking CSV ---
    ranking_csv = next(directory.glob("**/*_ranking_scores.csv"), None)
    ordered_model_stems = []
    if ranking_csv:
        import pandas as pd
        try:
            df = pd.read_csv(ranking_csv)
            # Reconstruct the model stem to match the file names
            csv_stem = ranking_csv.stem.replace('_ranking_scores', '')
            # Force seed/sample to be formatted as integers in the string
            df['model_stem'] = df.apply(lambda row: f"{csv_stem}_seed-{int(row['seed'])}_sample-{int(row['sample'])}", axis=1)
            
            # Sort by ranking_score, highest first
            df_sorted = df.sort_values(by='ranking_score', ascending=False)
            ordered_model_stems = df_sorted['model_stem'].tolist()
        except Exception as e:
            print(f"Warning: Could not parse AlphaFold 3 ranking CSV '{ranking_csv}': {e}")
    
    # --- Gather all model data ---
    model_data_map = {}
    cif_files = list(directory.glob("**/*_seed-*-*_model.cif"))
    if not cif_files: return []

    for cif_file in cif_files:
        conf_stem = cif_file.stem.replace('_model', '')
        summary_path = cif_file.parent / f"{conf_stem}_summary_confidences.json"
        conf_path = cif_file.parent / f"{conf_stem}_confidences.json"
        if not summary_path.exists() or not conf_path.exists(): continue
        
        try:
            summary_data = json.loads(summary_path.read_text())
            conf_data = json.loads(conf_path.read_text())
            plddt, chain_info, res_indices = _get_plddt_and_chains_from_cif(cif_file)
            
            # Use the conf_stem as the key for matching with the CSV
            model_data_map[conf_stem] = {
                'name': f"AlphaFold 3 (Sample {conf_stem.split('sample-')[-1]})",
                'model_class': 'alphafold-3',
                'ranking_score': summary_data.get('ranking_score'),
                'ptm': summary_data.get('ptm'), 'iptm': summary_data.get('iptm'),
                'plddt': plddt, 'avg_plddt': np.mean(plddt) if plddt else 0,
                'cif_path': cif_file, 'pae_path': conf_path,
                'chain_info': chain_info, 'res_indices': res_indices,
                'pae_matrix': np.array(conf_data['pae']),
            }
        except (IOError, json.JSONDecodeError): continue
    
    # --- Sort models based on the CSV order ---
    if ordered_model_stems:
        # Reorder the gathered models according to the CSV ranking
        for stem in ordered_model_stems:
            if stem in model_data_map:
                all_models.append(model_data_map[stem])

        # Add any models found that were not in the CSV (less ideal, but robust)
        for stem, model in model_data_map.items():
             if model not in all_models: # Check if model object is already added
                 all_models.append(model)
    else:
        # Fallback to original sorting if CSV is missing or fails to parse
        all_models = sorted(model_data_map.values(), key=lambda x: x.get('ranking_score', 0), reverse=True)

    return all_models

def parse_all_chai_outputs(directory: Path) -> List[Dict[str, Any]]:
    """Parses all Chai-1 model outputs and sorts by ranking score."""
    all_models = []
    score_files = sorted(list(directory.glob("scores.model_idx_*.npz")))
    if not score_files: return []

    for score_file in score_files:
        model_index = score_file.stem.split('_')[-1]
        cif_file = directory / f"pred.model_idx_{model_index}.cif"
        pae_file = directory / f"pae.model_idx_{model_index}.npz"
        if not cif_file.exists() or not pae_file.exists(): continue

        try:
            scores_data = np.load(score_file)
            plddt, chain_info, res_indices = _get_plddt_and_chains_from_cif(cif_file)
            pae_matrix = np.load(pae_file)['pae']
            if pae_matrix.ndim == 3: pae_matrix = pae_matrix[0]

            all_models.append({
                'name': f"Chai-1 (model {model_index})", 'model_class': 'chai-1',
                'ranking_score': scores_data.get('aggregate_score'),
                'ptm': scores_data.get('ptm'), 'iptm': scores_data.get('iptm'),
                'plddt': plddt, 'avg_plddt': np.mean(plddt) if plddt else 0,
                'cif_path': cif_file, 'pae_path': pae_file,
                'chain_info': chain_info, 'res_indices': res_indices,
                'pae_matrix': pae_matrix
            })
        except (IOError, KeyError): continue
            
    return sorted(all_models, key=lambda x: x.get('ranking_score', 0), reverse=True)

def parse_all_boltz_outputs(directory: Path) -> List[Dict[str, Any]]:
    """Finds and parses the Boltz-1 output."""
    # Boltz-1 typically produces one model, but we'll return a list for consistency.
    job_dir = next(directory.glob("boltz_results_*"), None)
    if not job_dir: return []
    
    pred_dir = next((job_dir / "predictions").iterdir(), None)
    if not pred_dir: return []
    
    cif_file = next(pred_dir.glob("*_model_0.cif"), None)
    if not cif_file: return []

    conf_file = pred_dir / f"confidence_{cif_file.stem}.json"
    pae_file = pred_dir / f"pae_{cif_file.stem}.npz"
    if not conf_file.exists(): return []

    try:
        conf_data = json.loads(conf_file.read_text())
        plddt, chain_info, res_indices = _get_plddt_and_chains_from_cif(cif_file)
        
        model = {
            'name': "Boltz 2", 'model_class': 'boltz-2',
            'ranking_score': conf_data.get('confidence_score'),
            'ptm': conf_data.get('ptm'), 'iptm': conf_data.get('iptm'),
            'plddt': plddt, 'avg_plddt': np.mean(plddt) if plddt else 0,
            'cif_path': cif_file, 'chain_info': chain_info, 'res_indices': res_indices
        }
        if pae_file.exists():
            model['pae_path'] = pae_file
            model['pae_matrix'] = np.load(pae_file)['pae']

        return [model]
    except (IOError, KeyError): return []


# --- PAE Viewer Generation ---
def generate_pae_viewer(model: Dict[str, Any], pae_output_dir: Path) -> str:
    """Generates a standalone PAE viewer for a given model and returns the relative path."""
    viewer_script_dir = Path(__file__).parent / "pae-viewer"
    export_script_path = viewer_script_dir / "resources" / "scripts" / "export_report.py"
    
    if not export_script_path.exists():
        print(f"Warning: PAE viewer export script not found at {export_script_path}")
        return "#"

    model_class = model.get('model_class')
    output_name = f"pae_report_{model['name'].replace(' ', '_').replace('(', '').replace(')', '')}"
    
    cmd = ["python3", str(export_script_path), "--output", output_name]
    
    # Extract labels from the model's chain info, which is derived from the output CIF.
    if 'chain_info' in model and model['chain_info']:
        labels = ";".join(model['chain_info'].keys())
        print(f"[DEBUG] {model['name']} chain labels resolved to: {labels}")
    else:
        print(f"Warning: No chain_info found for model {model['name']}. Cannot generate PAE viewer.")
        return "#"

    if model_class == 'alphafold-3':
        if 'cif_path' not in model or 'pae_path' not in model: return "#"
        cmd.extend(["--structure", str(model['cif_path']), "--scores", str(model['pae_path']), "--labels", labels])
    
    elif model_class == 'chai-1':
        model_index = re.search(r'model (\d+)', model['name']).group(1)
        chai_dir = model['cif_path'].parent
        # Do NOT pass explicit labels – let export_report.py derive order from CIF.
        cmd.extend(["--chai-input-dir", str(chai_dir), "--model-index", model_index])
        
    elif model_class == 'boltz-2': # Corrected from boltz-1
        # The export script expects the predictions directory (one level above the model subfolder)
        pred_dir = model['cif_path'].parent.parent
        boltz_model_name = model['cif_path'].stem.replace('_model_0', '')
        # Pass labels to ensure correct chain ordering
        cmd.extend(["--boltz-input-dir", str(pred_dir), "--boltz-model-name", boltz_model_name, "--labels", labels])
        
    else:
        return "#"

    try:
        print(f"[DEBUG] Running PAE viewer export command: {' '.join(cmd)}")
        # Run from the viewer script's directory to resolve its internal relative paths
        result = subprocess.run(cmd, check=True, capture_output=True, text=True, cwd=viewer_script_dir)
        print(result.stdout)
        if result.stderr:
            print('[STDERR]', result.stderr)
        
        # Find the output file path from the script's stdout
        match = re.search(r"Successfully created standalone report: (.*\.html)", result.stdout)
        if match:
            # Move the report to the central pae_viewers directory
            generated_file_path = Path(match.group(1).strip())
            destination_path = pae_output_dir / generated_file_path.name
            shutil.move(str(generated_file_path), str(destination_path))
            # Return a relative path for the link in the main report
            return f"pae_viewers/{destination_path.name}"
        else:
            print(f"Warning: Could not find PAE report path in stdout for {model['name']}:\n{result.stdout}")
            return "#"
            
    except subprocess.CalledProcessError as e:
        print(f"Warning: PAE viewer generation failed for {model['name']}:\n{e.stderr}")
        return "#"


# --- ipSAE Execution and Parsing ---

def run_and_parse_ipsae(model: Dict[str, Any]) -> Dict[str, Dict]:
    """Runs ipsae.py on a model and parses the resulting text file."""
    if 'pae_path' not in model or 'cif_path' not in model: return {}

    script_path = Path(__file__).parent / "ipsae.py"
    output_dir = model['cif_path'].parent
    
    # This is the filename that ipsae.py will create in its own directory
    cif_stem = model['cif_path'].stem
    expected_filename_in_script_dir = f"{cif_stem}_10_10.txt"
    expected_path_in_script_dir = script_path.parent / expected_filename_in_script_dir

    # This is the destination path in the model's directory
    final_output_path = output_dir / f"ipsae_{expected_filename_in_script_dir}"

    cmd = [
        "python", str(script_path),
        str(model['pae_path']), str(model['cif_path']),
        "10", "10", 
        "ignored_output_stem" # This argument is ignored by the script but required
    ]
    
    # Clean up old file if it exists
    if expected_path_in_script_dir.exists():
        expected_path_in_script_dir.unlink()

    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        # Move the output file from the script dir to the model's dir
        if expected_path_in_script_dir.exists():
            shutil.move(str(expected_path_in_script_dir), str(final_output_path))
        else:
             print(f"Warning: Expected ipSAE output not found at {expected_path_in_script_dir}")
             return {}
    except subprocess.CalledProcessError as e:
        print(f"Warning: ipsae.py failed for {model['name']}:\n{e.stderr}")
        return {}
    
    if not final_output_path.exists(): 
        print(f"Warning: ipSAE output file does not exist after moving: {final_output_path}")
        return {}

    interface_scores = {}
    with open(final_output_path, 'r') as f:
        for line in f:
            if "max" not in line: continue # We only want the 'max' summary line per pair
            parts = line.split()
            if len(parts) < 10: continue
            
            chain1, chain2, _, _, _, ipsae, _, _, _, _, pdockq, *_ = parts
            pair_key = f"{chain1}-{chain2}"
            interface_scores[pair_key] = {
                'ipsae': float(ipsae), 'pdockq': float(pdockq)
            }
            
    return interface_scores

# --- Plotting and Main Logic ---

def create_plddt_plot(all_models_data, best_model_names):
    """Creates an interactive pLDDT plot using Plotly, with togglable traces for all models."""
    fig = go.Figure()
    colors = {"AlphaFold 3": "#0054a6", "Boltz 2": "#f58220", "Chai-1": "#2ca02c"}

    # Use the chain info from the first "best" model for consistency in annotations
    first_best_model = next((m for m in all_models_data if m['name'] in best_model_names), None)

    for model in all_models_data:
        plddt_array = model.get('plddt', [])
        model_name = model.get('name', 'Unknown')
        
        line_color = 'grey' # Default color for non-best models in legend
        if "Boltz 2" in model_name: line_color = colors["Boltz 2"]
        elif "AlphaFold 3" in model_name: line_color = colors["AlphaFold 3"]
        elif "Chai-1" in model_name: line_color = colors["Chai-1"]
        
        # Set visibility: True for best models, 'legendonly' for others
        is_best = model_name in best_model_names
        
        fig.add_trace(go.Scatter(
            x=list(range(len(plddt_array))), y=plddt_array, name=model_name,
            mode='lines', line=dict(color=line_color, width=2.5),
            visible=True if is_best else 'legendonly',
            hovertemplate='Residue: %{x}<br>pLDDT: %{y:.2f}<extra></extra>'
        ))
    
    # Add confidence background rectangles for visual clarity
    fig.add_hrect(y0=90, y1=100, line_width=0, fillcolor="#c8e6c9", opacity=0.4, layer="below", annotation_text="Very High", annotation_position="right")
    fig.add_hrect(y0=70, y1=90, line_width=0, fillcolor="#bbdefb", opacity=0.4, layer="below", annotation_text="Confident", annotation_position="right")
    fig.add_hrect(y0=50, y1=70, line_width=0, fillcolor="#ffecb3", opacity=0.4, layer="below", annotation_text="Low", annotation_position="right")
    fig.add_hrect(y0=0, y1=50, line_width=0, fillcolor="#ffcdd2", opacity=0.4, layer="below", annotation_text="Very Low", annotation_position="right")

    # Add chain separator lines and labels based on the first best model
    if first_best_model and 'chain_info' in first_best_model and first_best_model['chain_info']:
        chain_info = first_best_model['chain_info']
        chain_names = list(chain_info.keys())
        
        # Calculate chain boundaries
        boundaries = np.cumsum([info['count'] for info in chain_info.values()]).tolist()
        boundaries.pop() # Remove the last boundary which is the total length

        # Add vertical lines
        for pos in boundaries:
            fig.add_vline(x=pos - 0.5, line=dict(color="grey", dash="dash", width=1), opacity=0.7)

        # Add chain labels as annotations in the middle of each chain's region
        last_pos = 0
        for chain_name in chain_names:
            if chain_name not in chain_info: continue
            chain_len = chain_info[chain_name]['count']
            mid_point = last_pos + (chain_len / 2)
            
            fig.add_annotation(
                x=mid_point,
                y=101, # Position the label just above the plot area
                text=f"<b>Chain {chain_name}</b>",
                showarrow=False,
                font=dict(size=11, color="#444"),
                xanchor='center'
            )
            last_pos += chain_len

    fig.update_layout(
        xaxis_title='<b>Residue Index</b>', yaxis_title='<b>pLDDT</b>',
        plot_bgcolor='white',
        paper_bgcolor='white',
        font=dict(family="-apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif", color="#212529"),
        xaxis=dict(gridcolor='#dee2e6', zeroline=False),
        yaxis=dict(
            gridcolor='#dee2e6',
            range=[40, 102],  # Zoom in on the most relevant pLDDT range
            zeroline=False
        ),
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02, # Position above the plot area
            xanchor="center",
            x=0.5,
            bgcolor='rgba(255, 255, 255, 0)',
            bordercolor='#ced4da',
            borderwidth=0
        ),
        margin=dict(l=60, r=20, t=80, b=80) # Adjust top margin for legend
    )
    return fig.to_html(full_html=False, include_plotlyjs=False)

# --- Main Logic ---

def run_report_generation(base_output_dir: Path):
    """
    Generates the full HTML report and accompanying files.
    This is the main entry point when called as a library function.
    """
    report_path = base_output_dir / "final_report.html"
    
    # Create a dedicated directory for PAE viewers
    pae_viewer_dir = base_output_dir / "pae_viewers"
    pae_viewer_dir.mkdir(exist_ok=True)

    # --- Data Gathering ---
    af_dir = base_output_dir / "alphafold3"
    boltz_dir = base_output_dir / "boltz"
    chai_dir = base_output_dir / "chai1"

    all_predictions = {}
    if af_dir.exists():
        all_predictions["AlphaFold 3"] = parse_all_af3_outputs(af_dir)
    if chai_dir.exists():
        all_predictions["Chai-1"] = parse_all_chai_outputs(chai_dir)
    if boltz_dir.exists():
        all_predictions["Boltz 2"] = parse_all_boltz_outputs(boltz_dir)


    # --- Generate PAE Viewers ---
    print("\nGenerating PAE viewers...")
    all_models_flat = []
    for method in all_predictions:
        for model in all_predictions[method]:
            model['pae_viewer_path'] = generate_pae_viewer(model, pae_viewer_dir)
            all_models_flat.append(model)


    best_models = [preds[0] for preds in all_predictions.values() if preds]
    best_models_names = [m['name'] for m in best_models]
    
    summary_metrics = ['avg_plddt', 'ptm', 'iptm', 'ranking_score']
    best_overall = {}
    for metric in summary_metrics:
        valid_models = [m for m in best_models if m.get(metric) is not None]
        if valid_models:
            best_overall[metric] = max(valid_models, key=lambda x: x[metric])['name']

    interface_data_by_pair = {}
    all_pairs = set()
    for model in best_models:
        model_interface_scores = run_and_parse_ipsae(model)
        for pair, scores in model_interface_scores.items():
            all_pairs.add(pair)
            if pair not in interface_data_by_pair:
                interface_data_by_pair[pair] = {}
            interface_data_by_pair[pair][model['name']] = scores

    for pair in all_pairs:
        for model_name in best_models_names:
            if model_name not in interface_data_by_pair.get(pair, {}):
                if pair not in interface_data_by_pair: interface_data_by_pair[pair] = {}
                interface_data_by_pair[pair][model_name] = None
    
    plddt_plot_html = create_plddt_plot(all_models_flat, best_models_names)

    try:
        plotly_js_path = Path(__file__).parent / "plotly.min.js"
        plotly_js = plotly_js_path.read_text(encoding="utf-8")
    except IOError:
        print("Warning: plotly.min.js not found. Plot will not be interactive.")
        plotly_js = None

    env = Environment(loader=FileSystemLoader(Path(__file__).parent))
    template = env.get_template("template.html")
    
    html_content = template.render(
        all_predictions=all_predictions,
        best_overall=best_overall,
        best_models_names=best_models_names,
        interface_data_by_pair=interface_data_by_pair,
        plddt_plot=plddt_plot_html,
        plotly_js=plotly_js
    )
    
    report_path.write_text(html_content, encoding="utf-8")
    print(f"\n✅ Report successfully generated at: {report_path}")

    zip_path = base_output_dir / "OmniFold_Report.zip"
    print(f"Creating distributable ZIP file at: {zip_path}")
    with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
        zipf.write(report_path, report_path.name)
        if pae_viewer_dir.exists():
            for pae_file in pae_viewer_dir.glob("*.html"):
                zipf.write(pae_file, f"{pae_viewer_dir.name}/{pae_file.name}")
    print("✅ ZIP file created successfully.")


def main():
    """CLI entry point for generating the HTML report independently."""
    parser = argparse.ArgumentParser(description="Generate a self-contained HTML report from existing model outputs.")
    parser.add_argument("--output_dir", type=Path, required=True, help="The main output directory containing the model subdirectories (alphafold3, boltz, chai1).")
    args = parser.parse_args()
    run_report_generation(args.output_dir)


if __name__ == "__main__":
    main()
