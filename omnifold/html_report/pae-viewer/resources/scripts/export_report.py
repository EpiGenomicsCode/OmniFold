#!/usr/bin/env python

"""
Create a standalone HTML report for the PAE viewer with all dependencies inlined.
This script bundles all CSS, JavaScript, and HTML templates into a single
HTML file, and embeds user-provided data to create a portable, offline-first
viewer.
"""

import argparse
import base64
import json
import mimetypes
import os
import re
import shutil
import subprocess
import webbrowser
from datetime import datetime
from pathlib import Path
from typing import Optional, Any, List
import numpy as np
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning
import warnings
from Bio.PDB import MMCIFParser, PDBIO
import io

# Suppress Biopython's PDBConstructionWarning
warnings.simplefilter('ignore', PDBConstructionWarning)


def resolved_path(value: str) -> Path:
    """Validate and resolve a file path."""
    path = Path(value).resolve()
    if not path.exists():
        raise argparse.ArgumentTypeError(f"File not found: {path}")
    return path


def inline_css(project_root: Path, css_paths: List[Path]) -> str:
    """Read CSS files, embed linked assets as data URIs, and return as a single string."""
    style_content = []
    # This regex now includes a negative lookahead `(?!data:)` to skip existing data URIs.
    url_pattern = re.compile(r'url\((.*?)\)')

    for css_path in css_paths:
        with open(css_path, "r", encoding="utf-8") as f:
            content = f.read()
            # Find all url() references that are not already data URIs
            content = url_pattern.sub(
                lambda match: embed_asset(match, css_path.parent),
                content
            )
            style_content.append(content)
    return "\n".join(style_content)


def embed_asset(match: re.Match, base_path: Path) -> str:
    """Convert a matched file path in CSS to a data URI."""
    url_raw = match.group(1)
    # Strip quotes before checking, per user feedback
    url_path = url_raw.strip("'\"")

    if url_path.startswith("data:"):
        return f'url({url_raw})'

    asset_path = (base_path / url_path).resolve()
    if not asset_path.is_file():
        # This can be noisy with bootstrap's embedded SVGs, return original
        return f'url({url_raw})'

    mime_type, _ = mimetypes.guess_type(asset_path)
    if mime_type is None:
        return f'url({url_raw})'

    with open(asset_path, "rb") as f:
        asset_bytes = f.read()

    encoded_asset = base64.b64encode(asset_bytes).decode("utf-8")
    return f'url("data:{mime_type};base64,{encoded_asset}")'


def inline_html_templates(templates_dir: Path) -> str:
    """Embed all .tpl files into script tags."""
    template_elements = []
    for tpl_file in sorted(templates_dir.glob("*.tpl")):
        with open(tpl_file, "r", encoding="utf-8") as f:
            content = f.read()
            template_id = tpl_file.name
            template_elements.append(
                f'<script type="text/template" id="{template_id}">\n{content}\n</script>'
            )
    return "\n".join(template_elements)


def inline_js_libs(project_root: Path, lib_paths: List[Path]) -> str:
    """Concatenate JS libraries, escaping literal </script> tags."""
    script_content = []
    for lib_path in lib_paths:
        with open(lib_path, "r", encoding="utf-8") as f:
            content = f.read()
            # Strip sourceMappingURL comments, per user feedback
            content = re.sub(r"//#.*sourceMappingURL=.*", "", content)
            content = re.sub(r"/\*#.*?sourceMappingURL=.*?\*/", "", content, flags=re.S)
            # Escape literal </script> tags to prevent premature HTML parsing termination
            content = content.replace("</script>", "<\\/script>")
            script_content.append(content)
    return "\n".join(script_content)


def bundle_js_app(project_root: Path, entry_point: Path) -> str:
    """Bundle the application's ES6 modules using esbuild."""
    output_file = project_root / "dist" / "app.bundle.js"
    output_file.parent.mkdir(exist_ok=True)

    # Per user feedback, create a more robust path to the esbuild executable
    esbuild_path = project_root / "node_modules" / "esbuild" / "bin" / "esbuild"

    # Check for local esbuild executable
    if not esbuild_path.is_file():
        # Fallback for different OS structures or if esbuild is in PATH
        esbuild_path_alt = project_root / "node_modules" / ".bin" / "esbuild"
        if esbuild_path_alt.is_file():
            esbuild_path = esbuild_path_alt
        else:
            print("Error: `esbuild` not found in `node_modules`.")
            print("Please run `npm install` to install dependencies.")
            exit(1)

    command = [
        str(esbuild_path),
        str(entry_point),
        "--bundle",
        "--format=iife",
        "--minify",
        f"--outfile={output_file}",
    ]
    
    try:
        subprocess.run(
            command,
            check=True,
            capture_output=True,
            text=True,
            cwd=project_root
        )
    except subprocess.CalledProcessError as e:
        print("Error during JavaScript bundling with esbuild:")
        print(e.stderr)
        exit(1)

    with open(output_file, "r", encoding="utf-8") as f:
        return f.read()


def convert_cif_to_pdb_string(cif_path: Path) -> str:
    """
    Reads a CIF file and converts it to a PDB formatted string using Biopython.
    """
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("structure", str(cif_path))
    
    # Use an in-memory string buffer to "write" the PDB file
    with io.StringIO() as buffer:
        pdb_io = PDBIO()
        pdb_io.set_structure(structure)
        pdb_io.save(buffer)
        return buffer.getvalue()


def get_chains_from_cif(cif_path: Path) -> list[str]:
    """
    Return the chain IDs **in the order they appear** in the first model of the mmCIF.

    We do *not* sort the chains alphabetically – the viewer expects the label list
    to correspond position-wise to the order in which chains are encountered by
    NGL when it parses the structure.  Keeping the original order avoids label →
    sequence mismatches.
    """
    parser = MMCIFParser(QUIET=True)
    try:
        structure = parser.get_structure("structure", str(cif_path))
        seen = set()
        ordered_chains: list[str] = []
        for chain in structure.get_chains():
            if chain.id not in seen:
                ordered_chains.append(chain.id)
                seen.add(chain.id)

        if not ordered_chains:
            raise ValueError("No chains parsed from CIF.")

        return ordered_chains
    except Exception as e:
        print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        print(f"FATAL ERROR: Biopython failed to parse chain IDs from the CIF file: {cif_path}")
        print(f"Encountered error: {e}")
        print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        raise

def get_plddt_from_cif(cif_path: Path) -> list[float]:
    """Parses a CIF file to get the pLDDT score for each C-alpha atom."""
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("structure", cif_path)
    plddts = []
    for atom in structure.get_atoms():
        if atom.name == 'CA':
            plddts.append(atom.get_bfactor())
    return plddts


def prepare_chai_session_data(input_dir: Path, model_index: int) -> dict:
    """
    Loads and processes data from a Chai model output directory.
    """
    print(f"Processing Chai model output from: {input_dir} (model index: {model_index})")

    # Define file paths
    scores_path = input_dir / f"scores.model_idx_{model_index}.npz"
    pae_path = input_dir / f"pae.model_idx_{model_index}.npz"
    structure_path = input_dir / f"pred.model_idx_{model_index}.cif"

    # Verify all files exist
    for p in [scores_path, pae_path, structure_path]:
        if not p.exists():
            raise FileNotFoundError(f"Required Chai input file not found: {p}")

    # 1. Load scores from .npz files
    scores_data = np.load(scores_path)
    pae_data = np.load(pae_path)
    pae = pae_data['pae'].tolist()
    ptm = float(scores_data['ptm'].item())
    iptm = float(scores_data['iptm'].item())

    # 2. Parse pLDDT and chains from the .cif file using Biopython
    plddt_scores = get_plddt_from_cif(structure_path)
    chain_labels = ";".join(get_chains_from_cif(structure_path))
    mean_plddt = np.mean(plddt_scores) if plddt_scores else 0

    # 3. Create the final scores JSON object
    final_scores = {
        "pae": pae,
        "plddt": plddt_scores,
        "ptm": ptm,
        "iptm": iptm,
        "meanPlddt": mean_plddt,
        "max_pae": np.max(pae) if pae else 0,
    }

    # 4. Prepare the session data dictionary, converting CIF to PDB for stability
    pdb_content = convert_cif_to_pdb_string(structure_path)
    session_data = {
        "structureFile": {
            "name": structure_path.with_suffix(".pdb").name,
            "content": pdb_content,
        },
        "scoresFile": {
            "name": "scores.json",
            "content": json.dumps(final_scores)
        },
        "chainLabels": chain_labels,
        "crosslinksFile": None
    }
    return session_data

def prepare_boltz_session_data(input_dir: Path, model_name: str, model_index: int) -> dict:
    """
    Loads and processes data from a Boltz model output directory.
    """
    print(f"Processing Boltz model output from: {input_dir} (model: {model_name}, index: {model_index})")

    # Define file paths
    base_name = f"{model_name}_model_{model_index}"
    pred_dir = input_dir / model_name
    scores_path = pred_dir / f"confidence_{base_name}.json"
    pae_path = pred_dir / f"pae_{base_name}.npz"
    plddt_path = pred_dir / f"plddt_{base_name}.npz"
    structure_path = pred_dir / f"{base_name}.cif"

    # Verify all files exist
    for p in [scores_path, pae_path, plddt_path, structure_path]:
        if not p.exists():
            raise FileNotFoundError(f"Required Boltz input file not found: {p}")

    # 1. Load scores
    with open(scores_path, "r") as f:
        scores_data = json.load(f)
    pae_data = np.load(pae_path)
    plddt_data = np.load(plddt_path)

    pae = pae_data['pae'].tolist()
    plddt_scores_raw = plddt_data['plddt']
    
    # Per analysis, Boltz pLDDT scores are 0-1, but the viewer expects 0-100.
    # We must scale them before passing them to the viewer.
    plddt_scores = (plddt_scores_raw * 100).tolist()

    ptm = scores_data['ptm']
    iptm = scores_data['iptm']
    
    # Still need chains from CIF, now with robust parsing
    chain_labels = ";".join(get_chains_from_cif(structure_path))
    mean_plddt = np.mean(plddt_scores) if plddt_scores else 0

    # 2. Create the final scores JSON object
    final_scores = {
        "pae": pae,
        "plddt": plddt_scores,
        "ptm": ptm,
        "iptm": iptm,
        "meanPlddt": mean_plddt,
        "max_pae": np.max(pae) if pae else 0,
    }

    # 3. Prepare the session data dictionary, converting CIF to PDB for stability
    pdb_content = convert_cif_to_pdb_string(structure_path)
    session_data = {
        "structureFile": {
            "name": structure_path.with_suffix(".pdb").name,
            "content": pdb_content,
        },
        "scoresFile": {
            "name": "scores.json",
            "content": json.dumps(final_scores)
        },
        "chainLabels": chain_labels,
        "crosslinksFile": None
    }
    return session_data


def get_session_data(
    structure_path: Path,
    chain_labels: str,
    scores_path: Optional[Path] = None,
    crosslinks_path: Optional[Path] = None,
) -> dict[str, Any]:
    """Prepare the session data dictionary from input files."""
    data = {
        "structureFile": {
            "name": structure_path.name,
            "content": structure_path.read_text(),
        },
        "chainLabels": chain_labels,
    }
    if scores_path:
        data["scoresFile"] = {
            "name": scores_path.name,
            "content": scores_path.read_text(),
        }
    if crosslinks_path:
        data["crosslinksFile"] = {
            "name": crosslinks_path.name,
            "content": crosslinks_path.read_text(),
        }
    return data


def create_standalone_html(
    project_root: Path,
    session_data: dict,
    output_handle: str
):
    """Assemble the final standalone HTML file."""
    # 1. Read the main template
    with open(project_root / "index.html", "r", encoding="utf-8") as f:
        html_template = f.read()

    # 2. Define asset paths
    css_files = sorted((project_root / "libs").glob("*.css")) + \
                sorted((project_root / "src" / "css").glob("*.css"))
    
    js_libs = [
        project_root / "libs" / "jquery.min.js",
        project_root / "libs" / "jquery-ui.min.js",
        project_root / "libs" / "bootstrap.bundle.min.js",
        project_root / "libs" / "ngl.js",
        project_root / "libs" / "chroma.min.js",
        project_root / "libs" / "FileSaver.min.js",
    ]
    
    js_app_entrypoint = project_root / "src" / "js" / "setup.js"
    templates_dir = project_root / "src" / "templates"
    
    # 3. Process and inline all assets
    print("Inlining CSS...")
    inlined_css = inline_css(project_root, css_files)
    
    print("Inlining HTML templates...")
    inlined_templates = inline_html_templates(templates_dir)
    
    print("Inlining JS libraries...")
    inlined_js_libs = inline_js_libs(project_root, js_libs)
    
    print("Bundling application JavaScript with esbuild...")
    bundled_js_app = bundle_js_app(project_root, js_app_entrypoint)

    # 4. Embed data
    print(f"[DEBUG] Final chain labels used in session_data: {session_data.get('chainLabels')}")
    session_json = json.dumps(session_data, indent=4)
    session_script = f'<script type="application/json" id="session-data">\n{session_json}\n</script>'

    # Per user feedback, escape literal closing tags in *both* inline blocks
    # to prevent premature HTML parsing errors.
    bundled_js_app  = bundled_js_app.replace("</script>", "<\\/script>")

    # 5. Assemble the final HTML
    # Remove original script and link tags
    final_html = re.sub(r'<link rel="stylesheet".*?>', '', html_template, flags=re.DOTALL)
    # Per user feedback, replace all icon links with a single data URI to prevent 404s
    final_html = re.sub(r'<link rel="icon".*?>', '', final_html, flags=re.DOTALL)
    final_html = re.sub(r'<script defer .*?></script>', '', final_html)
    final_html = re.sub(r'<script type="text/javascript".*?></script>', '', final_html)

    # Inject inlined content
    final_html = final_html.replace("</head>",
        '<link rel="icon" href="data:,">\n'
        f'<style>{inlined_css}</style>\n'
        f'{session_script}\n'
        f'</head>'
    )
    
    final_html = final_html.replace("</body>",
        f'{inlined_templates}\n'
        f'<script>{inlined_js_libs}</script>\n'
        f'<script>{bundled_js_app}</script>\n'
        f'</body>'
    )

    # Add a script to notify the parent window that it has loaded successfully
    messaging_script = """
    <script>
        window.addEventListener('load', () => {
            if (window.opener) {
                window.opener.postMessage('pae_viewer_loaded', '*');
            }
        });
    </script>
    """
    final_html = final_html.replace("</body>", f"{messaging_script}</body>")

    # 6. Save the output file
    timestamp = datetime.now().strftime("%Y%m%d%H%M%S")
    output_filename = f"{output_handle}_report_{timestamp}.html"
    output_path = project_root / output_filename
    
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(final_html)
        
    print(f"\nSuccessfully created standalone report: {output_path}")
    return output_path


def main():
    """Main function to parse arguments and run the export."""
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter,
        conflict_handler='resolve'  # Allows overriding args
    )

    # Standard workflow arguments
    parser.add_argument(
        "-s", "--structure", help="Path to structure file (e.g., PDB, CIF)",
        type=resolved_path
    )
    parser.add_argument(
        "-r", "--scores", help="Path to JSON file containing PAE scores",
        type=resolved_path
    )
    parser.add_argument(
        "-l", "--labels", help="Semicolon-separated list of chain labels",
        type=str
    )
    parser.add_argument(
        "-c", "--crosslinks", help="Path to CSV/TSV containing crosslinks",
        type=resolved_path
    )
    parser.add_argument(
        "-o", "--output", help="Handle for the output file name (e.g., 'MyComplex')",
        type=str, default="pae_viewer"
    )
    parser.add_argument(
        "--open", help="Open the generated report in the browser",
        action="store_true"
    )

    # Arguments for Chai model processing
    parser.add_argument(
        "--chai-input-dir",
        help="Path to a Chai model output directory. Using this option ignores other data inputs.",
        type=Path
    )
    parser.add_argument(
        "--model-index",
        help="The model index to process from the Chai directory (default: 0).",
        type=int,
        default=0
    )

    # Arguments for Boltz model processing
    parser.add_argument(
        "--boltz-input-dir",
        help="Path to a Boltz model 'predictions' directory.",
        type=Path
    )
    parser.add_argument(
        "--boltz-model-name",
        help="The model name prefix used in Boltz output files (e.g., 'Hemoglobin_tetramer_boltz_inference_generated_with_msas').",
        type=str
    )

    args = parser.parse_args()

    # Manually check for required arguments based on the workflow
    if args.chai_input_dir:
        # Chai workflow doesn't need the others
        pass
    elif args.boltz_input_dir:
        if not args.boltz_model_name:
            parser.error("--boltz-model-name is required when using --boltz-input-dir")
    elif not all([args.structure, args.scores, args.labels]):
        parser.error("For the standard workflow, the following arguments are required: -s/--structure, -r/--scores, -l/--labels")


    project_root = Path(__file__).parent.resolve().parents[1]

    if args.chai_input_dir:
        session_data = prepare_chai_session_data(args.chai_input_dir, args.model_index)
        output_handle = f"chai_report_model_{args.model_index}"
    elif args.boltz_input_dir:
        session_data = prepare_boltz_session_data(args.boltz_input_dir, args.boltz_model_name, args.model_index)
        output_handle = f"boltz_report_{args.boltz_model_name}_model_{args.model_index}"
    else:
        # This block is now safe because we've already validated the args
        session_data = get_session_data(
            args.structure, args.labels, args.scores, args.crosslinks
        )
        output_handle = args.output

    # ------------------------------------------------------------------
    # Optional override of chain labels
    # If the user supplied --labels we should respect it for ALL workflows,
    # even when using the specialized Chai / Boltz helpers which otherwise
    # infer labels alphabetically from the structure. This allows callers
    # (like OmniFold) to define an explicit ordering that matches the PAE
    # matrix.
    # ------------------------------------------------------------------
    if args.labels:
        print(f"Overriding chain labels via --labels: {args.labels}")
        session_data["chainLabels"] = args.labels

    output_path = create_standalone_html(project_root, session_data, output_handle)

    if args.open:
        webbrowser.open(output_path.as_uri())


if __name__ == "__main__":
    main() 