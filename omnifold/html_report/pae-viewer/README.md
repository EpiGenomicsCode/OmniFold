# PAE Viewer

This repository contains the source code for the main components of
the [PAE Viewer webserver](https://pae-viewer.uni-goettingen.de/).

Released under MIT license.

## Offline usage
If you would like to use the PAE Viewer offline in your browser and
start it with a command line interface, you can use the scripts
provided by the repository
[general-microbiology/pae-viewer](https://gitlab.gwdg.de/general-microbiology/pae-viewer). 
Python >=3.9 is required to run the offline version.

### Instructions
1. Download the PAE Viewer project files
   ([pae-viewer-main.zip](https://gitlab.gwdg.de/general-microbiology/pae-viewer/-/archive/main/pae-viewer-main.zip)) 
   and extract the archive.

2. Using the terminal, change the current working directory to the project root
   (`/your/download/path/pae-viewer-main`) and start a local HTTP server with Python.

   ```bash
   cd /your/download/path/pae-viewer-main
   python -m http.server 8000
   ``` 
3. Open http://localhost:8000/ in your browser to 
   start the offline version.

4. To start PAE Viewer in your browser with arguments supplied via CLI, run the
   `pae-viewer-main/resources/pae_viewer.py` script from the project directory. 
   For example:

   ```bash
   cd /your/download/path/pae-viewer-main
   python resources/scripts/pae_viewer.py \
     --structure resources/sample/GatA-GatB/fold_gata_gatb_model_0.cif \
     --labels 'GatA;GatB' \
     --scores resources/sample/GatA-GatB/fold_gata_gatb_full_data_0.json \
     --crosslinks resources/sample/GatA-GatB/GatA-GatB.csv
   ```
 
   For this, the local HTTP server must be running, as well. The Python
   script creates a HTML session file with all data embedded and opens
   it in the browser. It basically prefills the upload form and submits
   it. The session files are stored permanently in the
   `pae-viewer-main/standalone` directory, and can be revisited at any time.


## Dependencies
- [Bootstrap 5.2.3](https://getbootstrap.com/)
- [Chroma 2.4.2](https://gka.github.io/chroma.js/)
- [NGL Viewer 2.0.1](http://nglviewer.org/#ngl)
- [FileSaver.js 2.0.4](https://github.com/eligrey/FileSaver.js#readme)
