<div class="tab-pane fade show" id="upload" role="tabpanel" aria-labelledby="upload-tab">
  <form id="uploadStructure" class="needs-validation" novalidate>
    <div class="my-3">
      This form can be used to visualize the output structure and PAE of
      <a href="https://doi.org/10.1038/s41586-024-07487-w" target="_blank">
        AlphaFold 3</a>
      (use
      <a href="https://alphafoldserver.com/" target="_blank">
        AlphaFold Server</a>
      for an online version)
      and
      <a href="https://www.biorxiv.org/content/10.1101/2021.10.04.463034v2"
         target="_blank">
        AlphaFold-Multimer</a>
      (use
      <a href="https://colab.research.google.com/github/sokrypton/ColabFold/blob/v1.5.2/AlphaFold2.ipynb"
         target="_blank">
        ColabFold</a>
      or the
      <a href="https://colab.research.google.com/github/deepmind/alphafold/blob/main/notebooks/AlphaFold.ipynb"
         target="_blank">
        Colab notebook from Deepmind</a>
      for an online version), as well as crosslink data.
      You can also
      <a href="resources/sample/GatA-GatB.zip" download="GatA-GatB.zip">
        download some sample files (GatA-GatB.zip)</a>.
      Please have a look at the included ReadMe and the 'Help' section at the
      page bottom for instructions.
    </div>

    <div class="my-3">
      <label class="form-label mb-2" for="structureFile">
        Structure file
      </label>
      <input type="file" id="structureFile" class="form-control"
             name="structureFile" accept=".pdb,.mmcif,.cif,.mcif" required />
      <div id="structureFileExample" class="form-text">
        Structure file
        (<span class="font-monospace">.pdb,.mmcif,.cif,.mcif</span>),
        e.g.
        <span class="font-monospace">fold_*_model_?.cif</span>
        or
        <span class="font-monospace">*_unrelaxed_rank_?_model_?.pdb</span>.
      </div>
      <div id="structureFileFeedback" class="invalid-feedback"></div>
    </div>
    <div class="my-3">
      <label class="form-label mb-2" for="chainLabels">
        Chain labels (optional)
      </label>
      <div class="has-validation">
        <input type="text" id="chainLabels" class="form-control"
               name="chainLabels" />
        <div class="form-text">
          You can provide labels for the chains in the order they appear
          in the structure file by using a semicolon-separated list, e.g.
          <span class="font-monospace">GatA;GatB</span>.
        </div>
        <div id="chainLabelsFeedback" class="invalid-feedback"></div>
      </div>
    </div>
    <div class="my-3">
      <label class="form-label mb-2" for="scoresFile">Scores file</label>
      <input type="file" id="scoresFile" class="form-control"
             name="scoresFile" accept=".json" required />
      <div class="form-text">
        JSON file containing PAE, e.g.
        <span class="font-monospace">fold_*_full_data_?.json</span>
        or
        <span class="font-monospace">*_unrelaxed_rank_?_model_?_scores.json</span>.
        Use the conversion script to convert the
        <span class="font-monospace">.pkl</span> output of AlphaFold-Multimer
        to the required <span class="font-monospace">.json</span> format:
        <span class="font-monospace">
          <a href="resources/scripts/jsonify_scores.py" download>jsonify_scores.py</a>
        </span>
        (see 'Help' section for details).
      </div>
      <div id="scoresFileFeedback" class="invalid-feedback"></div>
    </div>
    <div class="my-3">
      <label class="form-label mb-2" for="scoresFile">Crosslinks file (optional)</label>
      <input type="file" id="crosslinksFile" class="form-control"
             name="crosslinksFile" accept=".csv,.pb" />
      <div class="form-text">
        CSV file containing crosslink data, e.g.
        <span class="font-monospace">GatA-GatB.csv</span>.
        Alternatively, a <span class="font-monospace">.pb</span> (pseudobond)
        file can be uploaded (see 'Help' section for limitations).
      </div>
      <div id="crosslinksFileFeedback" class="invalid-feedback"></div>
    </div>

    <div>
      <input id="uploadStructureSubmit" class="btn btn-primary float-end" type="submit" value="Upload" disabled/>
    </div>
  </form>
</div>
