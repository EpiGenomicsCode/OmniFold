<template id="structure-panel-demo-template">
  <div class="structure-panel">
    <div class="sp-description hidden"></div>
    <div class="sp-complex hidden">
      <div class="hidden">
        <strong>Source</strong>:
        <a href="http://doi.org/10.15252/msb.202311544 " target="_blank">O'Reilly <i>et al.</i> (2023)</a>
      </div>
      <div class="sp-chain-legend"></div>
      <div class="sp-metrics">
        <strong>Confidence metrics:</strong>
        <ul>
          <li class="sp-metrics-plddt">
            <span>
            Mean pLDDT<sup><a
              href="https://doi.org/10.1038/s41586-021-03819-2" target="_blank"
              title="Highly accurate protein structure prediction with AlphaFold (2021)">?</a></sup>:</span>
            <span class="sp-metric-value" data-type="plddt"></span>
          </li>
          <li class="sp-metrics-ptm">
            <span>pTM<sup><a
              href="https://doi.org/10.1038/s41586-021-03819-2" target="_blank"
              title="Highly accurate protein structure prediction with AlphaFold (2021)">?</a></sup>:</span>
            <span class="sp-metric-value" data-type="ptm"></span>
          </li>
          <li class="sp-metrics-iptm">
            <span>ipTM<sup><a
              href="https://doi.org/10.1101/2021.10.04.463034"
              target="_blank" title="Protein complex prediction with AlphaFold-Multimer (2022)">?</a></sup>:</span>
            <span class="sp-metric-value" data-type="iptm"></span></li>
        </ul>
      </div>
      <div class="sp-color-options" data-scheme="structure">
        <strong>Structure color scheme:</strong>
        <div>
          <label class="sp-centered-label">
            <input type="radio" name="sp-scheme" value="subunit" checked>
            by subunit
          </label>
        </div>
        <div>
          <label class="sp-centered-label">
            <input type="radio" name="sp-scheme" value="confidence">
            model confidence (pLDDT)
          </label>
        </div>
        <div class="sp-confidence-legend hidden">
          <div class="sp-centered-label sp-legend-label">
            <span class="sp-color-square" style="background-color: #0053D6"></span>
            Very high (pLDDT > 90)
          </div>
          <div class="sp-centered-label sp-legend-label">
            <span class="sp-color-square" style="background-color: #64CBF3"></span>
            Confident (90 > pLDDT > 70)
          </div>
          <div class="sp-centered-label sp-legend-label">
            <span class="sp-color-square" style="background-color: #FFDB12"></span>
            Low (70 > pLDDT > 50)
          </div>
          <div class="sp-centered-label sp-legend-label">
            <span class="sp-color-square" style="background-color: #FF7D45"></span>
            Very low (pLDDT < 50)
          </div>
        </div>
      </div>
      <div class="sp-color-options hidden" data-scheme="crosslinks">
        <strong>Crosslink color scheme:</strong>
        <div>
          <label class="sp-centered-label">
            <input type="radio" name="crosslink-scheme"
                   value="restraint" checked>
            <span class="sp-color-marker" style="background-color: blue">
            satisfied
          </span>
            /
            <span class="sp-color-marker" style="background-color: red">
            violated
          </span>
            restraint
          </label>
        </div>
        <div>
          <label class="sp-centered-label">
            <input type="radio" name="crosslink-scheme" value="connection">
            <span class="sp-color-marker" style="background-color: grey">
            intra
          </span>
            /
            <span class="sp-color-marker" style="background-color: orange">
            inter
          </span>
            protein links
          </label>
        </div>
        <div>
          <label class="sp-centered-label">
            <input type="radio" name="crosslink-scheme" value="">
            hide all crosslinks
          </label>
        </div>
        <p class="sp-restraint-note hidden">
          <small>
            <strong>Note:</strong>
            a Cα-Cα distance ≥ 30 Å is considered a restraint violation.
          </small>
        </p>
      </div>
      <div class="hidden">
        <a class="sp-crosslinks-download"
           href="" download>
          <img src="src/img/bootstrap-icons/download.svg"
               title="Download complex data">
          <span>
          Download structure, crosslink data and PyMOL/ChimeraX files
        </span>
        </a>
      </div>
    </div>
  </div>

</template>
