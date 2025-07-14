<div class="accordion-item">
  <h2 class="accordion-header" id="helpPaeViewerHeading">
    <button class="accordion-button collapsed" type="button"
            data-bs-toggle="collapse" data-bs-target="#helpPaeViewer"
            aria-expanded="true" aria-controls="helpPaeViewer">
      <span class="fw-bold">PAE Viewer</span>
    </button>
  </h2>
  <div id="helpPaeViewer" class="accordion-collapse collapse"
       aria-labelledby="helpPaeViewerHeading"
       data-bs-parent="#help">
    <div class="accordion-body">
      <p>
        The predicted aligned error (PAE) is a metric for the confidence
        in relative positions of the predicted structure.
      </p>
      <p>
        A high PAE at position (x, y) in the PAE matrix indicates
        a high expected position error for the residue at x, if
        the predicted and true structure are aligned at residue y.
        An extensive tutorial for the PAE can be found on the
        <a href="https://alphafold.ebi.ac.uk/entry/Q5VSL9" target="_blank">
          AlphaFold Protein Structure Database
        </a>
        page bottom for any entry.
      </p>
      <p>
        In contrast to the viewer of the AlphaFold Protein Structure
        Database, this PAE viewer allows to view the PAE of multimers,
        and integrates visualization of crosslink data, as well.
        The latter can be an important indicator for the reliability
        of the structure prediction.
      </p>
      <div class="help-text">
        <img class="m-3" src="src/img/pae-demo/screenshot-pae-viewer.png"
             style="width: 400px" alt="Screenshot of the PAE viewer.">
        <div>
          <p>
            As can be seen on the left, the PAE matrix is divided up into
            several sections corresponding to the different subunits,
            whose labels can be found on the axes, as well. The axes ticks
            indicate the position of the residues within the subunits.
            The ticks at the end of the divider lines display the length
            of the subunit's amino acid sequence.
          </p>
          <p>
            The circular markers correspond to cross-linked residues.
            They are colored to indicate violation (red) or satisfaction
            (blue) of crosslinker length restraints. In the case of the
            examples, a Cα-Cα distance ≥ 30 Å is considered a restraint
            violation. When a marker is clicked, the corresponding
            crosslink in the structure viewer is highlighted.
          </p>
        </div>
      </div>
      <div class="help-columns">
        <div class="card" style="width: 250px;">
          <img class="card-img-top"
               src="src/img/pae-demo/screenshot-pae-viewer-point-selection.png"
               alt="Screenshot of the PAE viewer with a point selected.">
          <div class="p-3">
            <p>
              If some point of the PAE matrix is clicked, the position is
              marked and the corresponding x and y coordinates are
              projected onto the diagonal, which helps to interpret
              the relative distance of the residues within a sequence.
              The projected x coordinate is color-coded
              <span class="sp-color-marker" style="background-color: cyan; color: black">cyan</span>,
              and y is color-coded
              <span class="sp-color-marker" style="background-color: orange; color: black">orange</span>.
              In the sequence viewer, the residue pair corresponding
              to that point is highlighted. For the structure viewer,
              the distance between the residues is also visualized.
            </p>
          </div>
        </div>
        <div class="card" style="width: 250px;">
          <img class="card-img-top"
               src="src/img/pae-demo/screenshot-pae-viewer-range-selection.png"
               alt="Screenshot of the PAE viewer with a rectangle selected.">
          <div class="p-3">
            <p>
              If a rectangle is selected by holding a click, the
              corresponding sequence ranges are projected onto
              the diagonal, using the same color code.
              Again, the corresponding residues are highlighted in
              the sequence viewer as well as the structure viewer.
            </p>
          </div>
        </div>
        <div class="card" style="width: 250px;">
          <img class="card-img-top"
               src="src/img/pae-demo/screenshot-pae-viewer-overlapping-selection.png"
               alt="Screenshot of the PAE viewer with a rectangle selected. The lines projected onto the diagonal are overlapping.">
          <div class="p-3">
            <p>
              If the projected x and y ranges are overlapping,
              the overlap is color-coded
              <span class="sp-color-marker" style="background-color: magenta">magenta</span>.
              A rectangular selection of this nature is hard to interpret
              in terms of the PAE of the relative position. For this
              reason, this special color code was introduced.
            </p>
          </div>
        </div>
        <div class="card" style="width: 250px;">
          <img class="card-img-top"
               src="src/img/pae-demo/screenshot-pae-viewer-regions.png"
               alt="Screenshot of the PAE viewer with highlighted sections corresponding to the protein subunits.">
          <div class="p-3">
            <p>
              By holding the Shift key while the cursor is hovering
              over the PAE matrix, or by using the checkbox below
              the matrix, an overlay can be toggled, which shows
              rectangular sections of the matrix which correspond
              to the PAE of the complex subunits and of their
              positions relative to each other. If clicked, the
              corresponding residues are highlighted in structure
              viewer, using the color code of the subunits.
            </p>
          </div>
        </div>
      </div>
      <p class="mt-4">
        For all selections the corresponding (mean) PAE value is displayed
        numerically under the PAE graph. In case of crosslink selections,
        the mean value of the PAEs for both orientations (x/y and y/x) is
        displayed.
      </p>
    </div>
  </div>
</div>