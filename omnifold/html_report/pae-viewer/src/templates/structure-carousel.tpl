<template id="structure-carousel-template">
  <div class="structure-carousel">
    <div class="sc-container">
      <div class="sc-box">
        <div class="sc-viewer">
          <div class="carousel slide"
               data-bs-interval="false" data-bs-touch="false">
            <select class="sc-options" autocomplete="off">
              <option selected value="index">Default style</option>
              <option value="electrostatic">Electrostatic</option>
              <option value="hydrophobicity">Hydrophobicity</option>
            </select>
            <div class="sc-fullscreen-toggle sc-transparent-box">
              <img class="sc-enter-fullscreen"
                   src="src/img/bootstrap-icons/arrows-fullscreen.svg">
              <img class="sc-exit-fullscreen hidden"
                   src="src/img/bootstrap-icons/fullscreen-exit.svg">
            </div>
            <div class="sc-inner"></div>
            <button class="sc-control-prev hidden" type="button"
                    data-bs-slide="prev">
              <span class="sc-control-prev-icon" aria-hidden="true"></span>
              <span class="visually-hidden">Previous</span>
            </button>
            <button class="sc-control-next hidden" type="button"
                    data-bs-slide="next">
              <span class="sc-control-next-icon" aria-hidden="true"></span>
              <span class="visually-hidden">Next</span>
            </button>
            <div class="sc-info sc-transparent-box"
                 data-bs-toggle="tooltip" data-bs-placement="bottom"
                 data-bs-trigger="click" data-bs-html="true"
                 title="<a href='https://nglviewer.org/' target='_blank'>NGL Viewer</a> is used for molecular visualization.">
              <img src="src/img/bootstrap-icons/question-circle.svg">
            </div>
          </div>
          <div class="sc-legend-scale hidden">
            <span class="sc-legend-caption caption-left">polar</span>
            <span class="sc-legend-caption caption-right">unpolar</span>
          </div>
        </div>
        <div class="sc-panel"></div>
      </div>
    </div>
    <div class="sc-container sc-placeholder hidden"></div>
  </div>
</template>
