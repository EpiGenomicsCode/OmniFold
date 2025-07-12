<template id="pae-viewer-template">
  <div class="pae-viewer">
    <svg class="pv-graph"
         xmlns="http://www.w3.org/2000/svg"
         xmlns:xlink="http://www.w3.org/1999/xlink"
         overflow="visible"
    >
      <defs>
        <linearGradient id="rangeMarkerLower"
                        gradientTransform="rotate(45 0.5 0.5)">
          <stop offset="0%"  stop-color="cyan" />
          <stop offset="50%" stop-color="cyan" />
          <stop offset="50%" stop-color="magenta" />
          <stop offset="100%" stop-color="magenta" />
        </linearGradient>
        <linearGradient id="rangeMarkerLowerReversed"
                        gradientTransform="rotate(45 0.5 0.5)">
          <stop offset="0%"  stop-color="orange" />
          <stop offset="50%" stop-color="orange" />
          <stop offset="50%" stop-color="magenta" />
          <stop offset="100%" stop-color="magenta" />
        </linearGradient>
        <linearGradient id="rangeMarkerUpper"
                        gradientTransform="rotate(45 0.5 0.5)">
          <stop offset="0%"  stop-color="magenta" />
          <stop offset="50%" stop-color="magenta" />
          <stop offset="50%" stop-color="orange" />
          <stop offset="100%" stop-color="orange" />
        </linearGradient>
        <linearGradient id="rangeMarkerUpperReversed"
                        gradientTransform="rotate(45 0.5 0.5)">
          <stop offset="0%"  stop-color="magenta" />
          <stop offset="50%" stop-color="magenta" />
          <stop offset="50%" stop-color="cyan" />
          <stop offset="100%" stop-color="cyan" />
        </linearGradient>
        <pattern id="stripes-template" patternUnits="userSpaceOnUse"
                 width="2%" height="2%" patternTransform="rotate(45)">
          <rect x="0" y="0" width="2%" height="2%" fill="red"></rect>
          <line x1="0" y1="0" x2="0" y2="2%"
                stroke="#00FF00" stroke-width="2%" />
        </pattern>
      </defs>

      <g class="pv-graph-area">
        <image
          class="pv-pae-matrix"
          width="100%"
          height="100%"
          style="image-rendering:pixelated"
        />
        <g class="pv-axes">
          <text class="pv-axis-x-label" x="50%" y="-25%" text-anchor="middle"
                font-family="Arial, Helvetica, sans-serif">
            Scored residue / atom
          </text>
          <text class="pv-axis-y-label" x="-35%" y="50%" text-anchor="middle"
                font-family="Arial, Helvetica, sans-serif"
                transform-origin="-35% 50%" transform="rotate(-90)">
            Aligned residue / atom
          </text>
          <line class="pv-axis-x" x1="-5%" y1="0" x2="105%" y2="0"/>
          <line class="pv-axis-y" x1="0" y1="-5%" x2="0" y2="105%"/>
          <line class="pv-axis-diagonal" x1="0" y1="0" x2="102%" y2="102%"/>
          <g class="pv-sequence-ticks"></g>
          <g class="pv-unit-ticks"></g>
          <g class="pv-unit-tick-labels"></g>
          <g class="pv-sequence-tick-labels"></g>
        </g>
        <g class="pv-dividers"></g>
        <g class="pv-interactive-layer">
          <rect class="pv-backdrop" x="0" y="0"
                width="100%" height="100%" opacity="0" />
          <g class="pv-selections"></g>
          <g class="pv-crosslinks"></g>
          <g class="pv-regions"></g>
        </g>
      </g>
    </svg>
    <div class="pv-color-legend">
      <div class="pv-color-scale"></div>
      <div class="pv-color-ticks">
        <span style="position: relative">&nbsp;</span> <!-- placeholder -->
        <span style="left: 0%">0</span>
      </div>
      <div class="pv-color-legend-label">
        Expected position error (Ångströms)
      </div>
    </div>
    <div class="pv-status">
      <div>
        <div><b>Cursor</b>:</div>
        <div class="pv-status-cursor"></div>
      </div>
      <div>
        <div><b>Selection</b>:</div>
        <div class="pv-status-selection"></div>
      </div>
    </div>
    <div class="pv-region-toggle">
      <div>
        <label>
          <input type="checkbox" name="toggleRegion">
          toggle region overlay
        </label>
      </div>
      <div>
        <small>
          (alternatively, hold the 'Shift' key when the cursor is within
          the PAE Viewer)
        </small>
      </div>
    </div>
    <div class="d-flex flex-column align-items-center">
      <div class="m-2">
        <a class="pv-download-matrix btn btn-primary hidden" download>
          Download PAE matrix image
        </a>
      </div>
      <div class="m-2">
        <a class="pv-download-svg btn btn-primary" download>
          Download graph as SVG
        </a>
      </div>
    </div>
  </div>
</template>
