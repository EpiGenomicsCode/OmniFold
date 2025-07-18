{{main-consent.tpl}}

<div class="toast-pae-max-warning toast align-items-center text-bg-warning border-0 position-fixed m-3 top-0 start-50 translate-middle-x"
     role="alert" aria-live="assertive" aria-atomic="true">
  <div class="d-flex">
    <div class="toast-body"></div>
    <button type="button" class="btn-close btn-close-white me-2 m-auto" data-bs-dismiss="toast" aria-label="Close"></button>
  </div>
</div>

<div class="toast-modification-warning toast align-items-center text-bg-warning border-0 position-fixed m-3 top-0 start-50 translate-middle-x"
     role="alert" aria-live="assertive" aria-atomic="true">
  <div class="d-flex">
    <div class="toast-body">
      The atom-wise PAE values of modified amino acids were replaced
      by their arithmetic mean to get a residue-wise PAE (note the
      different PAE matrix).
    </div>
    <button type="button" class="btn-close btn-close-white me-2 m-auto" data-bs-dismiss="toast" aria-label="Close"></button>
  </div>
</div>

{{main-news.tpl}}
{{main-compatibility.tpl}}
{{main-tabs.tpl}}

<div id="complex-sequence-viewer"></div>

<div id="split-view">
  <div id="test-carousel"></div>
  <div id="pae-viewer"></div>
</div>

{{main-help.tpl}}

{{structure-carousel.tpl}}
{{structure-viewer.tpl}}
{{structure-panel-demo.tpl}}
{{pae-viewer.tpl}}
{{complex-sequence-viewer.tpl}}
