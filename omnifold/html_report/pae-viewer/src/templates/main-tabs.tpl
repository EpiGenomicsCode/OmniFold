<div class="card mb-3">
  <div class="card-header">
    <ul class="nav nav-tabs card-header-tabs" id="tabs" role="tablist">
      <li class="nav-item" role="presentation">
        <button class="nav-link active" id="examples-tab" data-bs-toggle="tab"
                data-bs-target="#examples" type="button" role="tab"
                aria-controls="examples" aria-selected="true">
          Examples
        </button>
      </li>
      <li class="nav-item" role="presentation">
        <button class="nav-link" id="upload-tab" data-bs-toggle="tab"
                data-bs-target="#upload" type="button" role="tab"
                aria-controls="upload" aria-selected="false">
          Upload
        </button>
      </li>
      <li class="nav-item" role="presentation">
        <button class="nav-link" id="offline-tab" data-bs-toggle="tab"
                data-bs-target="#offline" type="button" role="tab"
                aria-controls="offline" aria-selected="false">
          Offline version
        </button>
      </li>
      <li class="nav-item" role="presentation">
        <button class="nav-link" id="citation-tab" data-bs-toggle="tab"
                data-bs-target="#citation" type="button" role="tab"
                aria-controls="citation" aria-selected="false">
          Citation
        </button>
      </li>
    </ul>
  </div>

  <div class="tab-content p-3" id="tabContents">
    {{main-tab-examples.tpl}}
    {{main-tab-upload.tpl}}
    {{main-tab-offline.tpl}}
    {{main-tab-citation.tpl}}
  </div>
</div>
