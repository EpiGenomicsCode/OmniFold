<template id="structure-viewer-template">
  <div class="structure-viewer">
    <div class="sv-loading">
      <!-- empty divs necessary for CSS animation, don't delete! -->
      <div></div>
      <div></div>
      <div></div>
      <div></div>
    </div>
    <div class="sv-error hidden">
      <img src="src/img/bootstrap-icons/exclamation-circle.svg">
      <p class="sv-error-message"></p>
    </div>
    <div class="sv-display"></div>
    <div class="sv-caption hidden">
      <a class="sv-transparent-box" href="" target="_blank"></a>
    </div>
    <div class="sv-screenshot sv-transparent-box">
      <input type="image" src="src/img/bootstrap-icons/camera.svg" />
    </div>
  </div>
</template>
