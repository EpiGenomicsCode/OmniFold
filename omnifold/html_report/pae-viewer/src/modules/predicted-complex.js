import {PaeViewer} from "./pae-viewer.js";
import {StructureCarousel} from "./structure-carousel.js";
import {ComplexSequenceViewer} from "./complex-sequence-viewer.js";

const id = new URLSearchParams(window.location.search).get('id');

const structureCarousel = new StructureCarousel(
    document.querySelector('#structure-carousel'),
    document.querySelector('#structure-carousel-template'),
    document.querySelector('#structure-viewer-template'),
    document.querySelector('#structure-panel-template')
);

const toggleFullscreenBox = structureCarousel.getElement().querySelector(
    '.sc-fullscreen-toggle'
);

const paeViewer = new PaeViewer(
    document.querySelector('#pae-viewer'),
    document.querySelector('#pae-viewer-template')
);

paeViewer.fetch(
    `examples/${id}/meta.json`, {'headers': {'Accept': 'application/json'}}
);

structureCarousel.addStructure({
    id: id,
    type: "Predicted Complex",
    source: "predicted-complex",
    description: "Complex predicted with AlphaFold-Multimer."
}, true).then(entry => {
    document.querySelector('.sp-complex').open = true;

    document.addEventListener('mousedown', event => {
        if (
            !entry.item.contains(event.target)
            && !toggleFullscreenBox.contains(event.target)
        ) {
            paeViewer.deselectAll(true);
        }
    });
});

const complexSequenceViewer = new ComplexSequenceViewer(
    document.querySelector('#complex-sequence-viewer'),
    document.querySelector('#complex-sequence-viewer-template')
);

complexSequenceViewer.fetch(
    `examples/${id}/meta.json`, {'headers': {'Accept': 'application/json'}}
);
