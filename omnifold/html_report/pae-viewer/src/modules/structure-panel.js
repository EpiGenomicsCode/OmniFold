import {
    defaultStyles,
    RepresentationStyle,
    StructureStyle,
    subunitColors,
    ternaryColoring,
    trimerScheme
} from "./structure-style.js";
import {createRandomId} from "./utils.js";

export class StructurePanel {
    #root;
    #template;
    #element;
    #description;
    #crosslinks;
    #complex;

    #viewer;
    #empty;

    constructor(root, template, viewer, visible = true) {
        this.#root = root;
        this.#template = template;
        this.#element = this.#template.content.querySelector(
            '.structure-panel'
        ).cloneNode(true);

        this.#root.appendChild(this.#element);

        this.#viewer = viewer;

        const description = this.#element.querySelector('.sp-description');
        const crosslinks = this.#element.querySelector('.sp-crosslinks');
        const complex = this.#element.querySelector('.sp-complex');

        const meta = this.#viewer.getMeta();

        // workaround
        if (meta.id === "169") {
            description.innerHTML = 'Complex predicted with '
                + '<a href="https://doi.org/10.1101/2023.06.07.544059" target="_blank">AlphaLink2</a>.';
        } else if (meta.id === "171" || meta.id === "172") {
            description.innerHTML = meta.description
                + ' (source: <a href="https://pubmed.ncbi.nlm.nih.gov/38405867/" target="_blank">PubMed</a>)';
        } else {
            description.textContent = meta.description;
        }

        if (meta.complex !== null) {
            complex.classList.remove('hidden');
            const sourceInfo = complex.querySelector('.sp-complex-source');

            if (sourceInfo) {
                sourceInfo.classList.toggle(
                    'hidden', meta.complex.id > 167
                );
            }

            new PanelControls(complex, this.#viewer, meta);
        }

        this.#empty = !(meta.description /* || ... */)

        this.setVisible(visible);
    }

    setVisible(visible) {
        this.#element.classList.toggle('hidden', !visible || this.#empty);
    }
}

class PanelControls {
    #viewer;

    constructor(element, viewer, meta) {
        const complex = meta.complex;

        element.classList.remove('hidden');
        this.#viewer = viewer;

        const chainLegend = element.querySelector('.sp-chain-legend');

        // hacky check if complex is user-uploaded
        const isUserUpload = /^\d+$/.test(meta.id);

        const tag = isUserUpload ? 'a' : 'span';

        chainLegend.innerHTML = complex.members.map((member, i) => {
            const color = subunitColors[i % subunitColors.length];

            const attributes = isUserUpload ?
                ` href="gene?id=${member.id}" target="_blank"` : '';

            return `<${tag} class="sp-color-marker sp-chain-marker"`
                + ` style="background-color: ${color}" ${attributes}>`
                + member.title
                + `</${tag}>`
        }).join('-');


        const metrics = element.querySelector('.sp-metrics');

        if (
            complex.plddt === null
            && complex.ptm === null
            && complex.iptm === null
        ) {
            metrics.classList.add('hidden');
        } else {
            if (complex.plddt === null) {
                metrics.querySelector(
                    '.sp-metrics-plddt'
                ).classList.add('hidden');
            } else {
                metrics.querySelector(
                    '[data-type="plddt"]'
                ).textContent = complex.plddt;
            }

            if (complex.ptm === null) {
                metrics.querySelector(
                    '.sp-metrics-ptm'
                ).classList.add('hidden');
            } else {
                metrics.querySelector(
                    '[data-type="ptm"]'
                ).textContent = complex.ptm;
            }

            if (complex.iptm === null) {
                metrics.querySelector(
                    '.sp-metrics-iptm'
                ).classList.add('hidden');
            } else {
                metrics.querySelector(
                    '[data-type="iptm"]'
                ).textContent = complex.iptm;
            }
        }

        metrics.querySelector('[data-type="iptm"]').textContent = complex.iptm;

        // used to make radio input group names unique
        const id = createRandomId();

        const structureSchemeOptions = element.querySelector(
            '.sp-color-options[data-scheme="structure"]'
        );

        for (const input of structureSchemeOptions.querySelectorAll('input')) {
            input.name = `structure-scheme-${id}`;
        }

        const crosslinkSchemeOptions = element.querySelector(
            '.sp-color-options[data-scheme="crosslinks"]'
        );

        if (complex.crosslinksUrl || complex.crosslinksFile) {
            crosslinkSchemeOptions.classList.remove('hidden');
        }

        for (const input of crosslinkSchemeOptions.querySelectorAll('input')) {
            input.name = `crosslink-scheme-${id}`;
        }

        const confidenceLegend = structureSchemeOptions.querySelector(
            '.sp-confidence-legend'
        );

        structureSchemeOptions.addEventListener('change', event => {
            confidenceLegend.classList.toggle(
                'hidden', event.target.value !== 'confidence'
            );

            this.#viewer.setStyle(defaultStyles.get(event.target.value));
        });

        crosslinkSchemeOptions.addEventListener('change', event => {
            this.#viewer.setCrosslinkScheme(event.target.value);
        });

        const downloadLink = element.querySelector('.sp-crosslinks-download');
        downloadLink.href = `examples/${complex.handle}.zip`;

        if (meta.description !== "User uploaded structure.") {
            const restraintNote = element.querySelector('.sp-restraint-note');

            if (restraintNote) {
                restraintNote.classList.remove('hidden');
            }
        }
    }
}