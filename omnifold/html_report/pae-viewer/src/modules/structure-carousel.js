import {createStructureMeta} from "./structure-meta.js";
import {StructureViewer} from "./structure-viewer.js";
import {createRandomId, positiveModulo} from "./utils.js";
import {defaultStyles} from "./structure-style.js";
import {StructurePanel} from "./structure-panel.js";

export class StructureCarousel {
    #root;
    #template;
    #element;
    #viewerTemplate;
    #panelTemplate;
    #panelRoot;
    #container;
    #carousel;
    #controls;
    #inner;
    #stylesDropDown;
    #panel;

    #gene;
    #entries;
    #currentIndex;

    constructor(root, template, viewerTemplate, panelTemplate, panelRoot = null, gene = null) {
        this.#root = root;
        this.#template = template;
        this.#element = this.#template.content.querySelector(
            '.structure-carousel'
        ).cloneNode(true);

        this.#viewerTemplate = viewerTemplate;
        this.#panelTemplate = panelTemplate;

        this.#container = this.#element.querySelector(
            '.sc-container:not(.sc-placeholder)'
        );

        this.#setupFullscreenControls();

        this.#carousel = this.#container.querySelector('.carousel');
        this.#carousel.id = createRandomId('carousel');

        this.#controls = this.#carousel.querySelectorAll(
            '[data-bs-target], [data-bs-slide], [data-bs-slide-to]'
        );

        for (const control of this.#controls) {
            control.setAttribute('data-bs-target', `#${this.#carousel.id}`);
        }
        new bootstrap.Carousel(this.#carousel);

        this.#stylesDropDown
            = this.#carousel.querySelector('.sc-options');

        this.#inner = this.#carousel.querySelector('.sc-inner');
        this.#panel = this.#container.querySelector('.sc-panel');

        this.#panelRoot = panelRoot ? panelRoot : this.#panel;

        this.#root.appendChild(this.#element);

        this.#entries = [];

        window.addEventListener('resize', this.resize.bind(this), false);
        this.resize();

        // create 'help' tooltip in bottom right corner
        const tooltip = new bootstrap.Tooltip(
            this.#carousel.querySelector('[data-bs-toggle="tooltip"]')
        );

        // if user clicks somewhere aside from the tooltip, hide it
        document.addEventListener('click', event => {
            if (!(this.#element.contains(event.target)
                || (tooltip.tip !== null
                    && tooltip.tip.contains(event.target)))) {
                tooltip.hide();
            }
        });

        this.#carousel.addEventListener('slide.bs.carousel', event => {
             this.setCurrentIndex(event.to);
        });

        this.#stylesDropDown.addEventListener('change', event => {
            this.#updateStyle();
        });

        this.#currentIndex = 0;

        if (gene !== null) {
            this.fetchStructures(gene);
        }
    }

    getElement() {
        return this.#element;
    }

    #setupFullscreenControls() {
        const toggleFullscreenBox = this.#container.querySelector(
            '.sc-fullscreen-toggle'
        );

        const enterFullscreenButton = toggleFullscreenBox.querySelector(
            '.sc-enter-fullscreen'
        );

        const exitFullscreenButton = toggleFullscreenBox.querySelector(
            '.sc-exit-fullscreen'
        );

        const placeholder = this.#element.querySelector('.sc-placeholder');

        const setFullscreen = fullscreen => {
            this.#container.classList.toggle('sc-pop-up', fullscreen)
            placeholder.classList.toggle('hidden', !fullscreen);
            enterFullscreenButton.classList.toggle('hidden', fullscreen);
            exitFullscreenButton.classList.toggle('hidden', !fullscreen);
            this.resize();
        }

        const toggleFullscreen = () => setFullscreen(
            !this.#container.classList.contains('sc-pop-up')
        );

        toggleFullscreenBox.addEventListener('click', toggleFullscreen);

        this.#container.addEventListener('click', event => {
            if (
                event.target === event.currentTarget
                && event.currentTarget.classList.contains('sc-pop-up')
            ) {
                setFullscreen(false);
            }
        });
    }

    #makeSquare() {
        if (this.#container.classList.contains('pop-up')) {
            this.#carousel.style.removeProperty('height');
        } else {
            this.#carousel.style.height = Math.round(
                this.#carousel.getBoundingClientRect().width
            ) + "px";
        }
    }

    async fetchStructures(id) {
        this.#gene = id;

        return fetch(
            `gene/structures?id=${this.#gene}`,
            {'headers': {'Accept': 'application/json'}}
        ).then(response => response.json())
         .then(this.addStructures.bind(this))
         .catch(() => console.error(
            `Couldn't fetch structures for gene '${id}'!`
          )
        );
    }

    async addStructures(metas) {
        for (const [i, meta] of metas.entries()) {
            await this.addStructure(meta, i === 0);
        }
    }

    async addStructure(meta, showPanel) {
        return createStructureMeta(meta)
            .then(structure => this.addStructureMeta(structure, showPanel))
            .catch(message => {
                console.error(message);
                this.#entries.push(null);
            });
    }

    addStructureMeta(structure, showPanel) {
        const item = document.createElement('div');
        item.classList.add('carousel-item');
        this.#inner.appendChild(item);

        const rect = this.#inner.getBoundingClientRect();

        const viewer = new StructureViewer(
            item, this.#viewerTemplate, structure,
            [rect.width, rect.height],
            structure.source === 'predicted-complex' ?
                defaultStyles.get('subunit') : this.#getStyle()
        );

        const entry = {
            'structure': structure,
            'item': item,
            'viewer': viewer,
            'panel': new StructurePanel(
                this.#panelRoot, this.#panelTemplate, viewer, showPanel
            )
        };

        this.#entries.push(entry);

        if (this.#entries.length === 1) {
            this.#entries[0].item.classList.add('active');
        } else if (this.#entries.length > 1) {
            for (const control of this.#controls) {
                control.classList.remove('hidden');
            }
        }

        return entry;
    }

    getGene() {
        return this.#gene;
    }

    previous() {
        this.setCurrentIndex(positiveModulo(
            ++this.#currentIndex, this.#entries.length
        ));
    }

    next() {
        this.setCurrentIndex(positiveModulo(
            ++this.#currentIndex, this.#entries.length
        ));
    }

    setCurrentIndex(index) {
        this.#currentIndex = index;
        this.#updateStyle();
        this.#updatePanel();
    }

    #updateStyle() {
        this.#stylesDropDown.disabled = true;

        if (this.getCurrent().structure.source === 'predicted-complex') {
            return;
        }

        this.getCurrent().viewer.setStyle(
            this.#getStyle()
        ).then(() => {
            if (this.getCurrent().structure.source !== 'predicted-complex') {
                this.#stylesDropDown.disabled = false;
            }
        }).catch(() => {
            if (this.getCurrent().structure.source !== 'predicted-complex') {
                this.#stylesDropDown.disabled = false;
            }
        });
    }

    #updatePanel() {
        for (const [i, structure] of this.#entries.entries()) {
            structure.panel.setVisible(i === this.#currentIndex);
        }
    }

    #getStyle() {
        return defaultStyles.get(this.#stylesDropDown.value);
    }

    getCurrent() {
        return this.#entries[this.#currentIndex];
    }

    resize() {
        if (!CSS.supports('aspect-ratio', '1 / 1')) {
            this.#makeSquare();
        }

        const rect = this.#inner.getBoundingClientRect();

        for (const structure of this.#entries) {
            structure.viewer.setSize(rect.width, rect.height);
        }
    }
}
