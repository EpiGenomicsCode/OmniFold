export class StructureStyle {
    #styles;

    constructor(name, styleParams) {
        this.name = name;
        this.#styles = styleParams.map(([type, params, modifier]) =>
            new RepresentationStyle(type, {...params, 'name': name}, modifier)
        );
    }

    /**
     * Creates representations for the passed component using all stored
     * styles.
     *
     * @param {NGL.Component} component
     * @param stage
     * @param meta
     * @param chooseColoring
     * @param visible
     */
    render(component, stage, meta, chooseColoring = true, visible = true) {
        return Promise.all(this.#styles.map(style =>
            style.render(component, stage, meta, chooseColoring, visible)
        ));
    }
}

export class RepresentationStyle {
    constructor(type, params = {}, modifier = null) {
        this.type = type;
        this.params = params;
        this.modifier = modifier;
    }

    /**
     * Creates representation for the component according to the styles
     * and its optional modifier function. Returns promise which
     * resolves with completion of the rendering (surface rendering
     * might take a while).
     *
     * @param {NGL.Component} component
     * @param stage
     * @param meta
     * @param chooseColoring
     * @param visible
     */
    render(component, stage, meta, chooseColoring = true, visible = true) {
        return new Promise(resolve => {
            const colorScheme = chooseColoring ?
                this.#chooseColoring(component, meta) : {};

            let dynamicColorScheme = {};

            if (Object.hasOwn(this.params, 'colorScheme')) {
                if (Array.isArray(this.params.colorScheme)) {
                    const colors = this.params.colorScheme;

                    const scheme = meta.complex.members
                        .map(subunit => subunit.uniprot)
                        .map(uniprot => meta.complex.chains[uniprot])
                        .map((chain, i) => [
                            colors[i % colors.length],
                            ':' + chain
                        ]);

                    dynamicColorScheme = {
                        'colorScheme':
                            NGL.ColormakerRegistry.addSelectionScheme(
                                scheme
                            )
                    };
                } else if (this.params.colorScheme === 'bfactor') {
                    // use pLDDT schema from AlphaFold DB
                    const scheme = NGL.ColormakerRegistry.addScheme(function (e) {
                        this.atomColor = function (atom) {
                            if (atom.bfactor > 90) {
                                return 0x0053D6;  // very high, dark blue
                            } else if (atom.bfactor > 70) {
                                return 0x64CBF3;  // confident, light blue
                            } else if (atom.bfactor > 50) {
                                return 0xFFDB12;  // low, yellow
                            } else {
                                return 0xFF7D45;  // very low, orange
                            }
                        };
                    });

                    dynamicColorScheme = {'color': scheme};
                }
            }

            const representation = component.addRepresentation(this.type, {
                name: this.name, ...colorScheme, ...this.params,
                ...dynamicColorScheme
            });

            representation.setVisibility(visible);

            if (this.modifier !== null) {
                this.modifier(component);
            }

            if (this.type !== 'surface') {
                resolve(representation);
                return;
            }

            const finishSurfaceRender = (delta, _) => {
                if (delta === -1) {
                    resolve(representation);
                }
            }

            stage.tasks.signals.countChanged.add(finishSurfaceRender);
        })
    }

    #chooseColoring(component, meta) {
        if (component.structure.modelStore.count > 1) {
            return {colorScheme: 'modelindex'};
        } else if (component.structure.chainStore.count > 1) {
            return {colorScheme: 'chainindex'};
        } else {
            return meta.source === 'alphafold' ?
                {colorScheme: 'bfactor'}
                : {colorScheme: 'atomindex', colorScale: 'RdYlBu'};
        }
    }
}

// modified from Okabe_Ito
export const subunitColors = [
    '#991999', // PyMol deeppurple (0.6, 0.1, 0.6)
    '#00BFBF', // PyMol teal (0, 0.75, 0.75)
    '#e9967a', // salmon
    '#009e73',
    '#f0e442',
    '#0072b2',
    '#d55e00',
    '#cc79a7'
];

export const ternaryColoring = [
    ['#991999', ':A'],  // PyMol deeppurple (0.6, 0.1, 0.6)
    ['#00BFBF', ':B'],  // PyMol teal (0, 0.75, 0.75)
    ['#e9967a', ':C'],  // salmon
];

export const trimerScheme = NGL.ColormakerRegistry.addSelectionScheme(
    ternaryColoring, 'Trimer'
);

export const surfaceDefaultStyle = {
    opacity: 0.2,
    opaqueBack: false,
    sele: 'protein'
};

export const crosslinkStyleSelected = {
    opacity: 0.9,
    lineOpacity: 0.9,
    labelBackgroundOpacity: 0.9
};

export const crosslinkStyleDeselected = {
    opacity: 0.2,
    lineOpacity: 0.2,
    labelBackgroundOpacity: 0.1
};

export const crosslinkStyleCommon = {
    ...crosslinkStyleSelected,
    labelSize: 2,
    labelBackground: true,
    labelColor: 'black',
    linewidth: 8,
    labelUnit: 'angstrom',
    useCylinder: true,
    radius: 1
};


export const defaultStyles = new Map([
    ['index', [['cartoon', {}]]],
    ['hydrophobicity', [
        ['ball+stick', {colorScheme: 'hydrophobicity', sele: 'protein'}],
        ['surface', {
            ...surfaceDefaultStyle,
            colorScheme: 'hydrophobicity',
            opacity: 0.7
        }]
    ]],
    ['electrostatic', [
        ['ball+stick', {colorScheme: 'electrostatic', sele: 'protein'}],
        ['surface', {
            ...surfaceDefaultStyle,
            colorScheme: 'electrostatic',
            opacity: 0.7
        }]
    ]],
    ['subunit', [
        ['cartoon', {'colorScheme': subunitColors}],
        ['surface', {
            'colorScheme': subunitColors,
            opacity: 0.1,
            opaqueBack: false
        }],
        ['base', {sele: 'nucleic', 'colorScheme': subunitColors}],
        ['ball+stick', {sele: 'ligand', 'colorScheme': subunitColors}]
    ]],
    ['confidence', [
        ['cartoon', {'colorScheme': 'bfactor'}],
        ['base', {sele: 'nucleic', 'colorScheme': 'bfactor'}],
        ['ball+stick', {sele: 'ligand', 'colorScheme': 'bfactor'}]
    ]],
    ['selection', [
        ['cartoon', {}],
        ['base', {sele: 'nucleic'}],
        ['ball+stick', {sele: 'ligand'}]
    ]]
].map(([name, styles]) => [name, new StructureStyle(name, styles)]));
