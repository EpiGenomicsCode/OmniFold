import {fetchDSV} from "./utils.js";

class StructureMeta {
    id;
    type;
    source;
    description;
    url;
    complex;
    crosslinks;

    constructor(id, type, source, description) {
        this.id = id;
        this.type = type;
        this.source = source;
        this.description = description;

        this.url = null;
        this.complex = null;
        this.crosslinks = [];
    }

    /**
     * Initializes structure by loading corresponding URL of the actual
     * structure file as well as complex data in the case of predicted
     * complexes.
     *
     * If successful, the returned promise resolves to the initialized
     * structure, otherwise rejects with null.
     *
     * @returns {Promise<StructureMeta|null>}
     */
    async init() {
        const load = new Map([
            ['pdb', this.loadFromPDB],
            ['alphafold', this.loadFromAlphaFold],
            ['predicted-complex', this.loadPredictedComplex]
        ]).get(this.source).bind(this);

        return load().then(() =>
            this.complex?.crosslinksFile ?
                this.#loadCrosslinks() : this
        );
    }

    /**
     * Finds most suitable structure file by trying to fetch PDB URLs
     * (starting with biological assemblies) and storing
     * the first successful hit.
     *
     * If successful, the returned promise resolves to the initialized
     * structure, otherwise rejects with null.
     *
     * @returns {Promise<StructureMeta|null>}
     */
    async loadFromPDB() {
        const urlPrefix = `https://files.rcsb.org/download/${this.id}`;
        const extensions = ['-assembly1.cif', '.pdb1.gz', '.cif.gz', '.pdb.gz'];
        const urls = extensions.map(extension => urlPrefix + extension);

        return new Promise(async (resolve, reject) => {
            for (const url of urls) {
                const ok = await fetch(url, {method: 'HEAD'}).then(
                    response => response.ok
                ).catch(() => false); // ignore unsuccessful attempt

                if (ok) {
                    this.url = url;
                    resolve(this);
                    return;
                }
            }

            reject(`Couldn't load '${this.id}' from PDB!`);
        });
    }

    /**
     * Loads structure file URL from AlphaFold Protein Structure
     * database.
     *
     * If successful, the returned promise resolves to the initialized
     * structure, otherwise rejects with null.
     *
     * @returns {Promise<StructureMeta|null>}
     */
    async loadFromAlphaFold() {
        return fetch(`https://alphafold.ebi.ac.uk/api/prediction/${this.id}`)
            .then(response => response.json())
            .then(data => {
                this.url = data[0].pdbUrl;
                return this;
            }).catch(error => {
                console.error(
                    `Couldn't load '${this.id}' from AlphaFold Protein`
                    + ` Structure Database!`
                );
                throw error;
            });
    }

    /**
     * Loads structure file URL and complex data from server.
     *
     * If successful, the returned promise resolves to the initialized
     * structure, otherwise rejects with null.
     *
     * @returns {Promise<StructureMeta|null>}
     */
    async loadPredictedComplex() {
        const directory = `examples/${this.id}/`;

        return fetch(
            directory + 'meta.json',
            {'headers': {'Accept': 'application/json'}}
        ).then(response => response.json()).then(complex => {
            this.url = directory + 'structure.pdb';
            this.complex = complex;

            // workaround to give unique IDs and names
            const uniprotCounts = new Map();

            for (const member of this.complex.members) {
                const count = (uniprotCounts.get(member.uniprot) ?? 0) + 1;
                uniprotCounts.set(member.uniprot, count);

                if (count > 1) {
                    member.uniprot = `${member.uniprot}#${count}`;
                    member.title = `${member.title}#${count}`;
                }
            }

            // mapping from position to chain letter (A, B, C...)
            this.complex.chains = getChainMapping(complex);

            return this;
        }).catch(error => {
            console.error(`Couldn't load '${this.id}' from server!`);
            throw error;
        });
    }

    async #loadCrosslinks() {
        return fetchDSV(
            this.complex.crosslinksFile,
            null, ','
        ).then(crosslinks => {
            this.crosslinks = processCrosslinks(
                crosslinks, this.complex.members, this.complex.chains
            );

            return this;
        });
    }
}

/**
 * Tries to initialize structure from passed structure meta
 * (containing id, type, source and description).
 *
 * If successful, the returned promise resolves to the initialized
 * structure, otherwise rejects with null.
 *
 * @param meta
 * @returns {Promise<StructureMeta|null>}
 */
export async function createStructureMeta(meta) {
    const structure = new StructureMeta(
        meta.id, meta.type, meta.source, meta.description
    );

    return structure.init();
}

export function getChainMapping(complex) {
    return complex.members.reduce((obj, k, i) => (
        {...obj, [k.uniprot]: String.fromCharCode(65 + i % 26)}
    ), {});
}


function checkCrosslinkTableHeaders(crosslinks) {
    for (const requiredHeader of [
        'Protein1', 'SeqPos1', 'Protein2', 'SeqPos2'
    ]) {
        if (!crosslinks[0].has(requiredHeader)) {
            throw {
                name: 'MissingHeader',
                message: `Missing required header '${requiredHeader}' in`
                    + ` crosslinks table!`,
                missingHeader: requiredHeader,
            }
        }
    }
}

function checkCrosslinkTableChains(crosslinks, chains) {
    const chainsFromCrosslinks = new Set(
        crosslinks.map(crosslink => crosslink.get('Protein1')),
        crosslinks.map(crosslink => crosslink.get('Protein2')),
    );

    const chainIdentifiers = Object.keys(chains);

    for (const chainFromCrosslink of chainsFromCrosslinks) {
        if (!chainIdentifiers.includes(chainFromCrosslink)) {
            throw {
                name: 'InvalidChainId',
                message: `Invalid chain ID '${chainFromCrosslink}' encountered`
                    + ` while processing crosslinks (valid chains:`
                    + ` '${chainIdentifiers.join(',')}')!`,
                invalidChain: chainFromCrosslink,
                validChains: chainIdentifiers
            }
        }
    }
}

function checkCrosslinkTablePositions(crosslinks, subunits) {
    const subunitLengths = Object.fromEntries(subunits.map(
        subunit => [subunit.id, subunit.length]
    ));

    for (const [rowIndex, crosslink] of crosslinks.entries()) {
        for (const i of [1, 2]) {
            const posKey = `SeqPos${i}`;
            const posString = crosslink.get(posKey);

            if (!/^\d+$/.test(posString)) {
                throw {
                    name: 'MalformedSequencePosition',
                    message: `Malformed sequence position '${posString}' for`
                        + ` '${posKey}' in crosslink '${rowIndex}'!`
                        + ` Must be a positive integer >= 1!`,
                    key: posKey,
                    value: posString
                }
            }

            const pos = parseInt(posString);

            if (pos < 1) {
                throw {
                    name: 'InvalidSequencePosition',
                    message: `Invalid sequence position '${pos}' for`
                        + ` '${posKey}' in crosslink '${rowIndex}'!`
                        + ` Must be a positive integer >= 1!`,
                    key: posKey,
                    value: pos
                }
            }

            const proteinKey = `Protein${i}`;
            const chain = crosslink.get(proteinKey);

            if (pos > subunitLengths[chain]) {
                throw {
                    name: 'InvalidSequencePosition',
                    message: `Invalid sequence position '${pos}' for`
                        + ` '${posKey}' in crosslink '${rowIndex}'!`
                        + ` The value can't be higher than the length of chain`
                        + ` '${chain}' (${subunitLengths[chain]})!`,
                    chain: chain,
                    length: subunitLengths[chain],
                    key: posKey,
                    value: pos
                }
            }
        }
    }
}

function checkCrosslinkTableRestraints(crosslinks) {
    if (!crosslinks[0].has('RestraintSatisfied')) {
        return;
    }

    for (const [rowIndex, crosslink] of crosslinks.entries()) {
        const value = crosslink.get('RestraintSatisfied').toLowerCase();

        if (value !== 'true' && value !== 'false') {
            throw {
                name: 'InvalidRestraint',
                message: `Invalid 'RestraintSatisfied' value '${value}'`
                    + ` in crosslink '${rowIndex}', must be either 'true' or`
                    + ` 'false' (case-insensitive)!`,
                value: value,
                rowIndex: rowIndex
            }
        }
    }
}

export function processCrosslinks(crosslinks, subunits, chains) {
    if (crosslinks.length === 0) return [];

    checkCrosslinkTableHeaders(crosslinks);
    checkCrosslinkTableChains(crosslinks, chains);
    checkCrosslinkTablePositions(crosslinks, subunits);
    checkCrosslinkTableRestraints(crosslinks);

    for (const [i, crosslink] of crosslinks.entries()) {
        crosslink.set('id', i);
        crosslink.set('atoms', getAtomsFromCrosslink(crosslink, chains));
    }

    return crosslinks;
}

export function getAtomsFromCrosslink(link, chains, atom = 'CA') {
    const atom1 = link.get('Atom1') ?? atom;
    const atom2 = link.get('Atom2') ?? atom;

    return [
        `${link.get('SeqPos1')}:${chains[link.get('Protein1')]}.${atom1}`,
        `${link.get('SeqPos2')}:${chains[link.get('Protein2')]}.${atom2}`,
    ];
}
