import {getFileExtension, readAsText, readDSV, splitLines} from "./utils.js";
import {processCrosslinks} from "./structure-meta.js";

export function readCrosslinks(complex) {
    if (complex.crosslinksUrl === null) {
        return([complex, []]);
    }

    const path = complex.crosslinksUrl;
    const extension = getFileExtension(path).toLowerCase();

    if (!['csv', 'pb'].includes(extension)) {
        throw {
            type: 'crosslinks-input-error',
            message: `Invalid extension '.${extension}' of crosslinks file!`
                + ` Only .csv and .pb are supported!`
        };
    }

    return readAsText(path).catch(() => {
        throw {
            type: 'crosslinks-input-error',
            message: "Couldn't read crosslinks file!"
        };
    }).then(text => {
        if (extension === 'csv') {
            return readCrosslinksFromCsv(text, complex);
        } else { // .pb
            return readCrosslinksFromPseudobonds(text, complex)
        }
    });
}

export function readCrosslinksFromCsv(text, complex) {
    let crosslinks = [];

    try {
        crosslinks = readDSV(text, null, ',');
    } catch (error) {
        if (error.name === 'MalformattedDsvRow') {
            throw {
                type: 'crosslinks-input-error',
                message:
                    `Crosslinks table is possibly malformatted,`
                    + ` row ${error.index} has ${error.row.length}`
                    + ` fields, but ${error.headers.length} headers`
                    + ` are defined!`
            };
        } else {
            throw {
                type: 'crosslinks-input-error',
                message: "Unknown error occurred while parsing" +
                    " the crosslinks table!"
            };
        }
    }

    try {
        crosslinks = processCrosslinks(
            crosslinks, complex.members, complex.chains
        );

        return [complex, crosslinks];
    } catch (error) {
        handleProcessingError(error);
    }
}

export function readCrosslinksFromPseudobonds(text, complex) {
    let crosslinks = [];

    for (const [i, line] of splitLines(text).entries()) {
        const trimmed = line.trim();

        if (!trimmed || trimmed.startsWith(';')) {
            continue;
        }

        const matches = [...trimmed.matchAll(/\/(\S+):(\d+)@(\S+)/g)];

        if (matches.length !== 2) {
            throw {
                type: 'crosslinks-input-error',
                message: `Pseudobond files must contain two atom specifiers!`
                    + ` However, ${atoms.length} was/were found for`
                    + ` pseudobond ${i} ('${trimmed}')!`,
            }
        }

        const atoms = matches.map(match => ({
            chain: match[1].toUpperCase(),
            residue: match[2],
            atom: match[3].toUpperCase()
        }));

        crosslinks.push(new Map([
            ['Protein1', atoms[0].chain],
            ['SeqPos1', atoms[0].residue],
            ['Atom1', atoms[0].atom],
            ['Protein2', atoms[1].chain],
            ['SeqPos2', atoms[1].residue],
            ['Atom2', atoms[1].atom],
        ]));
    }

    try {
        crosslinks = processCrosslinks(
            crosslinks, complex.members, complex.chains
        );

        return [complex, crosslinks];
    } catch (error) {
        if (error.name === 'InvalidChainId') {
            throw {
                type: 'crosslinks-input-error',
                message:
                    `An invalid chain '${error.invalidChain}' was referenced in`
                    + ` a pseudobond. Please note that user-defined chain`
                    + ` labels are not supported for .pb files.`
            };
        }
        try {
            handleProcessingError(error);
        } catch (error) {
            error.message += " Please note that validation feedback is"
                + " optimized for crosslink tables in CSV format, so error"
                + " messages might not be accurate for pseudobond files. Please"
                + " refer to the 'Help' section for limitations when using"
                + " pseudobond files.";

            throw error;
        }
    }
}

function handleProcessingError(error) {
    switch (error.name) {
        case 'InvalidChainId':
            throw {
                type: 'crosslinks-input-error',
                message:
                    `An invalid chain '${error.invalidChain}' was referenced in`
                    + ` the crosslink table (valid chains: `
                    + `'${error.validChains}').`
                    + ` You can use the optional 'Chain labels' input to change`
                    + ` the complex chain names.`
            };
        case 'MissingHeader':
            throw {
                type: 'crosslinks-input-error',
                message:
                    `The column '${error.missingHeader}' is missing from the`
                    + ` crosslinks table! See the 'Help' section for the table`
                    + ` structure.`
            };
        case 'MalformedSequencePosition':
        case 'InvalidSequencePosition':
        case 'InvalidRestraint':
            throw {
                type: 'crosslinks-input-error',
                message: error.message,
            };
        default:
            throw {
                type: 'crosslinks-input-error',
                message: "Unknown error occurred while processing the"
                    + " crosslinks table!"
            };
    }
}
