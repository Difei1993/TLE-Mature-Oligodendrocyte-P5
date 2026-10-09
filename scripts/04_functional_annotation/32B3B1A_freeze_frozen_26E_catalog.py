###############################################################################
# 32B-3B1A
# Locate, audit and freeze the exact 26E senescence/competing gene-set catalog.
#
# This stage performs no enrichment and no disease test.
###############################################################################

from __future__ import annotations

import hashlib
import io
import json
import re
import shutil
import zipfile
from datetime import datetime
from pathlib import Path

import pandas as pd


###############################################################################
# 0. Paths
###############################################################################

PROJECT = Path(
    __import__("os").environ["P5_WORK_ROOT"]
)

STAMP = datetime.now().strftime("%Y%m%d_%H%M%S")

OUTPUT_ROOT = (
    PROJECT
    / "99_frozen_releases"
    / f"32B3B1A_FROZEN_26E_CATALOG_AUDIT_{STAMP}"
)

OUTPUT_ROOT.mkdir(
    parents=True,
    exist_ok=False,
)

AUDIT_ROOT = OUTPUT_ROOT / "audit"
FREEZE_ROOT = OUTPUT_ROOT / "freeze"

AUDIT_ROOT.mkdir()
FREEZE_ROOT.mkdir()


###############################################################################
# 1. Exact expected 26E program universe
###############################################################################

PRIMARY_PROGRAMS = [
    "DDR_genome_instability",
    "Stable_cell_cycle_exit",
    "SASP_ECM_remodelling",
    "Lysosome_autophagy_proteostasis",
    "Anti_apoptosis_survival",
    "Mitochondrial_oxidative_stress",
    "Chromatin_remodelling",
    "Lamina_nuclear_integrity",
]

EXTERNAL_SIGNATURES = [
    "SAUL_SEN_MAYO",
    "FRIDMAN_SENESCENCE_UP",
    "REACTOME_CELLULAR_SENESCENCE",
    "REACTOME_SENESCENCE_ASSOCIATED_SECRETORY_PHENOTYPE_SASP",
]

DECONFOUNDED_PROGRAMS = [
    "DDR_genome_instability_deconfounded",
    "Stable_cell_cycle_exit_deconfounded",
    "SASP_ECM_remodelling_deconfounded",
    "Lysosome_autophagy_proteostasis_deconfounded",
    "Anti_apoptosis_survival_deconfounded",
    "Mitochondrial_oxidative_stress_deconfounded",
    "Chromatin_remodelling_deconfounded",
]

COMPETING_PROGRAMS = [
    "Apoptosis_cell_death",
    "IFN_response",
    "Immediate_early_response",
    "JAK_STAT_activation",
    "NFkB_inflammation_chemokine",
    "Phagolysosome_lipid_processing",
    "Proliferation",
]

EXPECTED_PROGRAMS = (
    PRIMARY_PROGRAMS
    + EXTERNAL_SIGNATURES
    + DECONFOUNDED_PROGRAMS
    + COMPETING_PROGRAMS
)

if len(EXPECTED_PROGRAMS) != 26:
    raise RuntimeError(
        f"Expected 26 programs, obtained {len(EXPECTED_PROGRAMS)}."
    )

if len(set(EXPECTED_PROGRAMS)) != 26:
    raise RuntimeError("Duplicate program names in expected 26E universe.")


###############################################################################
# 2. Frozen program metadata
###############################################################################

PROGRAM_METADATA: dict[str, dict[str, object]] = {}

for program in PRIMARY_PROGRAMS:
    PROGRAM_METADATA[program] = {
        "analysis_class": "Senescence-associated",
        "evidence_layer": (
            "Inverse supplemental dimension"
            if program == "Lamina_nuclear_integrity"
            else "Primary senescence dimension"
        ),
        "expected_direction": (
            "DOWN"
            if program == "Lamina_nuclear_integrity"
            else "UP"
        ),
        "direction_multiplier": (
            -1
            if program == "Lamina_nuclear_integrity"
            else 1
        ),
    }

for program in EXTERNAL_SIGNATURES:
    PROGRAM_METADATA[program] = {
        "analysis_class": "Senescence-associated",
        "evidence_layer": "External published signature",
        "expected_direction": "UP",
        "direction_multiplier": 1,
    }

for program in DECONFOUNDED_PROGRAMS:
    PROGRAM_METADATA[program] = {
        "analysis_class": "Senescence-associated",
        "evidence_layer": "Deconfounded sensitivity dimension",
        "expected_direction": "UP",
        "direction_multiplier": 1,
    }

for program in COMPETING_PROGRAMS:
    PROGRAM_METADATA[program] = {
        "analysis_class": "Competing",
        "evidence_layer": "Competing biological program",
        "expected_direction": "RAW",
        "direction_multiplier": 1,
    }


###############################################################################
# 3. Utilities
###############################################################################

def now_string() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def sha256_bytes(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for block in iter(
            lambda: handle.read(1024 * 1024),
            b"",
        ):
            digest.update(block)

    return digest.hexdigest()


def normalize_key(value: object) -> str:
    text = str(value).strip().lower()

    text = re.sub(
        r"[^a-z0-9]+",
        "_",
        text,
    )

    return text.strip("_")


EXPECTED_KEY_TO_NAME = {
    normalize_key(program): program
    for program in EXPECTED_PROGRAMS
}


def canonical_program_name(
    value: object,
) -> str | None:

    return EXPECTED_KEY_TO_NAME.get(
        normalize_key(value)
    )


def normalize_gene(value: object) -> str | None:
    gene = str(value).strip()

    gene = gene.strip(
        "\"'[](){}"
    )

    if not gene:
        return None

    if gene.lower() in {
        "nan",
        "na",
        "none",
        "null",
        "gene",
        "genes",
    }:
        return None

    gene = gene.upper()

    if not re.fullmatch(
        r"[A-Z0-9][A-Z0-9_.\-]*",
        gene,
    ):
        return None

    return gene


def split_gene_value(
    value: object,
) -> list[str]:

    if value is None:
        return []

    raw = str(value).strip()

    if not raw:
        return []

    raw = raw.strip(
        "[](){}"
    )

    parts = re.split(
        r"[;,|\s]+",
        raw,
    )

    genes = []

    for part in parts:
        gene = normalize_gene(part)

        if gene is not None:
            genes.append(gene)

    return genes


def canonical_table_hash(
    table: pd.DataFrame,
) -> str:

    canonical = (
        table[
            [
                "program",
                "gene",
            ]
        ]
        .drop_duplicates()
        .sort_values(
            [
                "program",
                "gene",
            ]
        )
        .reset_index(drop=True)
    )

    content = canonical.to_csv(
        index=False,
        lineterminator="\n",
    ).encode("utf-8")

    return sha256_bytes(content)


###############################################################################
# 4. Candidate parsers
###############################################################################

PROGRAM_COLUMNS = [
    "program",
    "gene_set",
    "geneset",
    "gene_set_name",
    "set_name",
    "signature",
    "module",
    "pathway",
    "name",
]

GENE_COLUMNS = [
    "gene",
    "gene_symbol",
    "symbol",
    "genes",
    "members",
    "gene_list",
    "gene_symbols",
    "geneset_genes",
]


def clean_column_name(value: object) -> str:
    return normalize_key(value)


def parse_tabular_bytes(
    content: bytes,
    suffix: str,
) -> pd.DataFrame:

    separator = (
        "\t"
        if suffix.lower() in {
            ".tsv",
            ".txt",
        }
        else ","
    )

    table = pd.read_csv(
        io.BytesIO(content),
        sep=separator,
        dtype=str,
        keep_default_na=False,
        low_memory=False,
    )

    original_columns = list(
        table.columns
    )

    clean_columns = [
        clean_column_name(column)
        for column in original_columns
    ]

    table.columns = clean_columns

    rows: list[dict[str, str]] = []

    # Long-format table.
    program_column = next(
        (
            column
            for column in PROGRAM_COLUMNS
            if column in table.columns
        ),
        None,
    )

    gene_column = next(
        (
            column
            for column in GENE_COLUMNS
            if column in table.columns
        ),
        None,
    )

    if (
        program_column is not None
        and gene_column is not None
    ):
        for row in table.itertuples(
            index=False
        ):
            row_dictionary = row._asdict()

            program = canonical_program_name(
                row_dictionary[
                    program_column
                ]
            )

            if program is None:
                continue

            for gene in split_gene_value(
                row_dictionary[
                    gene_column
                ]
            ):
                rows.append(
                    {
                        "program": program,
                        "gene": gene,
                    }
                )

    # Wide-format table: one program per column.
    for column in table.columns:
        program = canonical_program_name(
            column
        )

        if program is None:
            continue

        for value in table[column]:
            for gene in split_gene_value(
                value
            ):
                rows.append(
                    {
                        "program": program,
                        "gene": gene,
                    }
                )

    return (
        pd.DataFrame(
            rows,
            columns=[
                "program",
                "gene",
            ],
        )
        .drop_duplicates()
    )


def parse_gmt_bytes(
    content: bytes,
) -> pd.DataFrame:

    text = content.decode(
        "utf-8",
        errors="replace",
    )

    rows: list[dict[str, str]] = []

    for line in text.splitlines():
        fields = line.rstrip(
            "\n\r"
        ).split("\t")

        if len(fields) < 3:
            continue

        program = canonical_program_name(
            fields[0]
        )

        if program is None:
            continue

        for raw_gene in fields[2:]:
            gene = normalize_gene(
                raw_gene
            )

            if gene is not None:
                rows.append(
                    {
                        "program": program,
                        "gene": gene,
                    }
                )

    return (
        pd.DataFrame(
            rows,
            columns=[
                "program",
                "gene",
            ],
        )
        .drop_duplicates()
    )


def parse_candidate_bytes(
    content: bytes,
    suffix: str,
) -> pd.DataFrame:

    if suffix.lower() == ".gmt":
        return parse_gmt_bytes(
            content
        )

    return parse_tabular_bytes(
        content,
        suffix,
    )


###############################################################################
# 5. Locate candidate files
###############################################################################

SKIP_DIRECTORIES = {
    ".git",
    ".rproj.user",
    "__pycache__",
    "python_envs",
    "00_raw_data",
    "02_objects",
}

NAME_PATTERN = re.compile(
    r"(senesc|geneset|gene_set|gene-set|signature|catalog|program)",
    flags=re.IGNORECASE,
)

SUPPORTED_SUFFIXES = {
    ".csv",
    ".tsv",
    ".txt",
    ".gmt",
}

candidate_sources: list[
    dict[str, object]
] = []


for path in PROJECT.rglob("*"):

    if any(
        part.lower() in SKIP_DIRECTORIES
        for part in path.parts
    ):
        continue

    if not path.is_file():
        continue

    if path.stat().st_size <= 0:
        continue

    # Direct candidate files.
    if (
        path.suffix.lower()
        in SUPPORTED_SUFFIXES
        and NAME_PATTERN.search(
            path.name
        )
        and path.stat().st_size
        <= 50 * 1024 * 1024
    ):
        candidate_sources.append(
            {
                "source_type": "file",
                "source_path": str(path),
                "source_name": path.name,
                "suffix": path.suffix.lower(),
                "content": path.read_bytes(),
            }
        )

    # Search only exact catalog-like names inside ZIPs.
    if (
        path.suffix.lower() == ".zip"
        and path.stat().st_size
        <= 2 * 1024 * 1024 * 1024
    ):
        try:
            with zipfile.ZipFile(
                path
            ) as archive:

                for entry in archive.infolist():

                    if entry.is_dir():
                        continue

                    entry_name = Path(
                        entry.filename
                    ).name

                    entry_suffix = Path(
                        entry_name
                    ).suffix.lower()

                    if (
                        entry_suffix
                        not in SUPPORTED_SUFFIXES
                    ):
                        continue

                    if not NAME_PATTERN.search(
                        entry_name
                    ):
                        continue

                    if (
                        "SENESCENCE_GENESET_CATALOG"
                        not in entry_name.upper()
                        and "26E" not in entry_name.upper()
                    ):
                        continue

                    content = archive.read(
                        entry
                    )

                    candidate_sources.append(
                        {
                            "source_type":
                                "zip_entry",

                            "source_path":
                                (
                                    str(path)
                                    + "::"
                                    + entry.filename
                                ),

                            "source_name":
                                entry_name,

                            "suffix":
                                entry_suffix,

                            "content":
                                content,
                        }
                    )

        except zipfile.BadZipFile:
            continue


###############################################################################
# 6. Parse and audit candidates
###############################################################################

audit_rows = []
parsed_candidates = []

for source in candidate_sources:

    try:
        parsed = parse_candidate_bytes(
            source["content"],
            str(source["suffix"]),
        )

        parsed = parsed.loc[
            parsed["program"].isin(
                EXPECTED_PROGRAMS
            )
        ].drop_duplicates()

        counts = (
            parsed.groupby(
                "program"
            )["gene"]
            .nunique()
            .to_dict()
        )

        found_programs = sorted(
            counts
        )

        missing_programs = [
            program
            for program in EXPECTED_PROGRAMS
            if program not in counts
        ]

        small_programs = [
            program
            for program, count
            in counts.items()
            if count < 3
        ]

        complete = bool(
            not missing_programs
            and not small_programs
            and len(found_programs) == 26
        )

        content_hash = (
            canonical_table_hash(
                parsed
            )
            if len(parsed) > 0
            else None
        )

        audit_rows.append(
            {
                "source_type":
                    source["source_type"],

                "source_path":
                    source["source_path"],

                "source_name":
                    source["source_name"],

                "source_sha256":
                    sha256_bytes(
                        source["content"]
                    ),

                "parsed_rows":
                    len(parsed),

                "programs_found":
                    len(found_programs),

                "missing_programs":
                    ";".join(
                        missing_programs
                    ),

                "programs_with_lt3_genes":
                    ";".join(
                        small_programs
                    ),

                "complete_exact_26E":
                    complete,

                "canonical_content_sha256":
                    content_hash,

                "parse_error":
                    "",
            }
        )

        if complete:
            parsed_candidates.append(
                {
                    **source,
                    "table": parsed,
                    "canonical_hash":
                        content_hash,
                }
            )

    except Exception as error:
        audit_rows.append(
            {
                "source_type":
                    source["source_type"],

                "source_path":
                    source["source_path"],

                "source_name":
                    source["source_name"],

                "source_sha256":
                    sha256_bytes(
                        source["content"]
                    ),

                "parsed_rows": 0,
                "programs_found": 0,
                "missing_programs":
                    ";".join(
                        EXPECTED_PROGRAMS
                    ),
                "programs_with_lt3_genes":
                    "",
                "complete_exact_26E":
                    False,
                "canonical_content_sha256":
                    "",
                "parse_error":
                    repr(error),
            }
        )


audit_table = pd.DataFrame(
    audit_rows
)

audit_table.to_csv(
    AUDIT_ROOT
    / "32B3B1A_CATALOG_CANDIDATE_AUDIT.csv",
    index=False,
)


###############################################################################
# 7. Expected-program table
###############################################################################

expected_rows = []

for index, program in enumerate(
    EXPECTED_PROGRAMS,
    start=1,
):
    metadata = PROGRAM_METADATA[
        program
    ]

    expected_rows.append(
        {
            "order": index,
            "program": program,
            **metadata,
        }
    )

expected_table = pd.DataFrame(
    expected_rows
)

expected_table.to_csv(
    AUDIT_ROOT
    / "32B3B1A_EXPECTED_26E_PROGRAMS.csv",
    index=False,
)


###############################################################################
# 8. Resolve exact catalog
###############################################################################

unique_content_hashes = sorted(
    {
        candidate[
            "canonical_hash"
        ]
        for candidate
        in parsed_candidates
    }
)

if len(unique_content_hashes) == 0:
    status = (
        "BLOCKED_NO_COMPLETE_26E_GENESET_CATALOG"
    )

    selected_candidate = None

elif len(unique_content_hashes) > 1:
    status = (
        "BLOCKED_CONFLICTING_26E_GENESET_CATALOGS"
    )

    selected_candidate = None

else:
    status = (
        "READY_EXACT_26E_CATALOG_FROZEN"
    )

    same_content_candidates = [
        candidate
        for candidate
        in parsed_candidates
        if candidate[
            "canonical_hash"
        ] == unique_content_hashes[0]
    ]

    same_content_candidates.sort(
        key=lambda candidate: (
            "SENESCENCE_GENESET_CATALOG"
            not in str(
                candidate["source_name"]
            ).upper(),

            candidate["source_type"]
            != "file",

            len(
                str(
                    candidate[
                        "source_path"
                    ]
                )
            ),
        )
    )

    selected_candidate = (
        same_content_candidates[0]
    )


###############################################################################
# 9. Freeze selected catalog
###############################################################################

selected_metadata = None

if selected_candidate is not None:

    canonical_long = (
        selected_candidate[
            "table"
        ]
        .drop_duplicates()
        .sort_values(
            [
                "program",
                "gene",
            ]
        )
        .reset_index(drop=True)
    )

    canonical_long = canonical_long.merge(
        expected_table.drop(
            columns=[
                "order",
            ]
        ),
        on="program",
        how="left",
        validate="many_to_one",
    )

    canonical_long = canonical_long[
        [
            "program",
            "analysis_class",
            "evidence_layer",
            "expected_direction",
            "direction_multiplier",
            "gene",
        ]
    ]

    canonical_long.to_csv(
        FREEZE_ROOT
        / "32B3B1A_FROZEN_26E_GENESETS_LONG.csv",
        index=False,
    )

    gene_count_table = (
        canonical_long.groupby(
            [
                "program",
                "analysis_class",
                "evidence_layer",
                "expected_direction",
                "direction_multiplier",
            ],
            as_index=False,
        )["gene"]
        .nunique()
        .rename(
            columns={
                "gene":
                    "unique_gene_count"
            }
        )
    )

    gene_count_table.to_csv(
        FREEZE_ROOT
        / "32B3B1A_FROZEN_26E_GENE_COUNTS.csv",
        index=False,
    )

    source_copy_suffix = str(
        selected_candidate[
            "suffix"
        ]
    )

    source_copy_path = (
        FREEZE_ROOT
        / (
            "32B3B1A_SELECTED_SOURCE"
            + source_copy_suffix
        )
    )

    source_copy_path.write_bytes(
        selected_candidate[
            "content"
        ]
    )

    selected_metadata = {
        "source_type":
            selected_candidate[
                "source_type"
            ],

        "source_path":
            selected_candidate[
                "source_path"
            ],

        "source_name":
            selected_candidate[
                "source_name"
            ],

        "source_sha256":
            sha256_bytes(
                selected_candidate[
                    "content"
                ]
            ),

        "canonical_content_sha256":
            selected_candidate[
                "canonical_hash"
            ],

        "duplicate_identical_sources":
            len(
                same_content_candidates
            ),

        "program_count": 26,

        "total_program_gene_rows":
            len(canonical_long),
    }


###############################################################################
# 10. Lock and report
###############################################################################

lock = {
    "stage": "32B3B1A",
    "created_at": now_string(),
    "status": status,
    "expected_program_count": 26,
    "complete_candidates_found":
        len(parsed_candidates),
    "unique_canonical_catalogs":
        len(unique_content_hashes),
    "selected_catalog":
        selected_metadata,
    "new_disease_tests_performed":
        False,
    "enrichment_performed":
        False,
    "next_stage": (
        "32B3B1B_GENERAL_PATHWAY_AND_"
        "FROZEN_26E_OVERLAP"
        if status
        == "READY_EXACT_26E_CATALOG_FROZEN"
        else "RESOLVE_GENESET_CATALOG_PROVENANCE"
    ),
}

with (
    OUTPUT_ROOT
    / "32B3B1A_CATALOG_LOCK.json"
).open(
    "w",
    encoding="utf-8",
) as handle:
    json.dump(
        lock,
        handle,
        indent=2,
        ensure_ascii=False,
    )


report_lines = [
    "# 32B-3B1A Frozen 26E Gene-set Catalog Audit",
    "",
    f"## Status: {status}",
    "",
    f"- Candidate sources examined: {len(candidate_sources)}",
    f"- Complete 26E candidates: {len(parsed_candidates)}",
    (
        "- Unique normalized catalog contents: "
        f"{len(unique_content_hashes)}"
    ),
    "",
]

if selected_metadata is not None:
    report_lines.extend(
        [
            "## Selected source",
            "",
            f"- Source: `{selected_metadata['source_path']}`",
            (
                "- Source SHA-256: "
                f"`{selected_metadata['source_sha256']}`"
            ),
            (
                "- Canonical content SHA-256: "
                f"`{selected_metadata['canonical_content_sha256']}`"
            ),
            "- Programs: 26",
            "",
            "The exact frozen catalog is ready for formal program overlap.",
        ]
    )

else:
    report_lines.extend(
        [
            "## Blocking issue",
            "",
            (
                "No unique complete 26E gene-set catalog could be "
                "frozen. No enrichment or biological overlap should "
                "be run until the source catalog is resolved."
            ),
        ]
    )

(
    OUTPUT_ROOT
    / "32B3B1A_EXECUTIVE_SUMMARY.md"
).write_text(
    "\n".join(report_lines),
    encoding="utf-8",
)


###############################################################################
# 11. Manifest and ZIP
###############################################################################

output_files = [
    path
    for path in OUTPUT_ROOT.rglob("*")
    if path.is_file()
]

manifest_rows = []

for path in output_files:
    manifest_rows.append(
        {
            "relative_path":
                str(
                    path.relative_to(
                        OUTPUT_ROOT
                    )
                ),

            "size_bytes":
                path.stat().st_size,

            "sha256":
                sha256_file(path),
        }
    )

pd.DataFrame(
    manifest_rows
).to_csv(
    OUTPUT_ROOT
    / "32B3B1A_FILE_MANIFEST.csv",
    index=False,
)

zip_path = Path(
    str(OUTPUT_ROOT)
    + ".zip"
)

with zipfile.ZipFile(
    zip_path,
    mode="w",
    compression=zipfile.ZIP_DEFLATED,
    compresslevel=9,
) as archive:

    for path in sorted(
        OUTPUT_ROOT.rglob("*")
    ):
        if path.is_file():
            archive.write(
                path,
                arcname=path.relative_to(
                    OUTPUT_ROOT
                ),
            )

zip_hash = sha256_file(
    zip_path
)

Path(
    str(zip_path)
    + ".sha256.txt"
).write_text(
    f"{zip_hash}  {zip_path.name}\n",
    encoding="ascii",
)


###############################################################################
# 12. Console output
###############################################################################

print()
print("=" * 72)
print("32B-3B1A completed")
print("=" * 72)
print(f"Status: {status}")
print(f"Candidates examined: {len(candidate_sources)}")
print(f"Complete candidates: {len(parsed_candidates)}")
print(f"Unique catalog contents: {len(unique_content_hashes)}")
print()
print("Output:")
print(OUTPUT_ROOT)
print()
print("ZIP:")
print(zip_path)
print()
print("SHA256:")
print(zip_hash)
print("=" * 72)