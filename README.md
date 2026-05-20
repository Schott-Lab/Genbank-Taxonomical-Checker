# Genbank_Taxonomical_Check

A command-line tool that groups [GenBank](https://www.ncbi.nlm.nih.gov/genbank/) (`.gb`) records into separate files based on whether their taxonomic lineage matches terms in a configurable wordlist.

Built on [Biopython](https://biopython.org/), the script reads each record's taxonomy and organism annotations, assigns the record to the first matching term in the wordlist, and writes per-term GenBank files alongside a summary report. Records that match nothing are collected separately so none are silently dropped.

---

## Features

- **Wordlist-driven grouping** — each record is sorted into the first term it matches, in wordlist order.
- **Two matching styles** — plain terms match the term itself; prefixed terms (`Prefix;taxonname`) match only the part after the semicolon, letting you annotate entries while keeping the search precise.
- **No-match capture** — unmatched records are written to a dedicated `NO-MATCH.gb` file rather than discarded.
- **Per-file output folders** — results for each input live in their own directory, keeping batch runs tidy.
- **Summary report** — a plain-text report lists every term, whether it matched, and how many records it captured.
- **Readable console output** — matched terms are highlighted in bold so hits stand out at a glance.
- **Batch mode** — process a single file or every `.gb` file in the current directory.

---

## Requirements

- Python 3.6 or newer
- [Biopython](https://biopython.org/)

Install the dependency with:

```bash
pip install biopython
```

---

## Installation

No packaging step is required. Download `Genbank_Taxonomical_Check.py` and run it with the Python interpreter.

---

## Usage

```bash
python3 Genbank_Taxonomical_Check.py <input.gb>
python3 Genbank_Taxonomical_Check.py -a        # or --all
python3 Genbank_Taxonomical_Check.py -h        # or --help
```

| Argument | Description |
| --- | --- |
| `<input.gb>` | Process a single GenBank file. |
| `-a`, `--all` | Process every `.gb` file in the current directory. |
| `-h`, `--help` | Show usage and exit. |

The script requires exactly one argument.

---

## Configuring the wordlist

The grouping terms are defined by the `taxonomy_wordlist` list near the top of the script. This is the **only** part intended to be edited. The script ships with a default list covering major vertebrate groups (cartilaginous and bony fishes, amphibians, reptiles, birds, and mammals), with reptile clades broken out by family.

Two term formats are supported:

**Plain string** — searched directly against the lineage.

```python
"Bacteria"
```
Matches `Bacteria` in the lineage. Output file: `inputname_Bacteria.gb`

**Prefixed string** — only the part *after* the semicolon is searched. The prefix is purely for your own organization (e.g. grouping families under a clade).

```python
"Infraorder;taxonname"
```
Searches for `taxonname` in the lineage. Output file: `inputname_Infraorder_taxonname.gb` (the semicolon becomes an underscore in the filename).

Example:

```python
taxonomy_wordlist = [
    "Bacteria",
    "Infraorder1;taxonname1",
    "TagID2;taxonname2",
]
```

> **Note:** Matching is case-insensitive and compares against whole lineage components (each rank and the organism name), not substrings within a name. The wordlist must not be empty, or the script exits with an error.

---

## How matching works

For each record, the script builds a lineage set from the record's `taxonomy` annotation plus its `organism` name. Each wordlist term is checked in order, and the record is assigned to the **first** term whose search value appears as a component of that lineage. Once matched, the record is not considered for later terms. Records matching no term go to the no-match group.

---

## Output

For each input file, a folder named `<basename>_taxonsort/` is created containing:

| File | Contents |
| --- | --- |
| `<basename>_<term>.gb` | One file per wordlist term, holding the records assigned to it. (Semicolons in prefixed terms become underscores.) |
| `<basename>_NO-MATCH.gb` | Records that matched no term (written only when such records exist). |
| `result_report_<basename>.txt` | A text summary: each term, a `Y`/`-` match flag, and the record count. |

A per-term file is created for every wordlist entry; entries with no matches produce an empty file.

---

## Example console output

```
sample.gb
Y : Teleostei                           42
- : Coelacanthiformes
Y : Aves                                 7
...
- : NO-MATCH
```

Lines beginning with `Y :` (matches) are shown in bold in the terminal.

---

## License

MIT License.

## Author

Arshia Farajollahi
