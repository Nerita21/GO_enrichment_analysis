import csv
import gzip
import shutil
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path

def xml_converter(
    input_xml: str,
    output_tsv: str,
    record_tag: str,
    lookup_tsv: str | None = None,
    xml_key: str | None = None,
    tsv_key: str | None = None
) -> None:
    """
    Convert an XML file into a flattened TSV file.
    Optionally join with a lookup TSV on matching keys.

    Args:
        input_xml (str): Path to the XML file (can be .xml or .xml.gz).
        output_tsv (str): Path where the resulting TSV will be written.
        record_tag (str): XML tag name that defines one record.
        lookup_tsv (str | None): Optional path to lookup TSV file.
        xml_key (str | None): Optional key column name from XML for joining.
        tsv_key (str | None): Optional key column name from lookup TSV for joining.
    """

    # --- Helper: unzip gzipped XML if needed ---
    def _maybe_unzip(file_path):
        p = Path(file_path)
        if p.suffix == ".gz":
            out = p.with_suffix("")
            with gzip.open(p, "rb") as f_in, open(out, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
            return out
        return p

    # --- Helper: flatten XML element ---
    def _flatten_element(elem, prefix="", attr_prefix="@"):
        out = defaultdict(list)
        def _recurse(e, cur_prefix):
            for k, v in e.attrib.items():
                out[f"{cur_prefix}{attr_prefix}{k}"].append(v)
            children = list(e)
            text = (e.text or "").strip()
            if not children:
                out[cur_prefix.rstrip('.')].append(text)
            else:
                if text:
                    out[f"{cur_prefix.rstrip('.')}.#text"].append(text)
                for child in children:
                    _recurse(child, f"{cur_prefix}{child.tag}.")
        _recurse(elem, prefix)
        return out

    # --- Helper: iterate over XML records ---
    def _iter_records(xml_file, tag):
        context = ET.iterparse(xml_file, events=("end",))
        for _, elem in context:
            if elem.tag == tag:
                yield elem
                elem.clear()

    # --- Helper: load lookup TSV ---
    def _load_lookup_tsv(path, key_col):
        lookup = {}
        with open(path, newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                lookup[row[key_col]] = row
        return lookup, reader.fieldnames

    # --- Unzip XML if needed ---
    xml_path = _maybe_unzip(input_xml)

    # --- Discover XML fields ---
    fieldnames = []
    seen = set()
    for rec in _iter_records(xml_path, record_tag):
        flat = _flatten_element(rec)
        for k in flat.keys():
            if k not in seen:
                seen.add(k)
                fieldnames.append(k)
        rec.clear()

    # --- Optional join setup ---
    lookup = {}
    lookup_fields = []
    if lookup_tsv and xml_key and tsv_key:
        lookup, lookup_fields = _load_lookup_tsv(lookup_tsv, tsv_key)
        all_fields = fieldnames + [f for f in lookup_fields if f != tsv_key]
        do_join = True
    else:
        all_fields = fieldnames
        do_join = False

    # --- Write combined TSV ---
    with open(output_tsv, "w", newline="", encoding="utf-8") as f_out:
        writer = csv.DictWriter(f_out, fieldnames=all_fields, delimiter="\t")
        writer.writeheader()

        for rec in _iter_records(xml_path, record_tag):
            flat = _flatten_element(rec)
            row = {f: "|".join(flat.get(f, [])) for f in fieldnames}
            if do_join:
                key_val = row.get(xml_key, "")
                if key_val in lookup:
                    for k, v in lookup[key_val].items():
                        if k != tsv_key:
                            row[k] = v
            writer.writerow(row)
            rec.clear()

    print(f"Conversion complete. Output saved to: {output_tsv}")
    if do_join:
        print(f"Joined with lookup table: {lookup_tsv}")