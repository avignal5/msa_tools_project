from pathlib import Path
from typing import List, Dict, Optional


class GFFFeature:
    """
    Simple representation of a GFF3 feature.
    """

    def __init__(self, seqid: str, source: str, feature_type: str,
                 start: int, end: int, strand: str, attributes: Dict[str, str]):
        self.seqid = seqid
        self.source = source
        self.type = feature_type
        self.start = start
        self.end = end
        self.strand = strand
        self.attributes = attributes

    def __repr__(self):
        return f"<GFFFeature {self.type} {self.start}-{self.end}>"


def parse_attributes(attr_string: str) -> Dict[str, str]:
    attributes = {}
    for field in attr_string.split(";"):
        if "=" in field:
            key, value = field.split("=", 1)
            attributes[key] = value
    return attributes


def read_gff3(filepath: Path,
              feature_type: Optional[str] = None) -> List[GFFFeature]:
    """
    Parse a GFF3 file.

    Parameters
    ----------
    filepath : Path
    feature_type : str, optional
        Filter by feature type (e.g. 'gene', 'CDS')

    Returns
    -------
    List[GFFFeature]
    """

    filepath = Path(filepath)
    features = []

    with open(filepath) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue

            parts = line.strip().split("\t")
            if len(parts) != 9:
                continue

            seqid, source, ftype, start, end, score, strand, phase, attributes = parts

            if feature_type and ftype != feature_type:
                continue

            feature = GFFFeature(
                seqid=seqid,
                source=source,
                feature_type=ftype,
                start=int(start),
                end=int(end),
                strand=strand,
                attributes=parse_attributes(attributes),
            )

            features.append(feature)

    return features