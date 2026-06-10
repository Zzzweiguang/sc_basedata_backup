#!/usr/bin/env python3
from __future__ import annotations

import re
import xml.etree.ElementTree as ET
from pathlib import Path


BASE = Path("/home/lixiangyu/zr/Annotate/ANNOTATE_new/10_downstream/cell_type_rename/test_0530/mouse_literature_celltypes")
XML_DIR = BASE / "pmc_xml"
OUT_DIR = BASE / "article_text"


def clean(text: str) -> str:
    return re.sub(r"\s+", " ", text).strip()


def node_text(node: ET.Element) -> str:
    return clean(" ".join(t for t in node.itertext() if t and t.strip()))


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for xml_path in sorted(XML_DIR.glob("GSE*_PMC*.xml")):
        out_path = OUT_DIR / f"{xml_path.stem}.txt"
        if out_path.exists():
            continue
        try:
            root = ET.parse(xml_path).getroot()
        except ET.ParseError:
            print(f"SKIP non-XML/invalid XML: {xml_path}")
            continue
        chunks: list[str] = []
        for xpath in [
            ".//article-title",
            ".//abstract",
            ".//body//sec",
            ".//fig",
            ".//table-wrap",
            ".//supplementary-material",
        ]:
            for node in root.findall(xpath):
                text = node_text(node)
                if text:
                    chunks.append(text)
        out_path.write_text("\n\n".join(dict.fromkeys(chunks)) + "\n")
        print(out_path)


if __name__ == "__main__":
    main()
