#!/usr/bin/env python3
from __future__ import annotations

import csv
import re
import xml.etree.ElementTree as ET
from pathlib import Path


BASE = Path("/home/lixiangyu/zr/Annotate/ANNOTATE_new/10_downstream/cell_type_rename/test_0530/mouse_literature_celltypes")
GEO_DIR = BASE / "geo_text"
PUBMED_XML = BASE / "pubmed" / "pubmed_metadata.xml"


def clean(value: str | None) -> str:
    if not value:
        return ""
    return re.sub(r"\s+", " ", value).strip()


def parse_pubmed() -> dict[str, dict[str, str]]:
    root = ET.parse(PUBMED_XML).getroot()
    records: dict[str, dict[str, str]] = {}
    for article in root.findall(".//PubmedArticle"):
        pmid = clean(article.findtext("./MedlineCitation/PMID"))
        title_node = article.find("./MedlineCitation/Article/ArticleTitle")
        title = clean("".join(title_node.itertext()) if title_node is not None else "")
        journal = clean(
            article.findtext("./MedlineCitation/Article/Journal/Title")
            or article.findtext("./MedlineCitation/Article/Journal/ISOAbbreviation")
        )
        year = clean(
            article.findtext("./MedlineCitation/Article/Journal/JournalIssue/PubDate/Year")
            or article.findtext("./MedlineCitation/Article/ArticleDate/Year")
        )
        doi = ""
        pmcid = ""
        for aid in article.findall("./PubmedData/ArticleIdList/ArticleId"):
            id_type = aid.attrib.get("IdType", "")
            if id_type == "doi":
                doi = clean(aid.text)
            elif id_type == "pmc":
                pmcid = clean(aid.text)
        records[pmid] = {
            "pmid": pmid,
            "year": year,
            "journal": journal,
            "doi": doi,
            "pmcid": pmcid,
            "title": title,
            "pubmed_url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/" if pmid else "",
        }
    return records


def parse_geo_file(path: Path) -> dict[str, str]:
    text = path.read_text(errors="ignore")
    fields = {
        "GSE": path.stem,
        "Series_title": "",
        "PMIDs_from_GEO": "",
        "DOI_or_web_link_from_GEO": "",
        "Platform_organism": "",
        "Sample_organism": "",
        "Summary": "",
        "Overall_design": "",
        "GEO_url": f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={path.stem}",
    }
    pmids: list[str] = []
    platform_orgs: list[str] = []
    sample_orgs: list[str] = []
    summaries: list[str] = []
    designs: list[str] = []
    web_links: list[str] = []
    for line in text.splitlines():
        if line.startswith("!Series_title = "):
            fields["Series_title"] = clean(line.split("=", 1)[1])
        elif line.startswith("!Series_pubmed_id = "):
            pmids.append(clean(line.split("=", 1)[1]))
        elif line.startswith("!Series_web_link = "):
            web_links.append(clean(line.split("=", 1)[1]))
        elif line.startswith("!Series_platform_organism = "):
            platform_orgs.append(clean(line.split("=", 1)[1]))
        elif line.startswith("!Series_sample_organism = "):
            sample_orgs.append(clean(line.split("=", 1)[1]))
        elif line.startswith("!Series_summary = "):
            summaries.append(clean(line.split("=", 1)[1]))
        elif line.startswith("!Series_overall_design = "):
            designs.append(clean(line.split("=", 1)[1]))
    fields["PMIDs_from_GEO"] = ";".join(pmids)
    fields["DOI_or_web_link_from_GEO"] = ";".join(web_links)
    fields["Platform_organism"] = ";".join(dict.fromkeys(platform_orgs))
    fields["Sample_organism"] = ";".join(dict.fromkeys(sample_orgs))
    fields["Summary"] = " ".join(summaries)
    fields["Overall_design"] = " ".join(designs)
    return fields


def main() -> None:
    pubmed = parse_pubmed()
    rows = []
    for geo_path in sorted(GEO_DIR.glob("GSE*.txt")):
        row = parse_geo_file(geo_path)
        primary_pmids = [p for p in row["PMIDs_from_GEO"].split(";") if p]
        primary = pubmed.get(primary_pmids[0], {}) if primary_pmids else {}
        row.update(
            {
                "Primary_PMID": primary.get("pmid", ""),
                "Primary_title": primary.get("title", ""),
                "Primary_year": primary.get("year", ""),
                "Primary_journal": primary.get("journal", ""),
                "Primary_DOI": primary.get("doi", ""),
                "Primary_PMCID": primary.get("pmcid", ""),
                "Primary_PubMed_url": primary.get("pubmed_url", ""),
            }
        )
        rows.append(row)

    out = BASE / "mouse_geo_pubmed_metadata.tsv"
    fieldnames = [
        "GSE",
        "Series_title",
        "PMIDs_from_GEO",
        "Primary_PMID",
        "Primary_title",
        "Primary_year",
        "Primary_journal",
        "Primary_DOI",
        "Primary_PMCID",
        "DOI_or_web_link_from_GEO",
        "Platform_organism",
        "Sample_organism",
        "Summary",
        "Overall_design",
        "GEO_url",
        "Primary_PubMed_url",
    ]
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
