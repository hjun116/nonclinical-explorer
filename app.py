"""
Nonclinical Drug Explorer
─────────────────────────
APIs used (all free, no key required):
  • PubMed E-utilities  — https://eutils.ncbi.nlm.nih.gov/
  • ChEMBL REST API     — https://www.ebi.ac.uk/chembl/api/data/

Classification: keyword rule-based (no external API)
Deploy: Streamlit Cloud — push to GitHub, connect repo, done.
"""

import streamlit as st
import requests
import xml.etree.ElementTree as ET
import re
from typing import Optional

# ── Page config ───────────────────────────────────────────────
st.set_page_config(
    page_title="Nonclinical Drug Explorer",
    page_icon="🔬",
    layout="wide",
    initial_sidebar_state="collapsed",
)

# ── Minimal custom CSS ─────────────────────────────────────────
st.markdown("""
<style>
  .block-container { padding-top: 2rem; max-width: 960px; }
  .stTabs [data-baseweb="tab"] { font-size: 14px; }
  div[data-testid="metric-container"] { background: #f7f7f5; border-radius: 8px; padding: 12px; }
  .paper-card {
    border: 1px solid #e8e8e4;
    border-radius: 10px;
    padding: 1rem 1.25rem;
    margin-bottom: 10px;
    background: #fff;
  }
  .cat-badge {
    display: inline-block;
    font-size: 11px;
    padding: 2px 8px;
    border-radius: 4px;
    font-weight: 500;
  }
  .cat-safety  { background:#fde8e8; color:#8b1a1a; }
  .cat-pk      { background:#fef3d8; color:#7a4e00; }
  .cat-pd      { background:#ddeeff; color:#0a4080; }
  .cat-eff     { background:#e6f5d8; color:#2a6010; }
  .cat-other   { background:#f0efeb; color:#5a5a56; }
  .conf-low    { background:#fff3d0; color:#7a5000; font-size:11px; padding:2px 7px; border-radius:4px; }
  .disclaimer  { font-size:12px; color:#9a9a94; border-left:3px solid #e0e0d8; padding-left:10px; margin-bottom:1rem; }
</style>
""", unsafe_allow_html=True)


# ══════════════════════════════════════════════════════════════
# CHEMBL RESOLVER
# ══════════════════════════════════════════════════════════════
@st.cache_data(ttl=3600, show_spinner=False)
def resolve_chembl(query: str) -> Optional[dict]:
    """Resolve drug name → ChEMBL data including all synonyms."""
    base = "https://www.ebi.ac.uk/chembl/api/data"
    headers = {"Accept": "application/json"}

    # 1. Exact preferred name
    try:
        r = requests.get(
            f"{base}/molecule",
            params={"pref_name__iexact": query, "format": "json"},
            headers=headers, timeout=10
        )
        data = r.json()
        mol = (data.get("molecules") or [None])[0]

        # 2. Synonym search fallback
        if not mol:
            r = requests.get(
                f"{base}/molecule/search",
                params={"q": query, "format": "json"},
                headers=headers, timeout=10
            )
            data = r.json()
            mol = (data.get("molecules") or [None])[0]

        if not mol:
            return None

        synonyms = [s["synonym"] for s in (mol.get("molecule_synonyms") or []) if s.get("synonym")]
        props = mol.get("molecule_properties") or {}

        return {
            "chembl_id":  mol.get("molecule_chembl_id", ""),
            "inn":        mol.get("pref_name", query),
            "type":       mol.get("molecule_type", ""),
            "formula":    props.get("full_molformula", ""),
            "mw":         props.get("full_mwt", ""),
            "aliases":    list(dict.fromkeys([mol.get("pref_name", "")] + synonyms)),
        }
    except Exception:
        return None


# ══════════════════════════════════════════════════════════════
# PUBMED SEARCH
# ══════════════════════════════════════════════════════════════
NONCLINICAL_TERMS = (
    "preclinical[Title/Abstract] OR nonclinical[Title/Abstract] OR "
    "toxicology[Title/Abstract] OR pharmacokinetics[Title/Abstract] OR "
    "\"animal model\"[Title/Abstract] OR \"in vivo\"[Title/Abstract] OR "
    "ADME[Title/Abstract] OR \"drug metabolism\"[Title/Abstract] OR "
    "\"safety pharmacology\"[Title/Abstract]"
)


def build_pubmed_query(query: str, mode: str, aliases: list[str]) -> str:
    if mode == "Drug name":
        all_names = list(dict.fromkeys([query] + (aliases or [])))[:6]
        name_terms = " OR ".join(f'"{n}"[Title/Abstract]' for n in all_names)
        return f"({name_terms}) AND ({NONCLINICAL_TERMS})"
    elif mode == "Company":
        return f'"{query}"[Affiliation] AND ({NONCLINICAL_TERMS})'
    else:  # Indication
        return f'"{query}"[MeSH Terms] AND ({NONCLINICAL_TERMS})'


@st.cache_data(ttl=1800, show_spinner=False)
def search_pubmed(query: str, mode: str, aliases: list[str], max_results: int = 40, min_year: int = 2015) -> list[dict]:
    pm_query = build_pubmed_query(query, mode, aliases)
    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

    # Search — relevance-sorted, filtered to min_year onwards
    r = requests.get(f"{base}/esearch.fcgi", params={
        "db":       "pubmed",
        "term":     pm_query,
        "retmax":   max_results,
        "retmode":  "json",
        "sort":     "relevance",
        "datetype": "pdat",
        "mindate":  str(min_year),
        "maxdate":  "3000",
    }, timeout=15)
    ids = r.json().get("esearchresult", {}).get("idlist", [])
    if not ids:
        return []

    # Fetch
    r = requests.get(f"{base}/efetch.fcgi", params={
        "db": "pubmed", "id": ",".join(ids), "retmode": "xml"
    }, timeout=20)
    return parse_pubmed_xml(r.text)


def parse_pubmed_xml(xml_text: str) -> list[dict]:
    root = ET.fromstring(xml_text)
    papers = []
    for i, article in enumerate(root.findall(".//PubmedArticle")):
        pmid  = getattr(article.find(".//PMID"), "text", "")
        title = getattr(article.find(".//ArticleTitle"), "text", "No title") or "No title"

        abstract_parts = [el.text or "" for el in article.findall(".//AbstractText")]
        abstract = " ".join(abstract_parts).strip() or "No abstract available"

        author_nodes = article.findall(".//Author")[:3]
        authors = []
        for au in author_nodes:
            last = getattr(au.find("LastName"), "text", "")
            init = getattr(au.find("Initials"), "text", "")
            if last:
                authors.append(f"{last} {init}".strip())
        if len(article.findall(".//Author")) > 3:
            authors.append("et al.")

        journal = (
            getattr(article.find(".//ISOAbbreviation"), "text", "") or
            getattr(article.find(".//Title"), "text", "")
        )
        year = (
            getattr(article.find(".//PubDate/Year"), "text", "") or
            (getattr(article.find(".//MedlineDate"), "text", "") or "")[:4]
        )

        papers.append({
            "pmid":        pmid,
            "title":       title,
            "abstract":    abstract,
            "authors":     ", ".join(authors),
            "journal":     journal,
            "year":        year,
            "raw_index":   i,
            "category":    "other",
            "confidence":  "low",
            "key_findings": [],
            "source_db":   "PubMed",
        })
    return papers


# ══════════════════════════════════════════════════════════════
# EUROPE PMC SEARCH
# ══════════════════════════════════════════════════════════════
EPMC_NONCLINICAL_KW = (
    "(preclinical OR nonclinical OR toxicology OR pharmacokinetics OR "
    "\"animal model\" OR \"in vivo\" OR ADME OR \"drug metabolism\" OR "
    "\"safety pharmacology\")"
)


@st.cache_data(ttl=1800, show_spinner=False)
def search_europe_pmc(query: str, mode: str, aliases: list[str], max_results: int = 40, min_year: int = 2015) -> list[dict]:
    """Search Europe PMC — newest first, filtered to min_year onwards."""
    base = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"

    if mode == "Drug name":
        all_names = list(dict.fromkeys([query] + (aliases or [])))[:6]
        name_terms = " OR ".join(f'"{n}"' for n in all_names)
        q = f"({name_terms}) AND {EPMC_NONCLINICAL_KW} AND (PUB_YEAR:[{min_year} TO 3000])"
    elif mode == "Company":
        q = f'AFFILIATION:"{query}" AND {EPMC_NONCLINICAL_KW} AND (PUB_YEAR:[{min_year} TO 3000])'
    else:
        q = f'"{query}" AND {EPMC_NONCLINICAL_KW} AND (PUB_YEAR:[{min_year} TO 3000])'

    try:
        r = requests.get(base, params={
            "query":      q,
            "resultType": "core",
            "pageSize":   max_results,
            "format":     "json",
            "sort":       "P_PDATE_D desc",   # newest first
        }, timeout=15)
        results = r.json().get("resultList", {}).get("result", [])
        return parse_europe_pmc(results)
    except Exception:
        return []


def parse_europe_pmc(results: list[dict]) -> list[dict]:
    papers = []
    for i, item in enumerate(results):
        pmid     = item.get("pmid", "")
        title    = item.get("title", "No title").rstrip(".")
        abstract = item.get("abstractText", "") or "No abstract available"

        authors_list = (item.get("authorList") or {}).get("author", [])
        authors = [
            f'{a.get("lastName", "")} {a.get("initials", "")}'.strip()
            for a in authors_list[:3] if a.get("lastName")
        ]
        if len(authors_list) > 3:
            authors.append("et al.")

        journal = item.get("journalAbbreviation") or item.get("journalTitle", "")
        year    = str(item.get("pubYear", ""))

        papers.append({
            "pmid":        pmid,
            "title":       title,
            "abstract":    abstract,
            "authors":     ", ".join(authors),
            "journal":     journal,
            "year":        year,
            "raw_index":   i,
            "category":    "other",
            "confidence":  "low",
            "key_findings": [],
            "source_db":   "Europe PMC",
        })
    return papers


def merge_and_deduplicate(pubmed_papers: list[dict], epmc_papers: list[dict]) -> list[dict]:
    """Merge results: PubMed (relevance) first, Europe PMC (newest) fills gaps.
    Deduplication by PMID. Papers without PMID (preprints) are always included."""
    seen_pmids: set[str] = set()
    merged = []

    for p in pubmed_papers:
        pid = p.get("pmid", "")
        if pid:
            seen_pmids.add(pid)
        merged.append(p)

    for p in epmc_papers:
        pid = p.get("pmid", "")
        if not pid or pid not in seen_pmids:
            if pid:
                seen_pmids.add(pid)
            p["raw_index"] = len(merged)
            merged.append(p)

    return merged


# ══════════════════════════════════════════════════════════════
# KEYWORD CLASSIFIER
# ══════════════════════════════════════════════════════════════
RULES = {
    "safety": {
        "title": ["toxic","toxicolog","genotox","mutagenic","carcinogen",
                  "safety pharmacol","herg","hepatotox","nephrotox",
                  "repeated dose","acute tox","subchronic","micronucleus",
                  "ames test","reproductive tox","teratogen","noael","loael"],
        "abstract": ["toxic","safety","genotox","noael","loael","adverse",
                     "ld50","herg","hepatotox","nephrotox","organ tox"]
    },
    "pk": {
        "title": ["pharmacokinetic","adme","bioavailability","absorption",
                  "metabolism","excretion","half-life","clearance","auc",
                  "cmax","drug interaction","cyp","p-glycoprotein",
                  "plasma protein binding","tissue distribution"],
        "abstract": ["pharmacokinetic","adme","bioavailability","half-life",
                     "clearance","auc","cmax","metabolism","cyp","absorption"]
    },
    "pd": {
        "title": ["pharmacodynamic","mechanism of action","moa","binding affinity",
                  "ic50","selectivity","inhibition","agonist","antagonist",
                  "kinase inhibit","target engagement","receptor binding"],
        "abstract": ["pharmacodynamic","mechanism","binding affinity","ic50",
                     "selectivity","inhibit","receptor","agonist","antagonist"]
    },
    "efficacy": {
        "title": ["animal model","mouse model","rat model","in vivo","xenograft",
                  "tumor model","disease model","efficacy","antitumor",
                  "tumor growth","dose-response","ed50","survival benefit"],
        "abstract": ["animal model","in vivo","xenograft","efficacy","tumor",
                     "disease model","antitumor","survival","dose-response"]
    },
}

SPECIES_PATTERNS = [
    (r"\b(mice|mouse|murine)\b",                          "Mouse"),
    (r"\b(rat|rats|rodent)\b",                            "Rat"),
    (r"\b(monkey|cynomolgus|rhesus|primate|NHP)\b",       "Monkey"),
    (r"\b(dog|dogs|beagle|canine)\b",                     "Dog"),
    (r"\b(rabbit|rabbits)\b",                             "Rabbit"),
    (r"\b(hamster|guinea pig|minipig|swine|pig)\b",       "Other"),
]

# ── Category-specific extraction patterns ─────────────────────
SAFETY_PATTERNS = [
    ("NOAEL/NOEL",   r"(?:noael|noel|no.observed.adverse.effect.level)[^a-z\d]{0,15}(\d[\d.,]*\s*mg\/kg[^\s,;.]*)"),
    ("LOAEL",        r"(?:loael|lowest.observed.adverse.effect.level)[^a-z\d]{0,15}(\d[\d.,]*\s*mg\/kg[^\s,;.]*)"),
    ("LD50",         r"(?:ld50|lethal dose)[^a-z\d]{0,10}(\d[\d.,]*\s*(?:mg\/kg|μg\/kg)[^\s,;.]*)"),
    ("Duration",     r"(\d+[.\-]?(?:day|week|month)[s]?(?:\s+(?:repeat|repeated|oral|iv|dosing))?)"),
    ("Tox finding",  r"\b(hepatotox\w*|nephrotox\w*|cardiotox\w*|neurotox\w*|genotox\w*|hepatic\s+(?:injury|damage|failure))\b"),
]

PK_PATTERNS = [
    ("Cmax",   r"(?:cmax|c\s*max|peak\s+concentration)[^a-z\d]{0,10}(\d[\d.,]*\s*(?:ng|μg|µg)[\s·/]?(?:m[Ll]|L)[^\s,;.]*)"),
    ("AUC",    r"(?:auc(?:0.∞|0-inf|0-t|last)?)[^a-z\d]{0,10}(\d[\d.,]*\s*(?:ng|μg|µg)[\s·•/]?h[\s·•/]?(?:m[Ll]|L)[^\s,;.]*)"),
    ("t½",     r"(?:half.life|t½|t1\/2)[^a-z\d]{0,10}(\d[\d.,]*\s*h(?:ours?)?)"),
    ("F%",     r"(?:bioavailability|oral\s+(?:bioavailability|f))[^a-z\d]{0,10}(\d[\d.,]*\s*%)"),
    ("CL",     r"(?:clearance|cl)[^a-z\d]{0,10}(\d[\d.,]*\s*(?:ml|L)[\s·/](?:min|h)[\s·/]?(?:kg)?[^\s,;.]*)"),
]

PD_PATTERNS = [
    ("IC50",       r"(?:ic50)[^a-z\d]{0,10}(\d[\d.,]*\s*(?:n[Mm]|μ[Mm]|µ[Mm]|p[Mm])[^\s,;.]*)"),
    ("EC50",       r"(?:ec50)[^a-z\d]{0,10}(\d[\d.,]*\s*(?:n[Mm]|μ[Mm]|µ[Mm])[^\s,;.]*)"),
    ("Ki",         r"(?:\bki\b)[^a-z\d]{0,10}(\d[\d.,]*\s*(?:n[Mm]|μ[Mm]|µ[Mm])[^\s,;.]*)"),
    ("Selectivity",r"(?:selectivity|selective)[^a-z\d]{0,15}(\d[\d.,]*[\s-]?fold[^\s,;.]*)"),
]

EFFICACY_PATTERNS = [
    ("Model",        r"\b(xenograft|PDX|syngeneic|orthotopic|allograft|transgenic|knock.?out|induced\s+model)\b"),
    ("Tumor inhib.", r"(?:tumor|tumour)\s+(?:growth\s+)?(?:inhibition|reduction|suppression)[^a-z\d]{0,10}(\d[\d.,]*\s*%)"),
    ("Survival",     r"(?:survival|overall\s+survival)[^a-z\d]{0,10}(\d[\d.,]*\s*(?:%|days?|weeks?)[^\s,;.]*)"),
    ("ED50/TGI",     r"(?:ed50|tgi|tumor\s+growth\s+inhibition)[^a-z\d]{0,10}(\d[\d.,]*\s*(?:mg\/kg|%)[^\s,;.]*)"),
]

CAT_PATTERNS = {
    "safety":   SAFETY_PATTERNS,
    "pk":       PK_PATTERNS,
    "pd":       PD_PATTERNS,
    "efficacy": EFFICACY_PATTERNS,
}


def extract_species(text: str) -> Optional[str]:
    for pattern, label in SPECIES_PATTERNS:
        if re.search(pattern, text, re.I):
            return label
    return None


def extract_findings(text: str, category: str, title: str = "") -> list[dict]:
    """Extract category-specific key findings from abstract + title text."""
    findings = []

    # Species — extract for ALL categories using title + abstract combined
    full_text = (title + " " + text).strip()
    sp = extract_species(full_text)
    if sp:
        findings.append({"key": "Species", "val": sp})

    # Category-specific metrics (abstract text only)
    patterns = CAT_PATTERNS.get(category, [])
    for label, pat in patterns:
        m = re.search(pat, text, re.I)
        if m:
            val = (m.group(1) if m.lastindex else m.group(0)).strip()
            val = re.sub(r"\s+", " ", val).strip(" ,;.")
            if val:
                findings.append({"key": label, "val": val})

    return findings


def keyword_classify(paper: dict) -> dict:
    title    = paper["title"].lower()
    abstract = paper["abstract"].lower()
    for cat, rule in RULES.items():
        t_hit = any(kw in title    for kw in rule["title"])
        a_hit = any(kw in abstract for kw in rule["abstract"])
        if t_hit and a_hit: return {"category": cat, "confidence": "high"}
        if t_hit or  a_hit: return {"category": cat, "confidence": "low"}
    return {"category": "other", "confidence": "low"}


# ── In vivo detection ─────────────────────────────────────────
# Positive signals: words that strongly indicate animal experiment
IN_VIVO_POSITIVE = [
    r"\b(mice|mouse|murine|rats?|rodent)\b",
    r"\b(cynomolgus|rhesus|monkey|primate|NHP)\b",
    r"\b(beagle|canine|dogs?)\b",
    r"\b(rabbit|hamster|guinea\s+pig|minipig|swine)\b",
    r"\bin\s+vivo\b",
    r"\banimal\s+(model|study|experiment|stud)\b",
    r"\b(xenograft|allograft|syngeneic|orthotopic)\b",
    r"\b(oral\s+(?:dosing|administration|gavage)|intravenous(?:ly)?|subcutaneous(?:ly)?)\b",
    r"\b(repeated.dose|single.dose|dose.escalation)\s+(?:study|tox)",
    r"\b(pharmacokinetics?\s+in|pk\s+in)\s+(?:rats?|mice|dogs?|monkey)",
]

# Negative signals: words that suggest in vitro / cell-only study
IN_VITRO_NEGATIVE = [
    r"\bcell\s+line\b",
    r"\bin\s+vitro\b",
    r"\bcell\s+culture\b",
    r"\b(hela|mcf.7|a549|hek293|jurkat|thp.?1|hct116)\b",
    r"\b(transfect|siRNA|shRNA|CRISPR)\b",
]


def detect_in_vivo(paper: dict) -> dict:
    """
    Returns {"in_vivo": True/False, "in_vivo_confidence": "confirmed"/"likely"/"unlikely"}.

    Logic:
      - Any positive signal found → candidate
      - If species extracted (from key_findings) → confirmed
      - Negative-only (no positive) → unlikely (in vitro)
      - Mixed (positive + negative) → likely (both, but animal data present)
    """
    text = (paper.get("title", "") + " " + paper.get("abstract", "")).lower()

    pos_hits = sum(1 for pat in IN_VIVO_POSITIVE  if re.search(pat, text, re.I))
    neg_hits = sum(1 for pat in IN_VITRO_NEGATIVE if re.search(pat, text, re.I))

    # Species already extracted = strong confirmation
    has_species = any(f["key"] == "Species" for f in paper.get("key_findings", []))

    if has_species or pos_hits >= 2:
        confidence = "confirmed"
        in_vivo = True
    elif pos_hits == 1 and neg_hits == 0:
        confidence = "likely"
        in_vivo = True
    elif pos_hits >= 1 and neg_hits >= 1:
        confidence = "likely"   # mixed — animal data probably present
        in_vivo = True
    else:
        confidence = "unlikely"
        in_vivo = False

    return {"in_vivo": in_vivo, "in_vivo_confidence": confidence}


def keyword_classify_all(papers: list[dict]) -> list[dict]:
    for p in papers:
        r = keyword_classify(p)
        p["category"]     = r["category"]
        p["confidence"]   = r["confidence"]
        # Extract after category is known so patterns are category-specific
        p["key_findings"] = extract_findings(p["abstract"], p["category"], p.get("title", ""))
        # Attach in_vivo flag (uses key_findings, so must come after)
        iv = detect_in_vivo(p)
        p["in_vivo"]            = iv["in_vivo"]
        p["in_vivo_confidence"] = iv["in_vivo_confidence"]
    return papers


# ══════════════════════════════════════════════════════════════
# UI HELPERS
# ══════════════════════════════════════════════════════════════
CAT_LABELS = {
    "safety":   "Safety / Tox",
    "pk":       "PK / ADME",
    "pd":       "PD / MOA",
    "efficacy": "In vivo Efficacy",
    "other":    "Other",
}
CAT_CSS = {
    "safety": "cat-safety", "pk": "cat-pk",
    "pd": "cat-pd", "efficacy": "cat-eff", "other": "cat-other",
}


def render_paper(p: dict):
    cat_label = CAT_LABELS.get(p["category"], "Other")
    cat_class = CAT_CSS.get(p["category"], "cat-other")
    conf_tag  = '<span style="font-size:11px;background:#fff0e6;color:#c05000;padding:2px 8px;border-radius:4px;">⚠ Verify with source</span>' if p["confidence"] == "low" else ""

    # In vivo badge
    iv_conf = p.get("in_vivo_confidence", "unlikely")
    if iv_conf == "confirmed":
        iv_badge = '<span style="font-size:11px;background:#e6f5d8;color:#2a6010;padding:2px 8px;border-radius:4px;font-weight:500;">🐭 In vivo confirmed</span>'
    elif iv_conf == "likely":
        iv_badge = '<span style="font-size:11px;background:#fff3d0;color:#7a5000;padding:2px 8px;border-radius:4px;">🐭 In vivo likely</span>'
    else:
        iv_badge = '<span style="font-size:11px;background:#f5f5f5;color:#aaa;padding:2px 8px;border-radius:4px;">⚗ In vitro / unclear</span>'

    # Key findings as a small grid table
    findings = p.get("key_findings", [])
    if findings:
        rows = "".join(
            f'<tr>'
            f'<td style="font-size:11px;color:#888;padding:2px 16px 2px 0;white-space:nowrap;">{f["key"]}</td>'
            f'<td style="font-size:12px;font-weight:500;color:#1a1a18;padding:2px 0;">{f["val"]}</td>'
            f'</tr>'
            for f in findings
        )
        findings_html = f'<table style="border-collapse:collapse;margin-bottom:10px;">{rows}</table>'
    else:
        findings_html = '<div style="font-size:12px;color:#bbb;margin-bottom:10px;font-style:italic;">No structured data extracted from abstract</div>'


    pubmed_link = (
        f'<a href="https://pubmed.ncbi.nlm.nih.gov/{p["pmid"]}/" target="_blank" '
        f'style="font-size:12px;color:#1a1a18;border:1px solid #ddd;border-radius:5px;'
        f'padding:3px 9px;text-decoration:none;">Full record →</a>'
    ) if p["pmid"] else ""

    st.markdown(f"""
    <div class="paper-card">
      <div style="display:flex;justify-content:space-between;align-items:flex-start;gap:10px;margin-bottom:8px;">
        <div style="font-size:15px;font-weight:600;line-height:1.4;flex:1;color:#1a1a18;">{p["title"]}</div>
        <span class="cat-badge {cat_class}" style="flex-shrink:0;">{cat_label}</span>
      </div>
      <div style="display:flex;flex-wrap:wrap;gap:12px;margin-bottom:12px;">
        <span style="font-size:12px;color:#5a5a56;">
          👤 {p["authors"] or "No author info"}
        </span>
        <span style="font-size:12px;color:#5a5a56;">
          📅 {p["year"] or "—"}
        </span>
        <span style="font-size:12px;color:#5a5a56;font-style:italic;">
          📖 {p["journal"] or "—"}
        </span>
      </div>
      <div style="border-top:1px solid #f0efeb;padding-top:10px;margin-bottom:10px;">
        {findings_html}
      </div>
      <div style="font-size:11px;color:#b0a090;margin-bottom:8px;font-style:italic;">
        ⚠ Extracted values are based on abstract text only. Always refer to the full paper for complete and verified data.
      </div>
      <div style="display:flex;gap:8px;align-items:center;flex-wrap:wrap;">{iv_badge}{conf_tag}{pubmed_link}</div>
    </div>
    """, unsafe_allow_html=True)


def render_coverage(papers: list[dict]):
    counts = {k: 0 for k in ["safety", "pk", "pd", "efficacy"]}
    for p in papers:
        if p["category"] in counts:
            counts[p["category"]] += 1
    cols = st.columns(4)
    labels = {"safety": "Safety / Tox", "pk": "PK / ADME", "pd": "PD / MOA", "efficacy": "In vivo Efficacy"}
    for col, (cat, label) in zip(cols, labels.items()):
        col.metric(label, counts[cat])


def render_tab(papers: list[dict], cat: str):
    filtered = papers if cat == "all" else [p for p in papers if p["category"] == cat]
    if not filtered:
        st.info("No data confirmed within the current scope of searched public literature.  \n*Note: This reflects the limits of publicly available literature and does not necessarily indicate an absence of data.*")
        return
    for p in filtered:
        render_paper(p)


# ══════════════════════════════════════════════════════════════
# FILTERS
# ══════════════════════════════════════════════════════════════
def apply_filters(papers: list[dict], year: str, species: str, keyword: str, conf: str, in_vivo_only: bool = True) -> list[dict]:
    out = papers
    if in_vivo_only:
        out = [p for p in out if p.get("in_vivo", False)]
    if year:
        out = [p for p in out if p["year"] == year]
    if species:
        out = [p for p in out if any(
            f["key"] == "Species" and species.lower() in f["val"].lower()
            for f in p.get("key_findings", [])
        )]
    if keyword:
        kw = keyword.lower()
        out = [p for p in out if kw in p["title"].lower() or kw in p["abstract"].lower()]
    if conf == "High only":
        out = [p for p in out if p["confidence"] == "high"]
    elif conf == "Needs verification":
        out = [p for p in out if p["confidence"] == "low"]
    return out


def sort_papers(papers: list[dict], sort_by: str) -> list[dict]:
    if sort_by == "Newest first":
        return sorted(papers, key=lambda p: p["year"] or "0", reverse=True)
    elif sort_by == "Oldest first":
        return sorted(papers, key=lambda p: p["year"] or "0")
    return sorted(papers, key=lambda p: p["raw_index"])


# ══════════════════════════════════════════════════════════════
# MAIN APP
# ══════════════════════════════════════════════════════════════
def main():
    # ── Header ───────────────────────────────────────────────
    st.title("🔬 Nonclinical Drug Explorer")
    st.caption("Preclinical data explorer · PubMed (relevance) + Europe PMC (newest) · ChEMBL name normalization · Keyword classification")

    st.markdown("""
<div style="background:#fffbf0;border:1px solid #f0d080;border-radius:10px;padding:14px 18px;margin-bottom:1.2rem;">
  <div style="font-size:13px;font-weight:600;color:#7a5000;margin-bottom:8px;">⚠ Important — Please read before use</div>
  <div style="font-size:12px;color:#7a5000;line-height:1.9;">
    ✅ &nbsp;Data is sourced from <strong>publicly available literature only</strong> (PubMed / Europe PMC).<br>
    ✅ &nbsp;Results may reflect <strong>publication bias</strong> — positive findings are more likely to be published.<br>
    ✅ &nbsp;<strong>GLP nonclinical study reports</strong> are confidential and are <u>not</u> included in this tool.<br>
    ✅ &nbsp;Extracted values (NOAEL, Cmax, IC50, etc.) are parsed from <strong>abstract text only</strong> and may be incomplete.<br>
    ✅ &nbsp;<strong>Always verify all data by accessing the original publication</strong> via the provided links.
  </div>
</div>
""", unsafe_allow_html=True)

    # ── Sidebar ───────────────────────────────────────────────
    with st.sidebar:
        st.header("Settings")
        max_results = st.slider("Max papers per source", 10, 40, 30, step=10)
        min_year = st.slider("Published from (year)", 2000, 2024, 2015, step=1)
        st.divider()
        st.caption("**Study type**")
        in_vivo_only = st.toggle(
            "Animal studies only",
            value=True,
            help="Show only papers with confirmed or likely in vivo / animal data. Turn off to include in vitro / cell studies."
        )
        st.divider()
        st.caption("**Search criteria**")
        st.caption(f"• PubMed: relevance-sorted, {min_year}–present")
        st.caption(f"• Europe PMC: newest-first, {min_year}–present")
        st.caption("• Merged & deduplicated by PMID")
        st.caption("• Classification: keyword rule-based")

    # ── Search bar ────────────────────────────────────────────
    col1, col2, col3 = st.columns([3, 1, 1])
    with col1:
        query = st.text_input("", placeholder="Enter drug name (code name / INN / brand name all supported)",
                              label_visibility="collapsed")
    with col2:
        mode = st.selectbox("", ["Drug name", "Company", "Indication"],
                            label_visibility="collapsed")
    with col3:
        search = st.button("Search", type="primary", use_container_width=True)

    if not query or not search:
        if "results" not in st.session_state:
            st.markdown("---")
            st.markdown("**Try searching:** `imatinib` · `semaglutide` · `osimertinib` · `non-small cell lung cancer`")
        elif "results" in st.session_state:
            pass  # fall through to render cached results
        else:
            return

    # ── Run search ────────────────────────────────────────────
    if search and query:
        st.session_state.pop("results", None)

        chembl = None
        aliases = []

        with st.status("Searching...", expanded=True) as status:
            # 1. ChEMBL
            if mode == "Drug name":
                st.write("Looking up drug info in ChEMBL...")
                chembl = resolve_chembl(query)
                if chembl:
                    aliases = [a for a in chembl.get("aliases", []) if a and a.lower() != query.lower()]
                    st.write(f"Found: **{chembl['inn']}** ({chembl['chembl_id']}) — {len(aliases)} known aliases")

            # 2. PubMed (relevance-sorted) + Europe PMC (newest-first)
            st.write(f"Searching PubMed — relevance-sorted, {min_year}–present...")
            pubmed_papers = search_pubmed(query, mode, aliases, max_results, min_year)

            st.write(f"Searching Europe PMC — newest-first, {min_year}–present...")
            epmc_papers = search_europe_pmc(query, mode, aliases, max_results, min_year)

            papers = merge_and_deduplicate(pubmed_papers, epmc_papers)

            if not papers:
                status.update(label="No results", state="error")
                st.error(f'No public nonclinical literature found for **"{query}"** within the current search scope.')
                return

            st.write(f"PubMed: **{len(pubmed_papers)}** · Europe PMC unique: **{len(epmc_papers) - (len(papers) - len(pubmed_papers))}** added → **{len(papers)}** total")

            # 3. Classify (keyword rule-based)
            papers = keyword_classify_all(papers)

            status.update(label=f"Done — {len(papers)} papers loaded", state="complete")

        st.session_state["results"]  = papers
        st.session_state["chembl"]   = chembl
        st.session_state["query"]    = query
        st.session_state["mode"]     = mode
        st.session_state["min_year"] = min_year

    # ── Render results ────────────────────────────────────────
    if "results" not in st.session_state:
        return

    papers  = st.session_state["results"]
    chembl  = st.session_state.get("chembl")
    q_label = st.session_state.get("query", query)

    st.divider()

    # Drug header
    drug_name = chembl["inn"] if chembl else q_label
    st.subheader(drug_name)

    meta_cols = st.columns([1, 1, 1, 2])
    meta_cols[0].markdown("`PubMed + Europe PMC`")
    meta_cols[1].markdown(f"`{len(papers)} papers`")
    if chembl:
        meta_cols[2].markdown(f"`{chembl['chembl_id']}`")
        if chembl.get("formula"):
            meta_cols[3].markdown(f"`{chembl['formula']} · {chembl['mw']} g/mol`")

    # Search criteria notice
    min_year_used = st.session_state.get("min_year", 2015)
    st.caption(
        f"Search criteria: PubMed (relevance-sorted, {min_year_used}–present) "
        f"+ Europe PMC (newest-first, {min_year_used}–present) · merged & deduplicated by PMID · "
        "keyword rule-based classification · publication bias toward positive results may apply"
    )

    # ChEMBL aliases
    if chembl and aliases:
        with st.expander(f"Known aliases ({len(aliases)})", expanded=False):
            st.write(" · ".join(aliases[:15]))

    # Coverage
    st.markdown("**Coverage by category**")
    render_coverage(papers)
    st.divider()

    # ── Filters ───────────────────────────────────────────────
    with st.expander("Filters & Sort", expanded=False):
        fc1, fc2, fc3, fc4, fc5 = st.columns(5)
        years   = ["All"] + sorted({p["year"] for p in papers if p["year"]}, reverse=True)
        species_vals = list({f["val"] for p in papers for f in p.get("key_findings",[]) if f["key"]=="Species"})

        f_year    = fc1.selectbox("Year",       years)
        f_species = fc2.selectbox("Species",    ["All"] + sorted(species_vals))
        f_keyword = fc3.text_input("Keyword",   placeholder="title / abstract")
        f_conf    = fc4.selectbox("Confidence", ["All", "High only", "Needs verification"])
        f_sort    = fc5.selectbox("Sort",       ["Relevance", "Newest first", "Oldest first"])

    filtered = apply_filters(
        papers,
        year       = "" if f_year    == "All" else f_year,
        species    = "" if f_species == "All" else f_species,
        keyword    = f_keyword,
        conf       = "" if f_conf    == "All" else f_conf,
        in_vivo_only = in_vivo_only,
    )
    filtered = sort_papers(filtered, f_sort)

    total_iv  = sum(1 for p in papers if p.get("in_vivo", False))
    shown_msg = f"{len(filtered)} papers shown"
    if in_vivo_only:
        shown_msg += f" (animal studies only · {total_iv}/{len(papers)} total detected)"
    else:
        shown_msg += f" (all studies · {total_iv}/{len(papers)} with animal data)"
    st.caption(shown_msg)

    # ── Category tabs ─────────────────────────────────────────
    def cat_count(cat):
        pool = filtered if cat == "all" else [p for p in filtered if p["category"] == cat]
        return len(pool)

    tabs = st.tabs([
        f"All ({cat_count('all')})",
        f"Safety / Tox ({cat_count('safety')})",
        f"PK / ADME ({cat_count('pk')})",
        f"PD / MOA ({cat_count('pd')})",
        f"In vivo Efficacy ({cat_count('efficacy')})",
    ])
    cat_keys = ["all", "safety", "pk", "pd", "efficacy"]
    for tab, cat in zip(tabs, cat_keys):
        with tab:
            render_tab(filtered, cat)

    # ── Clinical banner ───────────────────────────────────────
    st.divider()
    inn = chembl["inn"] if chembl else q_label
    chembl_id = chembl["chembl_id"] if chembl else ""
    st.info(
        f'**Link to Clinical Trials App** — Track nonclinical → clinical development for **{inn}**  \n'
        f'Link key: INN = `{inn}`' + (f'  ·  ChEMBL = `{chembl_id}`' if chembl_id else "")
    )


if __name__ == "__main__":
    main()
