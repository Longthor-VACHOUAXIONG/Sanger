"""
Minimal i18n helper for the Sanger Streamlit app (English + Lao).

Font handling:
  * If fonts/NotoSansLao-Regular.ttf (and optionally -Bold.ttf) exist next to
    this file, they are embedded as base64 @font-face (works offline).
  * Otherwise Noto Sans Lao is loaded from Google Fonts.
"""

import base64
from pathlib import Path

import streamlit as st

_LANG = "en"
FONT_DIR = Path(__file__).parent / "fonts"

LANGUAGES = {"en": "English", "lo": "ລາວ"}

TRANSLATIONS = {
    "en": {
        "language": "Language",
        "sidebar_pipeline": "**Pipeline:** Sanger Sequencing Data Analysis",
        "steps_title": "**Steps:**",
        "steps": (
            "1. AB1 → FASTQ conversion\n"
            "2. Quality trimming (seqtk)\n"
            "3. Reverse complement\n"
            "4. Paired-end merge\n"
            "5. MAFFT alignment\n"
            "6. Consensus generation\n"
            "7. Quality report"
        ),
        "settings": "Settings",
        "trim_label": "Trim quality threshold",
        "trim_help": "Lower = more aggressive trimming. Default 0.05 ≈ Phred 13.",
        "links": "Links",
        "usage": "📖 Usage Guide",
        "cleanup_hint": "🔒 Use 'Clean up temp files' after downloading results.",
        "title": "🧬 Sanger Sequencing",
        "intro": "Upload your AB1 files and run the Sanger pipeline.",
        "mode_label": "Select mode:",
        "mode_single": "Single pair",
        "mode_multi": "Multiple samples (auto-detect pairs)",
        "upload_header": "### 📁 Upload AB1 Files",
        "fwd_label": "Forward reads (.ab1)",
        "rev_label": "Reverse reads (.ab1)",
        "mixed_label": "Upload ALL AB1 files (forward + reverse, mixed)",
        "mixed_help": "Pairs are auto-detected by filename (e.g., *_forward.ab1 / *_reverse.ab1)",
        "info_single": "⬆️ Upload forward and reverse AB1 files to start.",
        "info_multi": "⬆️ Upload AB1 files to start.",
        "run": "🚀 Run Workflow",
        "saved_pair": "📁 Saved {nf} forward + {nr} reverse files",
        "saved_n": "📁 Saved {n} files",
        "running": "🚀 Running workflow...",
        "setup": "Setting up pipeline...",
        "starting": "Starting...",
        "done": "✅ Workflow completed in {s}s",
        "failed": "❌ Workflow failed after {s}s",
        "batch_running": "🚀 Running batch processing...",
        "detect_pairs": "Auto-detecting pairs...",
        "batch_done": "✅ Batch completed in {s}s",
        "batch_failed": "❌ Batch failed after {s}s",
        "final_ok": "✅ Workflow completed in {s}s!",
        "final_fail": "❌ Workflow failed!",
        "full_output": "📋 Full workflow output",
        "results": "## 📊 Results",
        "m_consensus": "Consensus",
        "m_gc": "GC Content",
        "m_amb": "Ambiguity",
        "m_cov": "Coverage",
        "m_rating": "Rating",
        "detailed": "Detailed Statistics",
        "fwd_file": "Forward file",
        "rev_file": "Reverse file",
        "raw_len": "Raw length (F/R)",
        "trim_len": "Trimmed length (F/R)",
        "aln_len": "Alignment length",
        "valid_bases": "Valid bases",
        "amb_bases": "Ambiguous bases",
        "n_bases": "N bases",
        "fwd_q": "Forward avg quality",
        "rev_q": "Reverse avg quality",
        "consensus_seq": "🧬 Consensus Sequence",
        "download": "### ⬇️ Download Results",
        "samples": "Samples",
        "avg_consensus": "Avg Consensus",
        "avg_gc": "Avg GC",
        "avg_amb": "Avg Ambiguity",
        "avg_cov": "Avg Coverage",
        "quality_dist": "### Quality Distribution",
        "sample_results": "### Sample Results",
        "col_sample": "Sample",
        "col_consensus": "Consensus (bp)",
        "col_gc": "GC %",
        "col_amb": "Ambiguity %",
        "col_cov": "Coverage %",
        "col_fwdqv": "Fwd QV",
        "col_revqv": "Rev QV",
        "col_rating": "Rating",
        "individual": "### 🧬 Individual Consensus Sequences",
        "cleanup": "🗑️ Clean up temp files",
        "cleaned": "Cleaned up!",
        "footer": "Sanger Sequencing Data Analysis",
        "footer_help": "Upload AB1 files, run the pipeline, and download results.",
        "Excellent": "Excellent",
        "Good": "Good",
        "Fair": "Fair",
        "Poor": "Poor",
    },
    "lo": {
        "language": "ພາສາ",
        "sidebar_pipeline": "**Pipeline:** ການວິເຄາະຂໍ້ມູນ Sanger Sequencing",
        "steps_title": "**ຂັ້ນຕອນ:**",
        "steps": (
            "1. ແປງ AB1 → FASTQ\n"
            "2. ຕັດສ່ວນທີ່ມີຄຸນນະພາບຕ່ຳ (seqtk)\n"
            "3. Reverse complement\n"
            "4. ລວມ paired-end\n"
            "5. ຈັດລຽງດ້ວຍ MAFFT\n"
            "6. ສ້າງ consensus\n"
            "7. ລາຍງານຄຸນນະພາບ"
        ),
        "settings": "ການຕັ້ງຄ່າ",
        "trim_label": "ເກນຄຸນນະພາບສຳລັບການຕັດ",
        "trim_help": "ຄ່າຕ່ຳ = ຕັດຫຼາຍຂຶ້ນ. ຄ່າເລີ່ມຕົ້ນ 0.05 ≈ Phred 13.",
        "links": "ລິ້ງ",
        "usage": "📖 ຄູ່ມືການໃຊ້ງານ",
        "cleanup_hint": "🔒 ກົດ 'ລຶບໄຟລ໌ຊົ່ວຄາວ' ຫຼັງຈາກດາວໂຫຼດຜົນລັບແລ້ວ.",
        "title": "🧬 Sanger Sequencing",
        "intro": "ອັບໂຫຼດໄຟລ໌ AB1 ຂອງທ່ານ ແລ້ວແລ່ນ pipeline ຂອງ Sanger.",
        "mode_label": "ເລືອກໂໝດ:",
        "mode_single": "ຄູ່ດຽວ",
        "mode_multi": "ຫຼາຍຕົວຢ່າງ (ກວດຈັບຄູ່ອັດຕະໂນມັດ)",
        "upload_header": "### 📁 ອັບໂຫຼດໄຟລ໌ AB1",
        "fwd_label": "ລຳດັບ Forward (.ab1)",
        "rev_label": "ລຳດັບ Reverse (.ab1)",
        "mixed_label": "ອັບໂຫຼດໄຟລ໌ AB1 ທັງໝົດ (forward + reverse ປົນກັນ)",
        "mixed_help": "ກວດຈັບຄູ່ຈາກຊື່ໄຟລ໌ (ເຊັ່ນ *_forward.ab1 / *_reverse.ab1)",
        "info_single": "⬆️ ອັບໂຫຼດໄຟລ໌ AB1 forward ແລະ reverse ເພື່ອເລີ່ມຕົ້ນ.",
        "info_multi": "⬆️ ອັບໂຫຼດໄຟລ໌ AB1 ເພື່ອເລີ່ມຕົ້ນ.",
        "run": "🚀 ແລ່ນ Workflow",
        "saved_pair": "📁 ບັນທຶກແລ້ວ {nf} ໄຟລ໌ forward + {nr} ໄຟລ໌ reverse",
        "saved_n": "📁 ບັນທຶກແລ້ວ {n} ໄຟລ໌",
        "running": "🚀 ກຳລັງແລ່ນ workflow...",
        "setup": "ກຳລັງຕັ້ງຄ່າ pipeline...",
        "starting": "ກຳລັງເລີ່ມ...",
        "done": "✅ Workflow ສຳເລັດໃນ {s} ວິນາທີ",
        "failed": "❌ Workflow ລົ້ມເຫຼວຫຼັງຈາກ {s} ວິນາທີ",
        "batch_running": "🚀 ກຳລັງປະມວນຜົນເປັນຊຸດ...",
        "detect_pairs": "ກຳລັງກວດຈັບຄູ່ອັດຕະໂນມັດ...",
        "batch_done": "✅ ປະມວນຜົນເປັນຊຸດສຳເລັດໃນ {s} ວິນາທີ",
        "batch_failed": "❌ ປະມວນຜົນເປັນຊຸດລົ້ມເຫຼວຫຼັງຈາກ {s} ວິນາທີ",
        "final_ok": "✅ Workflow ສຳເລັດໃນ {s} ວິນາທີ!",
        "final_fail": "❌ Workflow ລົ້ມເຫຼວ!",
        "full_output": "📋 ຜົນລັບ workflow ທັງໝົດ",
        "results": "## 📊 ຜົນລັບ",
        "m_consensus": "Consensus",
        "m_gc": "ປະລິມານ GC",
        "m_amb": "ຄວາມບໍ່ແນ່ນອນ",
        "m_cov": "ການຄຸມ (Coverage)",
        "m_rating": "ລະດັບ",
        "detailed": "ສະຖິຕິລະອຽດ",
        "fwd_file": "ໄຟລ໌ Forward",
        "rev_file": "ໄຟລ໌ Reverse",
        "raw_len": "ຄວາມຍາວດິບ (F/R)",
        "trim_len": "ຄວາມຍາວຫຼັງຕັດ (F/R)",
        "aln_len": "ຄວາມຍາວການຈັດລຽງ",
        "valid_bases": "ເບສທີ່ຖືກຕ້ອງ",
        "amb_bases": "ເບສທີ່ບໍ່ແນ່ນອນ",
        "n_bases": "ເບສ N",
        "fwd_q": "ຄຸນນະພາບສະເລ່ຍ Forward",
        "rev_q": "ຄຸນນະພາບສະເລ່ຍ Reverse",
        "consensus_seq": "🧬 ລຳດັບ Consensus",
        "download": "### ⬇️ ດາວໂຫຼດຜົນລັບ",
        "samples": "ຈຳນວນຕົວຢ່າງ",
        "avg_consensus": "Consensus ສະເລ່ຍ",
        "avg_gc": "GC ສະເລ່ຍ",
        "avg_amb": "ຄວາມບໍ່ແນ່ນອນສະເລ່ຍ",
        "avg_cov": "Coverage ສະເລ່ຍ",
        "quality_dist": "### ການແຈກຢາຍຄຸນນະພາບ",
        "sample_results": "### ຜົນລັບຂອງແຕ່ລະຕົວຢ່າງ",
        "col_sample": "ຕົວຢ່າງ",
        "col_consensus": "Consensus (bp)",
        "col_gc": "GC %",
        "col_amb": "ຄວາມບໍ່ແນ່ນອນ %",
        "col_cov": "Coverage %",
        "col_fwdqv": "Fwd QV",
        "col_revqv": "Rev QV",
        "col_rating": "ລະດັບ",
        "individual": "### 🧬 ລຳດັບ Consensus ແຕ່ລະຕົວຢ່າງ",
        "cleanup": "🗑️ ລຶບໄຟລ໌ຊົ່ວຄາວ",
        "cleaned": "ລຶບສຳເລັດແລ້ວ!",
        "footer": "ການວິເຄາະຂໍ້ມູນ Sanger Sequencing",
        "footer_help": "ອັບໂຫຼດໄຟລ໌ AB1, ແລ່ນ pipeline, ແລ້ວດາວໂຫຼດຜົນລັບ.",
        "Excellent": "ດີເລີດ",
        "Good": "ດີ",
        "Fair": "ພໍໃຊ້",
        "Poor": "ບໍ່ດີ",
    },
}


def set_lang(lang: str) -> None:
    global _LANG
    _LANG = lang if lang in TRANSLATIONS else "en"


def get_lang() -> str:
    return _LANG


def t(key: str, **kwargs) -> str:
    """Translate a key (falls back to English, then to the key itself)."""
    text = TRANSLATIONS.get(_LANG, {}).get(key) or TRANSLATIONS["en"].get(key, key)
    return text.format(**kwargs) if kwargs else text


def _font_face_css() -> str:
    """Embed local font files if present, else use Google Fonts."""
    regular = FONT_DIR / "NotoSansLao-Regular.ttf"
    bold = FONT_DIR / "NotoSansLao-Bold.ttf"
    if regular.exists():
        def b64(p):
            return base64.b64encode(p.read_bytes()).decode()

        css = (
            "@font-face{font-family:'Noto Sans Lao';font-weight:400;"
            f"src:url(data:font/ttf;base64,{b64(regular)}) format('truetype');}}"
        )
        if bold.exists():
            css += (
                "@font-face{font-family:'Noto Sans Lao';font-weight:700;"
                f"src:url(data:font/ttf;base64,{b64(bold)}) format('truetype');}}"
            )
        return css
    return (
        "@import url('https://fonts.googleapis.com/css2?"
        "family=Noto+Sans+Lao:wght@400;600;700&display=swap');"
    )


def inject_font() -> None:
    """Apply Noto Sans Lao to the UI.

    Selectors are deliberately scoped to text elements so Streamlit's icon
    font and monospace code blocks are not overridden. The font is listed
    first, so Latin text also uses it, with the original sans-serif as fallback.
    """
    st.markdown(
        f"""
<style>
{_font_face_css()}
.stApp, .stApp p, .stApp label, .stApp li, .stApp a,
.stApp h1, .stApp h2, .stApp h3, .stApp h4,
.stApp button, .stApp input, .stApp textarea,
.stApp [data-testid="stMarkdownContainer"],
.stApp [data-testid="stMetricLabel"],
.stApp [data-testid="stMetricValue"],
.stApp [data-testid="stExpander"] summary,
.stApp [data-testid="stFileUploaderDropzoneInstructions"],
.stApp [data-testid="stSidebar"] * :not([data-testid="stIconMaterial"]) {{
    font-family: 'Noto Sans Lao', 'Source Sans', sans-serif;
}}
.stApp [data-testid="stIconMaterial"] {{
    font-family: 'Material Symbols Rounded', 'Material Icons' !important;
}}
.stApp code, .stApp pre {{
    font-family: 'Source Code Pro', monospace;
}}
</style>
""",
        unsafe_allow_html=True,
    )
