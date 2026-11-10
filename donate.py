"""
Donate section for the Sanger Streamlit app.

EDIT THE CONFIG BELOW with your own details, then put your QR images in the
same folder as this file (base project folder):
    bank_qr.png   (your bank QR, e.g. OnePay / BCEL One / other)
    usdt_qr.png   (optional - auto-generated from the address if missing)
    btc_qr.png    (optional - auto-generated from the address if missing)

Auto-generation needs:  pip install "qrcode[pil]"
"""

import io
from pathlib import Path

import streamlit as st

from i18n import t

# ---------------------------------------------------------------------------
# CONFIG - replace the placeholders
# ---------------------------------------------------------------------------
DONATE_DIR = Path(__file__).parent  # base project folder

BANK = {
    "image": DONATE_DIR / "bank_qr.png",
    "bank_name": "BCEL",          # e.g. BCEL
    "account_name": "LONGTHOR VACHOUAXIONG MR",
    "account_number": "0101210014123",
}

USDT = {
    "image": DONATE_DIR / "usdt_qr.png",
    "network": "TRC20 (Tron)",              # MUST match the address type
    "address": "TM2GbpiTN9ouNdnVnYgYViSLQtnRmPcbHh",
}

BTC = {
    "image": DONATE_DIR / "btc_qr.png",
    "network": "Bitcoin",
    "address": "1FcmyGvHmd5XKJb2mifbeS5jVWrWJbLeNX",
}

QR_WIDTH = 240


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _is_placeholder(value: str) -> bool:
    return (not value) or value.startswith("YOUR")


def _qr_image(image_path: Path, data: str | None):
    """Return a path/bytes for st.image: the file if it exists, else a generated QR."""
    if image_path.exists():
        return str(image_path)
    if data and not _is_placeholder(data):
        try:
            import qrcode

            buf = io.BytesIO()
            qrcode.make(data).save(buf, format="PNG")
            return buf.getvalue()
        except ImportError:
            return None
    return None


def _show_qr(img):
    if img is not None:
        st.image(img, width=QR_WIDTH)


def _crypto_tab(cfg: dict, coin_label: str):
    if _is_placeholder(cfg["address"]):
        st.info(t("donate_not_set"))
        return
    _show_qr(_qr_image(cfg["image"], cfg["address"]))
    st.markdown(f"**{t('donate_network')}:** {cfg['network']}")
    st.caption(t("donate_address"))
    st.code(cfg["address"], language="text")  # built-in copy button
    st.warning(t("donate_network_warning", coin=coin_label, network=cfg["network"]))


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------
def render_donate(expanded: bool = True):
    st.markdown(
        """
<style>
    div[data-testid="stExpander"]:has(.donate-marker) {
        border: 2px solid #e74c3c !important;
        border-radius: 10px;
        background: rgba(231, 76, 60, 0.05);
    }
</style>
""",
        unsafe_allow_html=True,
    )
    with st.expander(t("donate_title"), expanded=expanded):
        st.markdown('<span class="donate-marker"></span>', unsafe_allow_html=True)
        st.markdown(t("donate_intro"))
        tab_bank, tab_usdt, tab_btc = st.tabs(
            [t("donate_bank"), "USDT", "Bitcoin"]
        )

        with tab_bank:
            img = _qr_image(BANK["image"], None)
            if img is None and _is_placeholder(BANK["account_number"]):
                st.info(t("donate_not_set"))
            else:
                _show_qr(img)
                st.markdown(f"**{t('donate_bank_name')}:** {BANK['bank_name']}")
                st.markdown(f"**{t('donate_account_name')}:** {BANK['account_name']}")
                st.caption(t("donate_account_number"))
                st.code(BANK["account_number"], language="text")

        with tab_usdt:
            _crypto_tab(USDT, "USDT")

        with tab_btc:
            _crypto_tab(BTC, "Bitcoin")

        st.markdown(f"<small>{t('donate_thanks')}</small>", unsafe_allow_html=True)
