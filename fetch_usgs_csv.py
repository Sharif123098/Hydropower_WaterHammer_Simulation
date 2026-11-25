"""
Download real daily river discharge from USGS and save it as a CSV.

Setup (once):
    pip install dataretrieval pandas requests

Run:
    python fetch_usgs_csv.py

Output:
    usgs_daily_discharge.csv  with columns: date, site, discharge_cfs, discharge_m3s
"""

import io
import sys

import pandas as pd

# ---------------------------------------------------------------- settings
SITES = ["01646500", "01491000"]   # Potomac River (DC), Choptank River (MD)
START = "2010-01-01"
END = "2024-12-31"
PARAM = "00060"                    # discharge, cubic feet per second
OUT_FILE = "usgs_daily_discharge.csv"
CFS_TO_M3S = 0.0283168


# ---------------------------------------------------------------- download
def fetch_with_dataretrieval(site: str) -> pd.DataFrame:
    """Preferred route: the official USGS 'dataretrieval' package."""
    from dataretrieval import nwis

    df, _meta = nwis.get_dv(sites=site, parameterCd=PARAM, start=START, end=END)
    df = df.reset_index()

    # Value column is named like '00060_Mean'; the quality-code column ends in '_cd'.
    value_cols = [c for c in df.columns if c.startswith(PARAM) and not c.endswith("_cd")]
    if not value_cols:
        raise RuntimeError("no discharge column found in response")
    date_col = "datetime" if "datetime" in df.columns else df.columns[0]

    out = pd.DataFrame(
        {
            "date": pd.to_datetime(df[date_col]).dt.tz_localize(None).dt.date,
            "site": site,
            "discharge_cfs": pd.to_numeric(df[value_cols[0]], errors="coerce"),
        }
    )
    return out


def fetch_with_requests(site: str) -> pd.DataFrame:
    """Fallback route: the legacy USGS daily-values web service (RDB text)."""
    import requests

    url = "https://waterservices.usgs.gov/nwis/dv/"
    params = {
        "format": "rdb",
        "sites": site,
        "parameterCd": PARAM,
        "startDT": START,
        "endDT": END,
    }
    r = requests.get(url, params=params, timeout=60)
    r.raise_for_status()

    # RDB: '#' comment lines, then a header row, then a row of column widths.
    lines = [ln for ln in r.text.splitlines() if not ln.startswith("#")]
    df = pd.read_csv(io.StringIO("\n".join(lines)), sep="\t", skiprows=[1])

    value_cols = [c for c in df.columns if PARAM in c and not c.endswith("_cd")]
    if not value_cols:
        raise RuntimeError("no discharge column found in response")

    out = pd.DataFrame(
        {
            "date": pd.to_datetime(df["datetime"]).dt.date,
            "site": site,
            "discharge_cfs": pd.to_numeric(df[value_cols[0]], errors="coerce"),
        }
    )
    return out


def fetch_site(site: str) -> pd.DataFrame:
    try:
        return fetch_with_dataretrieval(site)
    except Exception as e:  # noqa: BLE001
        print(f"[{site}] dataretrieval failed ({e}); trying legacy web service...")
        return fetch_with_requests(site)


# ---------------------------------------------------------------- main
def main() -> None:
    frames = []
    for site in SITES:
        print(f"Downloading site {site} ...")
        frames.append(fetch_site(site))

    data = pd.concat(frames, ignore_index=True)
    data = data.dropna(subset=["discharge_cfs"]).sort_values(["site", "date"])
    data["discharge_m3s"] = (data["discharge_cfs"] * CFS_TO_M3S).round(4)

    data.to_csv(OUT_FILE, index=False)
    print(f"\nSaved {len(data)} rows to {OUT_FILE}")
    print(data.head())


if __name__ == "__main__":
    try:
        main()
    except Exception as e:  # noqa: BLE001
        print(f"Failed: {e}", file=sys.stderr)
        sys.exit(1)
