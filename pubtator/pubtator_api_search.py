#!/usr/bin/env python3
"""Count PubTator3 publications mentioning each taxon together with a condition term,
via the PubTator3 search API. Requests are rate limited (default 2.5/sec, limit is 3/sec).
Resumes from an existing --output file.

Usage:
    python3 pubtator_api_search.py --input taxa_union_wilcoxon.csv \
        --condition "ASD" --taxa-col cleaned_taxa --output pubtator_counts.csv
"""
import argparse
import os
import sys
import time

import pandas as pd
import requests

API_URL = "https://www.ncbi.nlm.nih.gov/research/pubtator3-api/search/"
DEFAULT_RATE = 2.5  # requests/second, safely under PubTator3's 3 req/sec cap
MAX_RETRIES = 3


def query_count(term_a, term_b, session, timeout=15):
    """Returns the total PubTator3 publication count for `"term_a" AND "term_b"` (or `"term_a"` alone if term_b is empty)."""
    query = f'"{term_a}" AND "{term_b}"' if term_b else f'"{term_a}"'
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            resp = session.get(API_URL, params={"text": query}, timeout=timeout)
            resp.raise_for_status()
            return resp.json().get("count", 0)
        except (requests.RequestException, ValueError) as e:
            if attempt == MAX_RETRIES:
                print(f"  WARNING: giving up on '{query}' after {MAX_RETRIES} attempts ({e})", file=sys.stderr)
                return None
            time.sleep(2 * attempt)
    return None


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", required=True, help="taxa_union_wilcoxon.csv-style CSV (taxa as index)")
    parser.add_argument("--taxa-col", default=None,
                         help="column holding the cleaned search term per taxon; if omitted, uses the index itself")
    parser.add_argument("--condition", required=True, help='disease/condition term to AND with each taxon, e.g. "ASD"; "" = taxon alone (popularity)')
    parser.add_argument("--output", required=True, help="output CSV path")
    parser.add_argument("--rate", type=float, default=DEFAULT_RATE, help=f"requests/second (default {DEFAULT_RATE}, cap 3)")
    args = parser.parse_args()

    if args.rate > 3:
        print("WARNING: PubTator3's own limit is 3 req/sec -- clamping rate to 3", file=sys.stderr)
        args.rate = 3.0
    delay = 1.0 / args.rate

    df = pd.read_csv(args.input, index_col=0)
    search_terms = df[args.taxa_col] if args.taxa_col else pd.Series(df.index, index=df.index)

    if os.path.exists(args.output):
        existing = pd.read_csv(args.output, index_col=0)
        done = set(existing.index)
        print(f"resuming: {len(done)} taxa already have results in {args.output}")
    else:
        existing = pd.DataFrame(columns=["search_term", "condition", "pubtator_count"])
        existing.index.name = df.index.name or "taxon"
        done = set()

    session = requests.Session()
    session.headers.update({"User-Agent": "mic_bc_benchmark-pubtator-prep/1.0 (research use)"})

    rows = []
    total = len(search_terms)
    for i, (taxon, term) in enumerate(search_terms.items(), 1):
        if taxon in done:
            continue
        count = query_count(term, args.condition, session)
        rows.append({"taxon": taxon, "search_term": term, "condition": args.condition, "pubtator_count": count})
        print(f"[{i}/{total}] {term!r} AND {args.condition!r} -> {count}")
        time.sleep(delay)

        # checkpoint every 25 queries in case of interruption
        if len(rows) % 25 == 0:
            pd.concat([existing, pd.DataFrame(rows).set_index("taxon")]).to_csv(args.output)

    if rows:
        new_df = pd.DataFrame(rows).set_index("taxon")
        pd.concat([existing, new_df]).to_csv(args.output)
    print(f"DONE -> {args.output}")


if __name__ == "__main__":
    main()
