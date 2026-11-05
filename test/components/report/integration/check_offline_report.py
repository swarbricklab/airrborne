#!/usr/bin/env python3
"""Open a copied release offline and compare plotted recovery with its table."""

import argparse
import csv
import json
from pathlib import Path
import shutil
import tempfile
from urllib.parse import unquote, urlsplit

from playwright.sync_api import sync_playwright


def settled(page) -> None:
    page.evaluate("async () => { await new Promise(requestAnimationFrame); await trust4_report.ready; }")
    assert page.locator("#report-error").inner_text() == ""


def check(release: Path) -> dict:
    """Check a recipient's copied files without relying on the source run tree."""
    with tempfile.TemporaryDirectory(prefix="airrborne-offline-") as temporary:
        copied = (Path(temporary) / "release").resolve()
        shutil.copytree(release, copied)
        with (copied / "metrics/headline.tsv").open(newline="") as handle:
            headline = {row["library_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
        assert headline, "The copied release has no headline libraries."
        errors, network = [], []
        with sync_playwright() as playwright:
            browser = playwright.chromium.launch()
            page = browser.new_page(offline=True)
            page.on("pageerror", lambda error: errors.append(str(error)))
            page.on("request", lambda request: network.append(request.url)
                    if urlsplit(request.url).scheme in {"http", "https"} else None)
            page.goto((copied / "report/index.html").as_uri())
            page.wait_for_function("window.trust4_report && trust4_report.render_count > 0")
            settled(page)
            observed_ids = page.evaluate("[...trust4_report.controller.libraries.keys()].sort()")
            assert observed_ids == sorted(headline), "Report libraries differ from the delivered headline."
            for figure, loci in (("signal-bcr", ("IGH", "IGK", "IGL")),
                                 ("signal-tcr", ("TRA", "TRB", "TRD", "TRG"))):
                observed = page.evaluate("""id => [...document.querySelectorAll(`#${id} .js-plotly-plot`)]
                    .flatMap(div => div.data.filter(trace => trace.meta?.kind === 'library')
                    .flatMap(trace => trace.customdata.map(value => [value[0], value[4], value[2]])))""", figure)
                expected = {(library, locus): float(row[f"{locus.lower()}_yield_per_million_read_pairs"])
                            for library, row in headline.items() for locus in loci
                            if row[f"{locus.lower()}_yield_per_million_read_pairs"] not in {"", "NA"}
                            and float(row[f"{locus.lower()}_yield_per_million_read_pairs"]) > 0}
                assert {(library, locus): value for library, locus, value in observed} == expected, figure

            # Follow one chapter link to exercise the rendered report's controls.
            page.locator('a[data-chapter="diversity"]').click()
            settled(page)
            assert page.locator("#diversity-section").is_visible()
            links = page.locator("a[href], link[href], script[src]")
            for href in links.evaluate_all("nodes => nodes.map(node => node.href || node.src)"):
                parsed = urlsplit(href)
                assert parsed.scheme == "file", f"Report link needs the network: {href}"
                destination = Path(unquote(parsed.path)).resolve()
                assert destination.is_relative_to(copied) and destination.is_file(), href
            assert not errors, errors
            assert not network, network
            browser.close()
        return {"scope": next(iter(headline.values()))["scope_id"],
                "libraries": len(headline), "browser": "chromium",
                "offline": True, "recovery_values_match": True}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--release", type=Path, required=True)
    arguments = parser.parse_args()
    print(json.dumps(check(arguments.release.resolve()), indent=2))
