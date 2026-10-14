#!/usr/bin/env python3
"""Headless UI walk-through of the Pan-enabled Forest Explorer (Playwright + Chromium); screenshots + console errors.
   LD_LIBRARY_PATH=<chromium libs> .uivenv/bin/python test/ui_walkthrough.py [base_url] [out_dir]
Steps: load → upload sample AOI → engine Pan·HBAD (params panel) → add harvest → run → results → engine Pan·PnET →
climate modal (preview) → phone width → sample bundle import."""
import sys, os, time, json
from playwright.sync_api import sync_playwright

BASE = sys.argv[1] if len(sys.argv) > 1 else "http://127.0.0.1:8088"
OUT = sys.argv[2] if len(sys.argv) > 2 else "/tmp/fe_ui"
HERE = os.path.dirname(os.path.abspath(__file__)); SAMPLES = os.path.join(HERE, "..", "samples")
os.makedirs(OUT, exist_ok=True)
log = []

def shot(pg, name):
    p = os.path.join(OUT, name + ".png"); pg.screenshot(path=p, full_page=False); log.append(("shot", p))

with sync_playwright() as p:
    b = p.chromium.launch(args=["--use-gl=angle", "--use-angle=swiftshader", "--enable-unsafe-swiftshader"])
    pg = b.new_page(viewport={"width": 1440, "height": 900})
    errs = []
    pg.on("console", lambda m: errs.append(f"{m.type}: {m.text}") if m.type in ("error", "warning") else None)
    pg.on("pageerror", lambda e: errs.append(f"pageerror: {e}"))
    pg.on("dialog", lambda d: d.accept())
    pg.goto(BASE + "/", wait_until="networkidle", timeout=120000)
    time.sleep(3); shot(pg, "01_initial")
    # AOI
    pg.set_input_files("#file", os.path.join(SAMPLES, "aois", "oconee_ga.geojson"))
    pg.wait_for_timeout(15000); shot(pg, "02_aoi")
    log.append(("status", pg.inner_text("#status")))
    # engine HBAD
    opts = pg.eval_on_selector_all("#engine option", "els => els.map(e => [e.value, e.textContent, e.disabled])")
    log.append(("engines", opts))
    pg.select_option("#engine", "pan-hbad"); pg.wait_for_timeout(8000); shot(pg, "03_engine_hbad_params")
    log.append(("params panel", pg.inner_text("#panParamsInfo") if pg.is_visible("#panParamsInfo") else "hidden"))
    pg.click("#addHarvest"); pg.wait_for_timeout(1500); shot(pg, "04_harvest_form")
    log.append(("actions", pg.inner_text("#actions")[:600]))
    pg.fill("#pCycles", "4"); pg.fill("#pPeriod", "10")
    pg.click("#runSim")
    t0 = time.time()
    while time.time() - t0 < 1500:
        txt = pg.inner_text("#status")
        if "error" in txt.lower() or (pg.is_visible("#simChart") and not pg.is_visible("#simEmpty")):
            break
        pg.wait_for_timeout(5000)
    pg.wait_for_timeout(4000); shot(pg, "05_results_hbad")
    log.append(("status after run", pg.inner_text("#status"), round(time.time() - t0)))
    log.append(("warnings", pg.inner_text("#simWarn") if pg.is_visible("#simWarn") else ""))
    # PnET + climate modal
    pg.click("#newScenario"); pg.wait_for_timeout(1000)
    pg.select_option("#engine", "pan-pnet"); pg.wait_for_timeout(5000)
    pg.click("#climateBtn"); pg.wait_for_timeout(60000); shot(pg, "06_climate_modal")
    pg.keyboard.press("Escape"); pg.wait_for_timeout(1000)
    # phone width
    pg.set_viewport_size({"width": 400, "height": 860}); pg.wait_for_timeout(2000); shot(pg, "07_phone")
    pg.set_viewport_size({"width": 1440, "height": 900}); pg.wait_for_timeout(1000)
    # sample bundle
    sel = pg.eval_on_selector_all("#bundleSample option", "els => els.map(e => e.value)")
    log.append(("bundle options", sel))
    if len(sel) > 1:
        pg.select_option("#bundleSample", sel[1]); pg.wait_for_timeout(30000); shot(pg, "08_bundle_imported")
        log.append(("bundle msg", pg.inner_text("#bundleMsg") if pg.is_visible("#bundleMsg") else ""))
        log.append(("chips", pg.inner_text("#scenarioChips")))
    b.close()
for e in log: print(*e)
print("CONSOLE (errors/warnings):", len(errs)); [print("  ", e[:300]) for e in errs[:40]]
