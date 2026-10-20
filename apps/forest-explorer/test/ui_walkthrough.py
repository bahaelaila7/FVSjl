#!/usr/bin/env python3
"""Headless UI walk-through of the Pan-enabled Forest Explorer (Playwright + Chromium); screenshots, checks, console errors.
   LD_LIBRARY_PATH=<chromium libs> PLAYWRIGHT_BROWSERS_PATH=<pw> FONTCONFIG_FILE=<fonts.conf> .uivenv/bin/python test/ui_walkthrough.py [base_url] [out_dir]
(on this machine: LD_LIBRARY_PATH=/workspace/fe_pan/.uitools/env/lib PLAYWRIGHT_BROWSERS_PATH=/workspace/fe_pan/.uitools/pw
 FONTCONFIG_FILE=/workspace/fe_pan/.uitools/env/etc/fonts/fonts.conf /workspace/fe_pan/.uivenv/bin/python test/ui_walkthrough.py http://127.0.0.1:8089 /tmp/fe_ui)
Checks (PASS/FAIL lines, exit code 1 on any FAIL):
  initial state (empty panels hidden, accordion defaults, keyboard toggle) · engine list incl. disabled HBAD3 placeholder and the engine
  descriptions (LANDIS-II attribution) · AOI · one scenario per engine (Biomass thin+plant, PnET grow, HBAD v1 harvest, HBAD v2 grow / thin /
  plant, FVSjl with a PPE thin) · PPE budget inside the plan panel · compare + downloaded report with a valid, non-blank embedded AOI map PNG
  and the engine descriptions · bundle sample · phone width · no console errors."""
import sys, os, time, json, re, base64
from playwright.sync_api import sync_playwright

BASE = sys.argv[1] if len(sys.argv) > 1 else "http://127.0.0.1:8088"
OUT = sys.argv[2] if len(sys.argv) > 2 else "/tmp/fe_ui"
os.makedirs(OUT, exist_ok=True)
ok = True

def check(name, cond, detail=""):
    global ok
    print(("PASS " if cond else "FAIL ") + name + (f" — {detail}" if detail else ""), flush=True)
    ok &= bool(cond)

# small AOI (~1.2 km box inside the Oconee sample): fast runs
RING = [[-83.405, 33.410], [-83.392, 33.410], [-83.392, 33.421], [-83.405, 33.421], [-83.405, 33.410]]
AOI = json.dumps({"type": "FeatureCollection", "features": [{"type": "Feature", "properties": {}, "geometry": {"type": "Polygon", "coordinates": [RING]}}]}).encode()

PNG_VARIANCE_JS = """async (u) => { const im = new Image(); im.src = u; await im.decode();
 const c = document.createElement('canvas'); c.width = im.width; c.height = im.height; const x = c.getContext('2d'); x.drawImage(im, 0, 0);
 const d = x.getImageData(0, 0, c.width, c.height).data; let n = 0, s = 0, s2 = 0; const set = new Set();
 for (let i = 0; i < d.length; i += 28) { const v = (d[i] + d[i + 1] + d[i + 2]) / 3; n++; s += v; s2 += v * v; set.add((d[i] >> 4) + ',' + (d[i + 1] >> 4) + ',' + (d[i + 2] >> 4)); }
 const m = s / n; return { w: im.width, h: im.height, variance: s2 / n - m * m, colors: set.size }; }"""

def shot(pg, name):
    pg.screenshot(path=os.path.join(OUT, name + ".png"), full_page=False)

with sync_playwright() as p:
    b = p.chromium.launch(args=["--use-gl=angle", "--use-angle=swiftshader", "--enable-unsafe-swiftshader"])
    ctx = b.new_context(viewport={"width": 1440, "height": 900}, accept_downloads=True)
    pg = ctx.new_page()
    errs = []
    pg.on("console", lambda m: errs.append(f"{m.type}: {m.text}") if m.type == "error" else None)
    pg.on("pageerror", lambda e: errs.append(f"pageerror: {e}"))
    pg.on("dialog", lambda d: d.accept())
    pg.goto(BASE + "/", wait_until="networkidle", timeout=120000)
    time.sleep(3); shot(pg, "01_initial")

    # ---- 1. initial state: nothing preallocated for empty content, accordion defaults
    vis = lambda sel: pg.is_visible(sel)
    check("initial: projection column, summary, chart and sim-layer groups hidden", not any(vis(s) for s in ["#planR", "#statsGrp", "#chartGrp", "#simLayerGrp", "#landscapeGrp"]))
    op = pg.evaluate("Object.fromEntries(['accSource','accAoi','accLayers','bundleGrp','planDet'].map(i => [i, document.getElementById(i).open]))")
    check("initial: accordion defaults (source+AOI open, layers+bundle+plan closed)", op == {"accSource": True, "accAoi": True, "accLayers": False, "bundleGrp": False, "planDet": False}, str(op))
    pg.focus("#accLayers > summary"); pg.keyboard.press("Enter")
    check("accordion toggles from the keyboard (Enter on the summary)", pg.evaluate("document.getElementById('accLayers').open") is True)
    pg.keyboard.press("Enter")
    check("PPE budget is no longer in the sidebar", pg.evaluate("!document.querySelector('aside #landscapeGrp')"))

    # ---- 2. engines
    opts = pg.eval_on_selector_all("#engine option", "els => els.map(e => [e.value, e.textContent, e.disabled])")
    ids = [o[0] for o in opts]
    check("engine list: FVSjl, Biomass, PnET, HBAD v1, HBAD v2, HBAD3 placeholder", ids == ["fvs", "pan-biomass", "pan-pnet", "pan-hbad", "pan-hbad2", "pan-hbad3"], str(opts))
    h3 = [o for o in opts if o[0] == "pan-hbad3"][0]
    check("HBAD3 placeholder is disabled and says 'coming'", h3[2] and "coming" in h3[1] and "HBAD3 (process-based)" in h3[1], str(h3))
    want = {"fvs": "Dixon 2002", "pan-biomass": "Scheller et al. 2007", "pan-pnet": "de Bruijn et al. 2014", "pan-hbad": "Reineke", "pan-hbad2": "overtopping"}
    # the engine text lives in a <details>; open it to check visibility
    pg.evaluate("document.getElementById('planDet').open = true")
    for eid, key in want.items():
        pg.select_option("#engine", eid); pg.wait_for_timeout(300)
        pg.evaluate("document.getElementById('engineHelp').open = true")
        t = pg.inner_text("#engineText")
        check(f"engine description visible: {eid}", pg.is_visible("#engineText") and key in t, t[:90])
    pg.select_option("#engine", "fvs")
    pg.evaluate("document.getElementById('planDet').open = false")

    # ---- 3. AOI
    pg.set_input_files("#file", files=[{"name": "small.geojson", "mimeType": "application/geo+json", "buffer": AOI}])
    pg.wait_for_function("document.getElementById('status').textContent.includes('plots') && !document.getElementById('status').textContent.includes('resolving')", timeout=120000)
    pg.wait_for_timeout(4000); shot(pg, "02_aoi")
    check("AOI: plan panel opened, summary + chart groups shown, projection column still hidden",
          pg.evaluate("document.getElementById('planDet').open") and vis("#statsGrp") and vis("#chartGrp") and not vis("#planR"), pg.inner_text("#status")[:80])
    check("PPE budget group appears inside the plan panel (FVSjl, AOI set)", pg.evaluate("!!document.querySelector('#planpanel #landscapeGrp')") and vis("#landscapeGrp"))

    def run(name, engine, actions=(), cycles=4, params=None):
        """new scenario → engine → actions via the UI buttons → run → wait for its result"""
        pg.click("#newScenario")
        pg.evaluate("plan.actions = []; renderActions();")
        pg.fill("#scenName", name)
        pg.select_option("#engine", engine); pg.wait_for_timeout(1500)
        if params:
            pg.select_option("#panSamples", params); pg.wait_for_timeout(8000)
        for a in actions:
            pg.click({"thin": "#addThin", "plant": "#addPlant", "harvest": "#addHarvest"}[a])
        pg.fill("#pCycles", str(cycles)); pg.fill("#pPeriod", "10"); pg.dispatch_event("#pCycles", "change"); pg.dispatch_event("#pPeriod", "change")
        return actions

    def go(name, engine, actions=(), cycles=4, params=None, mutate=None):
        run(name, engine, actions, cycles, params)
        if mutate: pg.evaluate(mutate);
        pg.click("#runSim")
        t0 = time.time()
        pg.wait_for_function(f"scenarios.some(s => s.name === {json.dumps(name)} && s.result) || /error/i.test(document.getElementById('status').textContent)", timeout=1800000)
        sc = pg.evaluate(f"(() => {{ const s = scenarios.find(s => s.name === {json.dumps(name)}); return s && s.result ? {{ n: s.result.ncycles, nplots: s.result.nplots_sim, eng: s.result.engine, carbon: s.result.cycles.map(c => c.carbon_live_ag), removed: s.result.cycles.map(c => c.removed_volume), ppe: !!s.result.ppe }} : null; }})()")
        good = bool(sc) and sc["n"] == cycles + 1 and sc["nplots"] > 0 and all(c > 0 for c in sc["carbon"][:1])
        check(f"scenario '{name}' ({engine}) ran", good, f"{sc} in {round(time.time()-t0)}s" if sc else pg.inner_text("#status")[:200])
        return sc

    # ---- 4. each engine runs a scenario
    sc = go("BS thin+plant", "pan-biomass", ["thin", "plant"])
    check("Biomass thin removed volume", sc and sum(sc["removed"]) > 0)
    go("PnET grow", "pan-pnet")
    sc = go("HBAD1 harvest", "pan-hbad", ["harvest"])
    check("HBAD v1 harvest removed volume", sc and sum(sc["removed"]) > 0)
    sc = go("HBAD2 grow", "pan-hbad2"); check("HBAD v2 engine id in result", sc and sc["eng"] == "pan-hbad2")
    check("HBAD v2 parameter set is the default v2j", "v2j" in pg.inner_text("#panParamsInfo"), pg.inner_text("#panParamsInfo")[:120])
    sc = go("HBAD2 thin", "pan-hbad2", ["harvest"], mutate="plan.actions[0].removal = {type:'residual_ba', value:60, direction:'below'}; plan.actions[0].cycle = 1; renderActions();")
    check("HBAD v2 thin removed volume", sc and sum(sc["removed"]) > 0)
    sc = go("HBAD2 plant", "pan-hbad2", ["plant"], mutate="plan.actions[0].species = 'PITA'; renderActions();")
    check("HBAD v2 plant scenario ran", sc is not None)
    pg.wait_for_timeout(2000); shot(pg, "03_after_sim")
    check("after a simulation: projection column + sim layer + chart shown", vis("#planR") and vis("#simLayerGrp") and vis("#chartGrp") and vis("#simChart"))

    # ---- 5. PPE budget as part of the scenario (FVSjl)
    pg.click("#newScenario"); pg.evaluate("plan.actions = []; renderActions();")
    pg.select_option("#engine", "fvs"); pg.wait_for_timeout(800)
    check("PPE group visible in plan panel for FVSjl", vis("#landscapeGrp") and pg.evaluate("!!document.querySelector('#planpanel #landscapeGrp')"))
    pg.evaluate("document.getElementById('landscapeGrp').open = true")
    pg.fill("#lsTarget", "500"); pg.click("#lsAddThin"); pg.wait_for_timeout(500)
    mode = pg.evaluate("plan.actions.map(a => [a.kind, a.mode, a.target_expr])")
    check("'Add as thin action' puts the PPE budget into the plan", mode == [["thin", "ppe", "500"]], str(mode))
    sel_mode = pg.evaluate("(document.querySelector('#actions select[data-f=mode]') || {}).value")
    check("thin action shows 'PPE budget picks stands'", sel_mode == "ppe", str(sel_mode))
    sc = go("FVS PPE", "fvs", [], cycles=3, mutate="plan.actions = [{kind:'thin',cycle:1,metric:'BA',target:80,direction:'below',species:'all',dbh_lo:0,dbh_hi:999,mode:'ppe',target_expr:'500',priority_expr:'BBA',credit_expr:'BBA'}]; renderActions();")
    check("FVSjl PPE-scheduled thin: result carries the budget flow", sc and sc["ppe"], str(sc and sc["ppe"]))
    pg.click("#runLandscape"); pg.wait_for_function("document.getElementById('lsResult').innerText.length > 10 || /error/.test(document.getElementById('lsInfo').textContent)", timeout=900000)
    check("PPE preview table", "when" in pg.inner_text("#lsResult"), pg.inner_text("#lsInfo")[:80])
    shot(pg, "04_ppe_in_plan")

    # ---- 6. compare + report
    n = pg.evaluate("scenarios.filter(s => s.result).length")
    check("scenarios stored", n >= 7, f"{n} scenarios")
    pg.click("#compareBtn"); pg.wait_for_timeout(1500); shot(pg, "05_compare")
    with pg.expect_download(timeout=60000) as dl:
        pg.click("#cmpDownload")
    path = os.path.join(OUT, "report.html"); dl.value.save_as(path)
    html = open(path, encoding="utf-8").read()
    m = re.search(r'<img id="aoiimg" src="data:image/png;base64,([A-Za-z0-9+/=]+)"', html)
    check("report embeds the AOI map as a PNG data URI", bool(m), f"{len(html)} bytes")
    if m:
        raw = base64.b64decode(m.group(1))
        check("embedded image has the PNG signature", raw[:8] == b"\x89PNG\r\n\x1a\n", f"{len(raw)} bytes")
        rp = ctx.new_page(); rp.set_content(html)
        ni = rp.evaluate("(() => { const i = document.getElementById('aoiimg'); return [i.complete, i.naturalWidth, i.naturalHeight]; })()")
        check("report <img> loads (naturalWidth > 0)", ni[0] and ni[1] > 0, str(ni))
        st = rp.evaluate(PNG_VARIANCE_JS, "data:image/png;base64," + m.group(1))
        check("embedded image is not blank (decoded: pixel variance > 50, > 50 colors)", st["variance"] > 50 and st["colors"] > 50, str(st))
        rp.screenshot(path=os.path.join(OUT, "06_report.png"), full_page=False)
        for k in ["Scheller", "de Bruijn", "Dixon", "Reineke", "overtopping", "vegetation: TreeMap 2022 pixels"]:
            check(f"report contains '{k}'", k in html)
        rp.close()
    pg.keyboard.press("Escape")

    # ---- 7. phone width, bundle sample
    pg.set_viewport_size({"width": 400, "height": 860}); pg.wait_for_timeout(1500); shot(pg, "07_phone")
    pg.set_viewport_size({"width": 1440, "height": 900}); pg.wait_for_timeout(800)
    pg.evaluate("document.getElementById('bundleGrp').open = true")
    sel = pg.eval_on_selector_all("#bundleSample option", "els => els.map(e => e.value)")
    check("sample bundles listed", len(sel) > 1, str(sel))
    if len(sel) > 1:
        pg.select_option("#bundleSample", sel[1]); pg.wait_for_timeout(20000); shot(pg, "08_bundle_imported")
        check("bundle imported", "scenario" in pg.inner_text("#bundleMsg").lower() or pg.evaluate("scenarios.length") > 0, pg.inner_text("#bundleMsg")[:100])
    b.close()

# the OSM/Esri basemap tiles can 4xx/5xx or be blocked: report them but only fail on JS/app errors
app_errs = [e for e in errs if not re.search(r"Failed to load resource|tile|\.png|arcgisonline|nationalmap|openstreetmap|net::", e)]
check("no console errors", not app_errs, f"{len(errs)} console errors total, {len(app_errs)} not network; " + " || ".join(e[:160] for e in app_errs[:5]))
for e in errs[:10]: print("  console:", e[:200])
print("ALL PASS" if ok else "SOME CHECKS FAILED")
sys.exit(0 if ok else 1)
