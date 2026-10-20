"""UI check of the stands layer and the "Delineated stands" vegetation source (public/stands.js).
  LD_LIBRARY_PATH=<chromium libs> PLAYWRIGHT_BROWSERS_PATH=<pw> FONTCONFIG_FILE=<fonts.conf> .uivenv/bin/python test/ui_stands.py [base_url] [out_dir]
Needs public/stands/baker_fl.geojson (git-ignored; copy it from the standmap export).
Part A (TreeMap source): stands layer toggle → popup → "Simulate this stand" sets the AOI.
Part B (stands source): vegetation source = stands → select stands by clicking and by a drawn polygon → the selection is the AOI → run Pan HBAD v2
and FVSjl → the projection map colors the STAND POLYGONS per cycle (values and rendered pixels change with the cycle) → stand popup.
Screenshots go to out_dir; PASS/FAIL lines; exit code 1 on any FAIL."""
import sys, os, time, hashlib, re
from playwright.sync_api import sync_playwright
base = sys.argv[1] if len(sys.argv) > 1 else 'http://127.0.0.1:8088'; out = sys.argv[2] if len(sys.argv) > 2 else '/tmp/ui_stands'
os.makedirs(out, exist_ok=True); ok = True

def check(name, cond, detail=""):
    global ok
    print(("PASS " if cond else "FAIL ") + name + (f" — {detail}" if detail else ""), flush=True); ok &= bool(cond)

with sync_playwright() as p:
    b = p.chromium.launch(args=["--use-gl=angle", "--use-angle=swiftshader", "--enable-unsafe-swiftshader"])
    pg = b.new_page(viewport={'width': 1500, 'height': 950})
    errs = []; pg.on('pageerror', lambda e: errs.append(str(e))); pg.on('console', lambda m: m.type == 'error' and errs.append(m.text))
    pg.on('dialog', lambda d: d.accept())
    pg.goto(base + '/', wait_until='networkidle'); time.sleep(3)

    # ---------------- Part A: TreeMap source, stands layer as before ----------------
    pg.evaluate("document.getElementById('accLayers').open = true; document.getElementById('stGrp').open = true")
    pg.evaluate("(() => { const e = document.getElementById('stOn'); e.checked = true; e.dispatchEvent(new Event('change')); })()")
    pg.wait_for_function("document.getElementById('stStatus').textContent.includes('stands')", timeout=120000); time.sleep(4)
    check('A: stands layer status', 'stands' in pg.inner_text('#stStatus'), pg.inner_text('#stStatus'))
    pg.screenshot(path=f'{out}/1_stands_status.png')
    pg.select_option('#stColor', 'BALIVE'); time.sleep(2); pg.select_option('#stColor', 'status'); time.sleep(1)
    # zoom in so single stands are clickable, then find a forest stand with TreeMap plots under the cursor
    pg.evaluate("map.jumpTo({center:[-82.286, 30.464], zoom:12.5})"); time.sleep(4)
    CAND = """() => { const c = map.getCanvas(), W = c.width / devicePixelRatio, H = c.height / devicePixelRatio, seen = new Map();
      for (let y = 60; y < H - 60; y += 14) for (let x = 60; x < W - 60; x += 14) {
        for (const f of map.queryRenderedFeatures([x, y], { layers: ['stands-fill'] })) { const pr = f.properties;
          if (pr.status === 1 && pr.tm_n > 0 && pr.area_ha > 8 && !seen.has(pr.stand_id)) seen.set(pr.stand_id, [x, y, pr.stand_id, pr.area_ha]); } }
      return [...seen.values()]; }"""
    cands = pg.evaluate(CAND)
    check('A: stands rendered and queryable', len(cands) >= 3, f'{len(cands)} candidate stands in view')
    box = pg.query_selector('canvas.maplibregl-canvas').bounding_box()
    x, y, sid, ha = cands[0]
    pg.mouse.click(box['x'] + x, box['y'] + y); time.sleep(2)
    check('A: popup (TreeMap source, click does not select)', pg.is_visible('.stpop') and pg.evaluate("aoiParts.length") == 0, pg.inner_text('.stpop')[:80].replace('\n', ' | '))
    pg.screenshot(path=f'{out}/2_popup.png')
    pg.click('#stSim'); pg.wait_for_function("document.getElementById('cards').innerText.length > 5", timeout=180000); time.sleep(3)
    check('A: "Simulate this stand" resolved the stand as AOI', pg.evaluate("aoiParts.length") == 1, pg.inner_text('#status')[:100])
    pg.screenshot(path=f'{out}/3_stand_as_aoi.png')

    # ---------------- Part B: vegetation source = delineated stands ----------------
    pg.select_option('#vegSource', 'stands'); time.sleep(3)
    opts = pg.eval_on_selector_all('#vegSource option', 'e => e.map(x => [x.value, x.textContent, x.disabled])')
    gen = [o for o in opts if o[0] == 'stands-gen'][0]
    check('B: generative-stands source is a disabled placeholder', gen[2] and 'coming' in gen[1] and 'TreeMap-independent' in gen[1], str(gen))
    check('B: switching source cleared the previous AOI, stand tools shown, patch tools hidden',
          pg.evaluate("aoiParts.length") == 0 and pg.is_visible('#standTools') and not pg.is_visible('#aoiTools'))
    check('B: stands layer switched on automatically', pg.evaluate("document.getElementById('stOn').checked"))
    pg.evaluate("map.jumpTo({center:[-82.286, 30.464], zoom:12.5})"); time.sleep(4)
    # (1) polygon selection: draw a polygon around the view centre
    pg.click('#stSelPoly'); time.sleep(0.5)
    cx, cy = box['x'] + box['width'] * 0.5, box['y'] + box['height'] * 0.5
    for dx, dy in [(-120, -90), (120, -90), (120, 90), (-120, 90)]:
        pg.mouse.click(cx + dx, cy + dy); time.sleep(0.45)
    pg.keyboard.press('Enter'); time.sleep(2.5)
    n_poly = pg.evaluate("document.getElementById('stSelCount').textContent")
    k = int(re.match(r'(\d+)', n_poly).group(1))
    check('B: "select stands inside a drawn polygon" selected stands', k >= 1, f'{n_poly}; {pg.inner_text("#status")[:90]}')
    pg.wait_for_function("aoiStandIds.length > 0 && document.getElementById('cards').innerText.length > 5", timeout=180000)
    check('B: polygon selection became the AOI (stands in the AOI list)', pg.evaluate("aoiStandIds.length") == k and 'Stand' in pg.inner_text('#patchList'))
    pg.screenshot(path=f'{out}/4_polygon_selection.png')
    pg.click('#stSelClear'); time.sleep(2)
    check('B: clear selection', pg.evaluate("aoiParts.length") == 0)
    # (2) click three stands
    cands = pg.evaluate(CAND)
    cands.sort(key=lambda c: (c[0] - 750) ** 2 + (c[1] - 475) ** 2)
    for x, y, sid, ha in cands[:3]:
        pg.mouse.click(box['x'] + x, box['y'] + y); time.sleep(0.7)
    pg.wait_for_function("aoiStandIds.length === 3 && document.getElementById('cards').innerText.length > 5", timeout=180000); time.sleep(3)
    check('B: three stands selected by clicking (toggle)', pg.evaluate("aoiStandIds.length") == 3, str(pg.evaluate("aoiStandIds")))
    # (opening the plan panel resized the map: look the stand up again before clicking it)
    sid0 = cands[0][2]
    pos = [c for c in pg.evaluate(CAND) if c[2] == sid0]
    box = pg.query_selector('canvas.maplibregl-canvas').bounding_box()
    check('B: selected stand still in view', bool(pos), f'stand {sid0}')
    x, y = pos[0][0], pos[0][1]; pg.mouse.click(box['x'] + x, box['y'] + y); time.sleep(0.8)
    check('B: clicking a selected stand again deselects it', pg.evaluate("document.getElementById('stSelCount').textContent").startswith('2'), pg.evaluate("document.getElementById('stSelCount').textContent"))
    pg.mouse.click(box['x'] + x, box['y'] + y); time.sleep(0.8)
    pg.wait_for_function("aoiStandIds.length === 3 && document.getElementById('cards').innerText.length > 5", timeout=180000); time.sleep(2)
    pg.screenshot(path=f'{out}/5_three_stands_selected.png')
    pg.evaluate("document.getElementById('planDet').open = true")

    def run_engine(name, engine, cycles=4):
        pg.click('#newScenario'); pg.evaluate("plan.actions = []; renderActions();")
        pg.fill('#scenName', name); pg.select_option('#engine', engine); pg.wait_for_timeout(1500)
        pg.fill('#pCycles', str(cycles)); pg.dispatch_event('#pCycles', 'change')
        pg.click('#runSim')
        pg.wait_for_function(f"scenarios.some(s => s.name === '{name}' && s.result) || /error/i.test(document.getElementById('status').textContent)", timeout=1800000)
        pg.wait_for_timeout(3000)
        r = pg.evaluate(f"""(() => {{ const s = scenarios.find(s => s.name === '{name}'); return s && s.result ? {{ n: s.result.ncycles, stands: (s.result.stands || []).length, nsim: s.result.nplots_sim, image: !!s.result.image,
              sims: (s.result.stands || []).map(t => t.nplots_sim), ba0: (s.result.stands || []).map(t => t.map.ba[0]), baN: (s.result.stands || []).map(t => t.map.ba[t.map.ba.length - 1]) }} : null; }})()""")
        check(f'B: {name} ({engine}) ran per stand', bool(r) and r['stands'] == 3 and r['n'] == cycles + 1 and all(s > 0 for s in r['sims']) and not r['image'], str(r) if r else pg.inner_text('#status')[:200])
        check(f'B: {name}: per-stand values differ between stands and grow over time', bool(r) and len(set(r['ba0'])) > 1 and all(b1 > a for a, b1 in zip(r['ba0'], r['baN'])), str(r and (r['ba0'], r['baN'])))
        # the map layer colors the stand polygons per cycle
        check(f'B: {name}: stand polygon layer drawn, raster not', pg.evaluate("!!map.getLayer('simstands-fill') && !map.getLayer('simlayer')"))
        clip = {'x': box['x'] + 100, 'y': box['y'] + 80, 'width': box['width'] - 200, 'height': box['height'] - 160}
        v0 = pg.evaluate("(() => { const f = map.querySourceFeatures('simstands'); return f.map(x => [x.properties.stand_id, x.properties.v]).sort(); })()")
        h0 = hashlib.md5(pg.screenshot(clip=clip)).hexdigest(); pg.screenshot(path=f'{out}/6_{engine}_cycle0.png')
        last = r['n'] - 1
        pg.evaluate(f"onSlider({last})"); pg.wait_for_timeout(2500)
        v1 = pg.evaluate("(() => { const f = map.querySourceFeatures('simstands'); return f.map(x => [x.properties.stand_id, x.properties.v]).sort(); })()")
        h1 = hashlib.md5(pg.screenshot(clip=clip)).hexdigest(); pg.screenshot(path=f'{out}/7_{engine}_cycle{last}.png')
        check(f'B: {name}: stand colors change per cycle (feature values and rendered pixels)', v0 != v1 and h0 != h1 and len({v for _, v in v0}) > 1, f'cycle0 {v0[:3]} → last {v1[:3]}')
        return r

    run_engine('Stands HBAD2', 'pan-hbad2')
    # stand popup on the projected map
    pt = pg.evaluate("""() => { const c = map.getCanvas(); for (let y = 100; y < c.height / devicePixelRatio - 100; y += 12) for (let x = 100; x < c.width / devicePixelRatio - 100; x += 12)
        if (map.queryRenderedFeatures([x, y], { layers: ['simstands-fill'] }).length) return [x, y]; return null; }""")
    pg.mouse.click(box['x'] + pt[0], box['y'] + pt[1]); time.sleep(1.5)
    check('B: projected stand popup (basal area, carbon, volume)', pg.is_visible('.stpop') and 'basal area' in pg.inner_text('.stpop'), pg.inner_text('.stpop')[:90].replace('\n', ' | '))
    check('B: sidebar legend/metric group shown for the stand layer', pg.is_visible('#simLayerGrp'))
    pg.select_option('#simMapMetric', 'ba'); pg.wait_for_timeout(1500)
    check('B: metric switch recolors stands', pg.evaluate("map.querySourceFeatures('simstands').length") >= 3)
    run_engine('Stands FVSjl', 'fvs')
    pg.screenshot(path=f'{out}/8_final.png')
    b.close()
real = [e for e in errs if not re.search(r'Failed to load resource|net::|tile', e)]
check('no console errors', not real, f'{len(errs)} total, {len(real)} non-network: ' + ' || '.join(e[:150] for e in real[:4]))
print('ALL PASS' if ok else 'SOME CHECKS FAILED'); sys.exit(0 if ok else 1)
