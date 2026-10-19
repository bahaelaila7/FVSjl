"""UI check of the stands layer (public/stands.js): toggle on → Baker County stands render → color by basal area → click a stand → popup →
"Simulate this stand" sets the AOI (the AOI pipeline resolves it). Screenshots to out_dir.
  LD_LIBRARY_PATH=<chromium libs> .uivenv/bin/python test/ui_stands.py [base_url] [out_dir]"""
import sys, os, time
from playwright.sync_api import sync_playwright
base = sys.argv[1] if len(sys.argv) > 1 else 'http://127.0.0.1:8088'; out = sys.argv[2] if len(sys.argv) > 2 else '/tmp/ui_stands'
os.makedirs(out, exist_ok=True); log = []
with sync_playwright() as p:
    b = p.chromium.launch(args=["--use-gl=angle", "--use-angle=swiftshader", "--enable-unsafe-swiftshader"])
    pg = b.new_page(viewport={'width': 1500, 'height': 950})
    errs = []; pg.on('pageerror', lambda e: errs.append(str(e))); pg.on('console', lambda m: m.type == 'error' and errs.append(m.text))
    pg.goto(base + '/', wait_until='networkidle'); time.sleep(3)
    pg.evaluate("document.getElementById('stGrp').open = true")
    pg.evaluate("(() => { const e = document.getElementById('stOn'); e.checked = true; e.dispatchEvent(new Event('change')); })()")
    pg.wait_for_function("document.getElementById('stStatus').textContent.includes('stands')", timeout=120000); time.sleep(4)
    log.append('status: ' + pg.inner_text('#stStatus'))
    pg.screenshot(path=f'{out}/1_stands_status.png')
    pg.select_option('#stColor', 'BALIVE'); time.sleep(3); pg.screenshot(path=f'{out}/2_stands_ba.png')
    pg.select_option('#stColor', 'status'); time.sleep(1)
    # click the centre of a rendered forest stand: query a feature under the canvas centre region
    pt = pg.evaluate("""() => { const c = map.getCanvas(); for (let dy = 0; dy < 200; dy += 10) for (let dx = 0; dx < 200; dx += 10) {
        const x = c.width / 2 / devicePixelRatio + dx, y = c.height / 2 / devicePixelRatio + dy;
        const f = map.queryRenderedFeatures([x, y], { layers: ['stands-fill'] }); if (f.length && f[0].properties.status === 1) return [x, y, f[0].properties.stand_id]; } return null; }""")
    log.append(f'click target {pt}')
    box = pg.query_selector('canvas.maplibregl-canvas').bounding_box()
    pg.mouse.click(box['x'] + pt[0], box['y'] + pt[1]); time.sleep(2)
    pg.screenshot(path=f'{out}/3_popup.png')
    log.append('popup: ' + pg.inner_text('.stpop')[:300].replace('\n', ' | '))
    pg.click('#stSim'); pg.wait_for_function("!document.getElementById('status') || !document.getElementById('status').textContent.includes('resolving')", timeout=180000); time.sleep(6)
    log.append('app status: ' + pg.evaluate("(document.getElementById('status')||{}).textContent || ''")[:200])
    log.append('cards: ' + pg.evaluate("(document.getElementById('cards')||{}).innerText || ''")[:300].replace('\n', ' | '))
    pg.screenshot(path=f'{out}/4_stand_as_aoi.png')
    log.append('page errors: ' + (' || '.join(errs) if errs else 'none'))
    b.close()
print('\n'.join(log))
