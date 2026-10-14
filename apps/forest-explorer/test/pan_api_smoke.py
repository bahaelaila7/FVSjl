#!/usr/bin/env python3
"""End-to-end API smoke test of the Pan-enabled Forest Explorer (stdlib only).
   python3 test/pan_api_smoke.py [base_url=http://127.0.0.1:8088]
Uses the shipped samples (samples/index.json): Oconee AOI, the three GA+FL parameter sets, the hot & dry climate."""
import json, sys, time, urllib.request, os

BASE = sys.argv[1] if len(sys.argv) > 1 else "http://127.0.0.1:8088"
HERE = os.path.dirname(os.path.abspath(__file__)); SAMPLES = os.path.join(HERE, "..", "samples")
ok = True

def req(path, body=None, ctype="application/json", timeout=3600):
    data = None if body is None else (body if isinstance(body, bytes) else json.dumps(body).encode())
    r = urllib.request.Request(BASE + path, data=data, headers={"Content-Type": ctype} if data else {})
    try:
        with urllib.request.urlopen(r, timeout=timeout) as f:
            return f.status, json.loads(f.read() or b"null")
    except urllib.error.HTTPError as e:
        return e.code, json.loads(e.read() or b"null")

def check(name, cond, detail=""):
    global ok
    print(("PASS " if cond else "FAIL ") + name + (f" — {detail}" if detail else ""))
    ok &= bool(cond)

st, meta = req("/api/meta")
check("meta has pan block", st == 200 and "pan" in meta, json.dumps(meta.get("pan", {}))[:160])
for _ in range(120):                       # wait for the worker warm-up
    if req("/api/meta")[1]["pan"]["state"] in ("ready", "failed", "stopped"): break
    time.sleep(10)
st, idx = req("/api/pan/samples"); check("samples index", st == 200 and len(idx["params"]) >= 3)
ids = {}
for p in idx["params"]:
    y = open(os.path.join(SAMPLES, p["file"]), "rb").read()
    t = time.time(); st, r = req("/api/pan/params", y, "text/yaml")
    check(f"upload {p['file']}", st == 200 and r["plugin"] == p["plugin"], f"{r.get('name', r)} plugin={r.get('plugin')} warnings={r.get('warnings')} ({time.time()-t:.0f}s)")
    if st == 200: ids.setdefault(p["plugin"], r["params_id"])
st, r = req("/api/pan/params", b"pan_params: {plugin: nope}\n", "text/yaml"); check("bad YAML rejected", st == 400, r.get("error", "")[:100])

geom = json.load(open(os.path.join(SAMPLES, "aois", "oconee_ga.geojson")))["features"][0]["geometry"]
st, aoi = req("/api/aoi", {"geometry": geom}); check("AOI resolves", st == 200 and aoi["nplots"] > 0, f"{aoi.get('nplots')} plots, {aoi.get('acres')} ac")

def sim(engine, plan, climate=None, pid=None):
    body = {"geometry": geom, "plan": plan, "engine": engine}
    if pid: body["params_id"] = pid
    if climate is not None: body["climate"] = climate
    t = time.time(); st, r = req("/api/simulate", body); return st, r, time.time() - t

harvest = {"kind": "harvest", "cycle": 1, "selection": {"species": "all", "age_lo": 0, "age_hi": 999, "dbh_lo": 0, "dbh_hi": 999},
           "removal": {"type": "residual_ba", "value": 60, "direction": "below"}, "stands": {"rank": "all", "area_pct": 100}, "plant_after": None}
thin_plant = [{"kind": "thin", "cycle": 1, "metric": "BA", "target": 70, "direction": "below", "dbh_lo": 0, "dbh_hi": 999, "species": "all", "mode": "rule"},
              {"kind": "harvest", "cycle": 2, "selection": {"species": "all"}, "removal": {"type": "clearcut", "value": 0},
               "stands": {"rank": "oldest", "area_pct": 30}, "plant_after": {"species": "PITA", "tpa": 600, "survival": 80}}]
hot = json.load(open(os.path.join(SAMPLES, "climates", "hot_dry_ssp585.json")))
runs = [("pan-hbad", {"name": "hbad harvest", "cycles": 4, "period": 10, "actions": [harvest]}, None, ids.get("hbad")),
        ("pan-biomass", {"name": "bs thin+plant", "cycles": 4, "period": 10, "actions": thin_plant}, None, ids.get("biomass")),
        ("pan-pnet", {"name": "pnet baseline", "cycles": 4, "period": 10, "actions": []}, None, ids.get("pnet")),
        ("pan-pnet", {"name": "pnet hot dry", "cycles": 4, "period": 10, "actions": []}, hot, ids.get("pnet")),
        ("fvs", {"name": "fvs grow", "cycles": 3, "period": 10, "actions": []}, None, None)]
for eng, plan, cl, pid in runs:
    st, r, dt = sim(eng, plan, cl, pid)
    good = st == 200 and r.get("ncycles", 0) > 0
    d = (f"{r.get('nplots_sim')} plots, {dt:.0f}s; carbon {[c['carbon_live_ag'] for c in r.get('cycles', [])]}; removed vol {[c['removed_volume'] for c in r.get('cycles', [])]}; "
         f"warnings {[w[:80] for w in r.get('warnings') or []][:2]}") if good else str(r)[:300]
    check(f"simulate {eng} '{plan['name']}'", good, d)
    if good: check(f"  {eng} map image + species×DBH", r.get("image") is not None and len(r["dist"][0]) > 0, f"dist rows {len(r['dist'][0])}")
st, r = req("/api/simulate", {"geometry": geom, "plan": {"name": "x", "cycles": 2, "period": 10, "actions": []}, "engine": "pan-hbad", "params_id": ids.get("biomass")})
check("params/engine mismatch rejected", st == 400, str(r)[:160])
st, r = req("/api/pan/climate_preview", {"geometry": geom, "climate": hot, "start_year": 2025, "years": 80})
good = st == 200 and len(r.get("vars", [])) >= 4
check("climate preview", good, (f"L4s {r['l4s']}; " + "; ".join(f"{v['name']} {v['baseline_annual'][0]}→{v['modified_annual'][-1]} {v['annual_units']}" for v in r["vars"])) if good else str(r)[:300])
# Region 4 coverage: default parameter sets (Region 4) on an AOI outside GA/FL
g2 = json.load(open(os.path.join(SAMPLES, "aois", "uwharrie_nc.geojson")))["features"][0]["geometry"]
st, r = req("/api/simulate", {"geometry": g2, "plan": {"name": "nc", "cycles": 3, "period": 10, "actions": []}, "engine": "pan-hbad"})
good = st == 200 and r.get("nplots_sim", 0) > 0.7 * r.get("nplots_aoi", 1)
check("Region 4 default set covers an NC AOI", good, f"{r.get('nplots_sim')} of {r.get('nplots_aoi')} plots; params {r.get('engine_info', {}).get('params_name') if isinstance(r, dict) else r}; {(r.get('warnings') or [''])[0][:120] if isinstance(r, dict) else ''}")
print("ALL PASS" if ok else "SOME CHECKS FAILED")
sys.exit(0 if ok else 1)
