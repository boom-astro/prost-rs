# BOOM Integration

prost can be used with the BOOM alert broker to associate transients
with host galaxies from Legacy Survey cross-matches.

## Querying BOOM for Galaxy Candidates

```python
import requests

BOOM_URL = "https://api.kaboom.caltech.edu"

# Authenticate
session = requests.Session()
resp = session.post(f"{BOOM_URL}/auth", data={
    "username": "YOUR_USERNAME",
    "password": "YOUR_PASSWORD",
})
token = resp.json()["access_token"]
session.headers["Authorization"] = f"Bearer {token}"

# Get transient coordinates
pipeline = [
    {"$match": {"objectId": "ZTF20abjbgjj"}},
    {"$group": {
        "_id": "$objectId",
        "ra": {"$first": "$candidate.ra"},
        "dec": {"$first": "$candidate.dec"},
    }},
]
data = session.post(f"{BOOM_URL}/queries/pipeline", json={
    "catalog_name": "ZTF_alerts",
    "pipeline": pipeline,
}).json()["data"]

ra, dec = data[0]["ra"], data[0]["dec"]
```

## Cross-matching with Legacy Survey

```python
# Cone search LS_DR10 for galaxies within 60"
matches = session.post(f"{BOOM_URL}/queries/cone_search", json={
    "catalog_name": "LS_DR10",
    "radius": 60.0,
    "unit": "Arcseconds",
    "object_coordinates": {"ZTF20abjbgjj": [ra, dec]},
    "projection": {
        "ra": 1, "dec": 1, "type": 1,
        "shape_r": 1, "shape_e1": 1, "shape_e2": 1, "z": 1,
    },
}).json()["data"]["ZTF20abjbgjj"]
```

## Running Association

```python
from prost import Transient, GalaxyCandidate, associate_host

transient = Transient(ra=ra, dec=dec, name="ZTF20abjbgjj")

galaxies = []
for doc in matches:
    if doc.get("type") == "PSF":
        continue  # Skip stars
    try:
        g = GalaxyCandidate.from_tractor(
            ra=doc["ra"], dec=doc["dec"],
            shape_r=doc["shape_r"],
            shape_e1=doc.get("shape_e1", 0.0),
            shape_e2=doc.get("shape_e2", 0.0),
            redshift=doc.get("z"),
            objtype=doc.get("type"),
            objname=str(doc.get("_id", "")),
        )
        galaxies.append(g)
    except ValueError:
        continue  # Skip invalid shapes

result = associate_host(transient, galaxies)

if result.best_host:
    h = result.best_host
    print(f"Host: {h.galaxy.objname}")
    print(f"Separation: {h.separation_arcsec:.1f}\"")
    print(f"FO: {h.fractional_offset:.2f}")
    print(f"P(host): {h.posterior:.3f}")
else:
    print(f"No host found (P(none) = {result.p_none:.3f})")
```

!!! note "BOOM already runs DLR association"
    BOOM's enrichment pipeline includes a built-in DLR + Bayesian host
    association (in `src/utils/host.rs`). The `host_galaxy` field on the
    aux document contains the result. Use prost for offline analysis,
    batch reprocessing, or when you need redshift-aware scoring.
