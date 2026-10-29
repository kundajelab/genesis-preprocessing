# Workflow maps

One animated map per workflow in `docs/images/`: `genesis_metro_map.svg` (DAP-seq) and
`genesis_{atac,scatac,scatac_reads,rna,tss}_metro_map.svg`.

Routes are data dependencies, not run order; crossings without a station do not join channels.
Motion stops under reduced-motion settings and in print, and the colors adapt to dark mode.

To update a map, edit `render_metro.py` with its workflow change and run
`pixi run uv run --frozen --project genesis_tools python docs/diagrams/render_metro.py`.
`tests/verify_metro.py` checks that the maps are current and that every process has a station.
