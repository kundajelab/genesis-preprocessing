"""Check reproducible map generation, accessibility and workflow coverage."""

from __future__ import annotations

import hashlib
import re
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
NS = {"s": "http://www.w3.org/2000/svg"}


def main() -> None:
    with tempfile.TemporaryDirectory(prefix="genesis-maps-") as folder:
        out = Path(folder)
        subprocess.run(
            [sys.executable, str(ROOT / "docs/diagrams/render_metro.py"), "--output-dir", folder],
            check=True,
        )
        assert len(list(out.glob("*.svg"))) == 7
        for generated in out.glob("*.svg"):
            source = ROOT / "docs/images" / generated.name
            assert generated.read_bytes() == source.read_bytes(), f"Stale map: {source}"
            svg = ET.parse(source).getroot()
            ids = [e.attrib["id"] for e in svg.iter() if "id" in e.attrib]
            assert len(ids) == len(set(ids))
            assert svg.attrib["aria-labelledby"] == "title description"
            assert svg.find("s:title", NS) is not None
            assert svg.find("s:desc", NS) is not None
            assert not svg.findall(".//s:script", NS) and not svg.findall(".//s:image", NS)
            text = source.read_text()
            assert "prefers-reduced-motion:reduce" in text
            assert "prefers-color-scheme:dark" in text
        for stem, workflow in [
            ("genesis", "workflows/dap-seq.nf"),
            ("genesis_atac", "workflows/atac.nf"),
        ]:
            static = ET.parse(out / f"{stem}_metro_map.svg").getroot()
            animated = ET.parse(out / f"{stem}_metro_map_animated.svg").getroot()

            def paths(svg: ET.Element) -> list[str]:
                return [p.attrib["d"] for p in svg.findall("s:path", NS) if "id" in p.attrib]

            assert paths(static) == paths(animated)
            code = (ROOT / workflow).read_text()
            expected = set(re.findall(r"include\s*\{\s*(\w+)\s*\}\s*from '\.\./modules/", code))
            stations = {g.attrib["id"] for g in static.findall("s:g", NS)}
            assert expected <= stations, expected - stations
            assert (
                hashlib.sha256(code.encode()).hexdigest()
                in (out / f"{stem}_metro_map.svg").read_text()
            )
    print("PASS: seven deterministic accessible maps; actual DAP/ATAC process coverage")


if __name__ == "__main__":
    main()
