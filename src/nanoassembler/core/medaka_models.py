"""Translate Dorado/Guppy basecalling model ids into Medaka model names.

Dorado stamps the basecalling model into every FASTQ header, e.g.
    dna_r10.4.1_e8.2_400bps_sup@v5.0.0
Medaka names the matching consensus model
    r1041_e82_400bps_sup_v5.0.0
so most of the work is a mechanical rewrite; the rest is picking the closest
available model when the exact one is not shipped with the installed Medaka.
"""
from __future__ import annotations

import re
import subprocess
from functools import lru_cache
from typing import Optional

_SPEEDS = ("fast", "hac", "sup")

# R9.4.1 predates the systematic naming, so the mechanical rewrite misses.
# These map to the MinION/GridION models; PromethION runs should be switched to
# the r941_prom_* equivalent with the per-sample override in the UI.
_ALIASES = {
    "r941_450bps_fast": "r941_min_fast_g507",
    "r941_450bps_hac": "r941_min_hac_g507",
    "r941_450bps_sup": "r941_min_sup_g507",
    "r103_450bps_fast": "r103_fast_g507",
    "r103_450bps_hac": "r103_hac_g507",
    "r103_450bps_sup": "r103_sup_g507",
}


def _apply_alias(name: str) -> tuple[str, str]:
    """Return (name, note) after applying the legacy-chemistry alias table."""
    base = re.sub(r"_(v[\d.]+|g\d+)$", "", name)
    target = _ALIASES.get(base)
    if not target:
        return name, ""
    return target, (
        f"{name} has no direct Medaka equivalent; using {target} "
        f"(MinION/GridION - override per sample for PromethION)"
    )


def dorado_to_medaka(model_id: str) -> Optional[str]:
    """Rewrite a Dorado model id into the Medaka naming convention."""
    if not model_id:
        return None
    name = model_id.strip()
    name = re.sub(r"^dna_", "", name)
    # r10.4.1 -> r1041, e8.2 -> e82 (digits keep their order, dots go away).
    # "_" is a word character, so \b cannot be used to close these patterns.
    edge = r"(?<![A-Za-z0-9])"
    name = re.sub(edge + r"r(\d+)\.(\d+)\.(\d+)(?![\d.])", r"r\1\2\3", name)
    name = re.sub(edge + r"r(\d+)\.(\d+)(?![\d.])", r"r\1\2", name)
    name = re.sub(edge + r"e(\d+)\.(\d+)(?![\d.])", r"e\1\2", name)
    name = name.replace("@v", "_v")
    name = re.sub(r"_+", "_", name).strip("_")
    if not re.search(r"(fast|hac|sup)", name):
        return None
    return name


def parse_components(medaka_name: str) -> dict[str, str]:
    """Split a Medaka-style name into chemistry / speed / version parts."""
    out: dict[str, str] = {}
    m = re.match(
        r"^(?P<chem>.+?)_(?P<speed>fast|hac|sup)(?:_variant)?(?:_(?P<rest>.*))?$",
        medaka_name,
    )
    if not m:
        return out
    out["chem"] = m.group("chem")
    out["speed"] = m.group("speed")
    rest = m.group("rest") or ""
    v = re.search(r"v[\d.]+", rest)
    out["version"] = v.group(0) if v else ""
    out["rest"] = rest
    return out


# Models that need information a FASTQ does not carry (dwell times / move
# tables) or that are trained for a specific protocol. They are never picked as
# an automatic fallback, only if the user names one explicitly.
_SPECIAL = ("_rl_lstm", "dwells", "bacterial_methylation", "joint_")


def is_general_purpose(name: str) -> bool:
    return not any(tok in name for tok in _SPECIAL)


def to_variant(name: str) -> str:
    """Consensus model name -> the matching variant-calling model name."""
    if "_variant" in name:
        return name
    return re.sub(r"_(v[\d.]+|g\d+)$", r"_variant_\1", name)


@lru_cache(maxsize=4)
def available_models(medaka_bin: str) -> tuple[str, ...]:
    """Ask the installed Medaka which models it can actually run."""
    try:
        proc = subprocess.run(
            [medaka_bin, "tools", "list_models"],
            capture_output=True,
            text=True,
            timeout=120,
        )
    except (OSError, subprocess.SubprocessError):
        return ()
    text = proc.stdout + proc.stderr
    models: list[str] = []
    for line in text.splitlines():
        if ":" in line and "models" in line.split(":")[0].lower():
            line = line.split(":", 1)[1]
        for token in re.split(r"[,\s]+", line):
            token = token.strip().strip("'\"[]")
            if re.match(r"^r(?:\d+|na\d+)[\w.]*_(?:fast|hac|sup)[\w.]*$", token):
                models.append(token)
    return tuple(dict.fromkeys(models))


def resolve(
    model_id: str, medaka_bin: Optional[str] = None, kind: str = "consensus"
) -> tuple[Optional[str], str, str]:
    """Map a basecalling model to an installed Medaka model.

    ``kind`` is "consensus" or "variant"; Medaka ships a separate, differently
    named network for variant calling.
    Returns (model_name, source, note) where source is one of
    'auto' (exact match), 'fallback' (closest available) or 'none'.
    """
    candidate = dorado_to_medaka(model_id)
    alias_note = ""
    if candidate:
        candidate, alias_note = _apply_alias(candidate)
        if kind == "variant":
            candidate = to_variant(candidate)
    if not candidate:
        return None, "none", f"could not interpret basecall model {model_id!r}"
    if not medaka_bin:
        # No Medaka to interrogate; hand back the mechanical translation so the
        # UI can still show what would be used.
        return candidate, "auto", "medaka not installed; model not verified"

    installed = available_models(medaka_bin)
    if not installed:
        return candidate, "auto", "could not list medaka models; using derived name"
    if candidate in installed:
        return candidate, "auto", alias_note

    want = parse_components(candidate)
    if not want:
        return None, "none", f"no medaka model matches {candidate!r}"

    def version_key(name: str) -> tuple:
        v = parse_components(name).get("version", "")
        parts = re.findall(r"\d+", v)
        return tuple(int(p) for p in parts) if parts else (0,)

    # Same chemistry and speed, highest version we have.
    want_variant = "_variant" in candidate
    pool = [
        m
        for m in installed
        if is_general_purpose(m) and ("_variant" in m) == want_variant
    ]
    same = [
        m
        for m in pool
        if parse_components(m).get("chem") == want.get("chem")
        and parse_components(m).get("speed") == want.get("speed")
    ]
    if same:
        best = max(same, key=version_key)
        return best, "fallback", f"exact model {candidate} unavailable; using {best}"

    # Same chemistry, step down the accuracy ladder.
    order = {s: i for i, s in enumerate(_SPEEDS)}
    same_chem = [m for m in pool if parse_components(m).get("chem") == want.get("chem")]
    if same_chem:
        best = max(
            same_chem,
            key=lambda m: (order.get(parse_components(m).get("speed", ""), -1), version_key(m)),
        )
        return best, "fallback", f"no {want.get('speed')} model for this chemistry; using {best}"

    return None, "none", f"no medaka model matches {candidate!r}"
