"""AI provider abstraction: mock, OpenAI-compatible, Ollama (spec 6.5).

- No API key / no model service -> system still works via the deterministic
  mock provider (rule-based review, clearly labelled model_name="mock").
- Vision is attempted only when provider.supports_vision AND screenshots
  exist AND (local provider OR explicit allow_remote_screenshot_transfer).
  If images are not sent, visual review is never claimed.
"""
from __future__ import annotations

import base64
import json
import os
import re
from pathlib import Path
from typing import Any, Optional, Protocol

import httpx

from .prompts import build_review_prompt, PROMPT_VERSION


class AIProviderError(Exception):
    pass


class AIProvider(Protocol):
    name: str
    model: str
    supports_vision: bool
    is_local: bool

    def review_variant(self, task_ctx: dict, evidence: dict,
                       rule_hits: list[dict],
                       image_paths: list[str]) -> dict: ...


def _extract_json(text: str) -> dict:
    """Best-effort JSON extraction; raises AIProviderError on garbage."""
    try:
        return json.loads(text)
    except json.JSONDecodeError:
        pass
    m = re.search(r"\{.*\}", text, re.DOTALL)
    if m:
        try:
            return json.loads(m.group(0))
        except json.JSONDecodeError:
            pass
    raise AIProviderError(f"Model output was not valid JSON: {text[:300]!r}")


class MockAIProvider:
    """Deterministic rule/evidence-based reviewer (no external calls).

    Produces the same structured output as a real model so the entire
    downstream pipeline is testable offline. Results are ALWAYS labelled
    with model_name='mock-reviewer' - never presented as a real LLM verdict.
    """
    name = "mock"
    model = "mock-reviewer-v1"
    supports_vision = False
    is_local = True

    def review_variant(self, task_ctx: dict, evidence: dict,
                       rule_hits: list[dict], image_paths: list[str]) -> dict:
        ev = evidence
        thr = ev.get("_thresholds") or {}
        min_dp = thr.get("min_dp", 10)
        min_alt_ad = thr.get("min_alt_ad", 3)
        alt = ev.get("alt_count") or 0
        ref = ev.get("ref_count") or 0
        total_reads = ev.get("total_reads") or 0
        sri = ev.get("strand_ratio_imbalance")
        mq = ev.get("mean_mapq_alt")
        bq = ev.get("mean_baseq_alt")
        cluster = ev.get("alt_start_cluster_frac")
        ad_ok = ev.get("ad_consistency")
        vaf_dev = ev.get("vaf_deviation")

        supporting, conflicting, artifacts, reasons = [], [], [], []
        strong_strand_bias = False
        high_severity_signals = [h for h in rule_hits if h.get("severity") == "high"]
        if alt >= 5 and (sri is None or sri < 0.9) and (mq is None or mq >= 30):
            supporting.append(
                f"{alt} usable ALT reads with balanced strand (ratio {sri}) and mean MAPQ {mq}")
        if ev.get("homopolymer_len_spanning") and ev["homopolymer_len_spanning"] >= 6 and \
                task_ctx.get("variant_type") in ("INS", "DEL"):
            artifacts.append("indel in/near homopolymer - representation risk")
        if sri is not None and sri >= 0.9 and alt > 0:
            conflicting.append(f"ALT reads one-sided on strand (imbalance {sri:.2f})")
            reasons.append("STRAND_BIAS")
            if alt >= 5:
                # a strong, one-sided ALT pile with several independent reads
                # is the classic strand-bias artifact
                strong_strand_bias = True
        if cluster is not None and cluster > 0.5 and alt > 0:
            conflicting.append(f"{cluster:.0%} of ALT reads share the same alignment start")
            reasons.append("START_POS_CLUSTER")
        if ad_ok is False:
            conflicting.append("VCF AD disagrees with recomputed pileup counts")
            reasons.append("AD_PILEUP_MISMATCH")
        if vaf_dev is not None and vaf_dev > 0.35:
            conflicting.append(f"VAF deviates {vaf_dev:.2f} from genotype expectation")
            reasons.append("VAF_GT_MISMATCH")
        if mq is not None and mq < 20:
            conflicting.append(f"low mean ALT MAPQ {mq}")
            reasons.append("LOW_MQ")
        if bq is not None and bq < 13:
            conflicting.append(f"low mean ALT base quality {bq}")
            reasons.append("LOW_ALT_BQ")

        if total_reads == 0 or (alt + ref) == 0:
            recommendation = "INSUFFICIENT_EVIDENCE"
            summary = "No usable reads overlap the variant locus."
            reasons.append("NO_COVERAGE")
        elif total_reads < min_dp or alt < min_alt_ad:
            # below the configured depth/support floor: not enough data to
            # judge authenticity either way (NOT an artifact verdict)
            recommendation = "INSUFFICIENT_EVIDENCE"
            summary = (f"Coverage/support below the configured floor "
                       f"(DP {total_reads} < {min_dp}, ALT reads {alt} < {min_alt_ad}).")
            reasons.append("LOW_ALT_AD" if alt < min_alt_ad else "LOW_DP")
        elif len(conflicting) >= 2 and alt <= 6:
            recommendation = "SUSPECTED_FALSE_POSITIVE"
            summary = (f"Only {alt} ALT reads; multiple artifact signals: "
                      + "; ".join(conflicting))
        elif len(conflicting) >= 1 and alt <= 3:
            recommendation = "SUSPECTED_FALSE_POSITIVE"
            summary = f"Weak ALT support ({alt} reads) plus artifact signals: " + "; ".join(conflicting)
        elif strong_strand_bias:
            recommendation = "SUSPECTED_FALSE_POSITIVE"
            summary = (f"{alt} ALT reads but all on one strand (imbalance {sri:.2f}) - "
                       "classic strand-bias artifact.")
        elif len(conflicting) >= 2:
            recommendation = "INCONCLUSIVE"
            summary = "Evidence conflicts: " + "; ".join(conflicting)
        elif len(supporting) == 0 and alt < 3:
            recommendation = "INSUFFICIENT_EVIDENCE"
            summary = f"Only {alt} ALT reads available - below any defensible support level."
            reasons.append("LOW_ALT_AD")
        elif supporting and not artifacts:
            recommendation = "LIKELY_TRUE_VARIANT"
            summary = f"{alt} ALT vs {ref} REF reads; " + "; ".join(supporting)
        elif supporting and artifacts:
            # representation-risk artifacts (e.g. indel in a homopolymer)
            # prevent a confident "true variant" call
            recommendation = "INCONCLUSIVE"
            summary = (f"{alt} ALT vs {ref} REF reads, but representation-risk signals "
                       f"are present: " + "; ".join(artifacts))
        elif supporting and high_severity_signals:
            # the quantitative support is there, but a high-severity signal
            # (e.g. caller/pileup discrepancy) contradicts it -> not a
            # confident "true variant"
            recommendation = "INCONCLUSIVE"
            summary = (f"{alt} ALT vs {ref} REF reads, but high-severity signal(s) "
                       f"{', '.join(h.get('rule_id', '?') for h in high_severity_signals)} "
                       f"contradict the quantitative support: " + "; ".join(conflicting))
        else:
            recommendation = "INCONCLUSIVE"
            summary = "Evidence does not clearly resolve variant authenticity."

        hit_ids = [h.get("rule_id") for h in rule_hits if h.get("rule_id")]
        reasons = list(dict.fromkeys(reasons + [r for r in hit_ids if r not in reasons]))

        return {
            "ai_recommendation": recommendation,
            "evidence_summary": summary,
            "visual_findings": {},
            "artifact_signals": artifacts,
            "supporting_evidence": supporting,
            "conflicting_evidence": conflicting,
            "limitations": ["Mock deterministic reviewer - no vision model was used."],
            "reason_codes": reasons,
        }


class OpenAICompatibleProvider:
    name = "openai_compatible"
    is_local = False

    def __init__(self, cfg):
        from ..config import AIProviderConfig
        if not isinstance(cfg, AIProviderConfig):
            raise TypeError("expected AIProviderConfig")
        self.cfg = cfg
        self.model = cfg.model
        self.supports_vision = cfg.supports_vision
        if not cfg.base_url:
            raise AIProviderError("openai_compatible provider requires base_url")
        if cfg.api_key_env and not os.environ.get(cfg.api_key_env):
            raise AIProviderError(
                f"API key env var {cfg.api_key_env} is not set. Set it in the "
                "environment - keys are never written to config files or logs.")

    def _send(self, system: str, user: str, images_b64: list[tuple[str, str]]) -> str:
        content: list[dict[str, Any]] = [{"type": "text", "text": user}]
        for mime, b64 in images_b64:
            content.append({
                "type": "image_url",
                "image_url": {"url": f"data:{mime};base64,{b64}"},
            })
        payload = {
            "model": self.cfg.model,
            "temperature": self.cfg.temperature,
            "max_tokens": self.cfg.max_output_tokens,
            "messages": [
                {"role": "system", "content": system},
                {"role": "user", "content": content if images_b64 else user},
            ],
        }
        headers = {}
        if self.cfg.api_key_env and os.environ.get(self.cfg.api_key_env):
            headers["Authorization"] = f"Bearer {os.environ[self.cfg.api_key_env]}"
        url = self.cfg.base_url.rstrip("/") + "/chat/completions"
        last_err = None
        for attempt in range(self.cfg.retries + 1):
            try:
                resp = httpx.post(url, json=payload, headers=headers,
                                  timeout=self.cfg.timeout_s)
                resp.raise_for_status()
                data = resp.json()
                return data["choices"][0]["message"]["content"]
            except (httpx.HTTPError, KeyError, ValueError) as e:
                last_err = e
        raise AIProviderError(f"OpenAI-compatible API call failed after retries: {last_err}")

    def review_variant(self, task_ctx, evidence, rule_hits, image_paths):
        images_b64 = []
        if self.supports_vision and image_paths and self.cfg.allow_remote_screenshot_transfer:
            images_b64 = [encode_image(p) for p in image_paths]
        system, user = build_review_prompt(task_ctx, evidence, rule_hits,
                                           images_available=bool(images_b64))
        text = self._send(system, user, images_b64)
        out = _extract_json(text)
        if images_b64:
            out["_vision_used"] = True
        return out


class OllamaProvider:
    name = "ollama"
    is_local = True

    def __init__(self, cfg):
        from ..config import AIProviderConfig
        if not isinstance(cfg, AIProviderConfig):
            raise TypeError("expected AIProviderConfig")
        self.cfg = cfg
        self.model = cfg.model
        self.supports_vision = cfg.supports_vision
        self.base_url = (cfg.base_url or "http://127.0.0.1:11434").rstrip("/")

    def review_variant(self, task_ctx, evidence, rule_hits, image_paths):
        images_b64 = []
        if self.supports_vision and image_paths:
            # local provider: screenshots stay on the local machine
            images_b64 = [encode_image(p) for p in image_paths]
        system, user = build_review_prompt(task_ctx, evidence, rule_hits,
                                           images_available=bool(images_b64))
        payload = {
            "model": self.model,
            "format": "json",
            "stream": False,
            "options": {
                "temperature": self.cfg.temperature,
                "num_predict": self.cfg.max_output_tokens,
            },
            "messages": [
                {"role": "system", "content": system},
                {"role": "user", "content": user,
                 **({"images": [b for _, b in images_b64]} if images_b64 else {})},
            ],
        }
        last_err = None
        for _ in range(self.cfg.retries + 1):
            try:
                resp = httpx.post(self.base_url + "/api/chat", json=payload,
                                  timeout=self.cfg.timeout_s)
                resp.raise_for_status()
                out = _extract_json(resp.json()["message"]["content"])
                if images_b64:
                    out["_vision_used"] = True
                return out
            except (httpx.HTTPError, KeyError, ValueError) as e:
                last_err = e
        raise AIProviderError(f"Ollama API call failed after retries: {last_err}")


def encode_image(path: str | Path, max_bytes: int = 8 * 1024 * 1024) -> tuple[str, str]:
    p = Path(path)
    if not p.exists():
        raise AIProviderError(f"Image not found: {p}")
    if p.stat().st_size < 100:
        raise AIProviderError(f"Image suspiciously small, refusing to send: {p}")
    if p.stat().st_size > max_bytes:
        raise AIProviderError(f"Image exceeds {max_bytes} bytes: {p}")
    mime = "image/png" if p.suffix == ".png" else "image/jpeg"
    b64 = base64.b64encode(p.read_bytes()).decode("ascii")
    return mime, b64


def make_provider(cfg) -> AIProvider:
    from ..config import AIProviderConfig
    if not isinstance(cfg, AIProviderConfig):
        raise TypeError("expected AIProviderConfig")
    if cfg.provider == "mock":
        return MockAIProvider()
    if cfg.provider == "openai_compatible":
        return OpenAICompatibleProvider(cfg)
    if cfg.provider == "ollama":
        return OllamaProvider(cfg)
    raise AIProviderError(f"Unknown provider: {cfg.provider}")
