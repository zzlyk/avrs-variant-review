# SECURITY

Genetic data is sensitive health information. AVRS is **local-first** by
design; this document describes what the software enforces, what it does
*not* enforce, and what the deploying organisation must still do.

> AVRS does **not** claim HIPAA, CLIA, CAP, ISO 15189 or GDPR compliance.
> Compliance is a property of the deployment and the operating procedures,
> not of a Python package.

---

## 1. Data-handling defaults

| Control | Implementation | Where |
|---|---|---|
| Local-first processing | BAM/VCF/FASTQ never leave the machine unless an AI provider is explicitly configured with remote transfer | `security.local_first`, `ai.allow_remote_screenshot_transfer` |
| No screenshots to third parties by default | Remote transfer requires an explicit `true`; local providers (Ollama) keep images on the host | `avrs/ai/providers.py` |
| No API keys in config or logs | Keys are read from an environment variable *at call time* and never logged or persisted | `ai.api_key_env`, `avrs/ai/providers.py` |
| Pseudonymous sample IDs | Sample identifiers are operator-controlled; the demo uses `DEMO_SAMPLE_001`; `security.redact_names_in_requests` is available for request shaping | `config/security` |
| Raw data is immutable | The original BAM/VCF/FASTQ are never rewritten; AVRS only reads them and records their checksums | `avrs/variants`, `avrs/review` |
| No arbitrary command execution from data | VCF/annotation text is treated as untrusted data; it is never passed to a shell | `avrs/ai/prompts.py` |
| Bounded model context | Only one variant plus filtered evidence is sent per request; the full VCF is never uploaded | `security.max_llm_context_variants` |
| Audit trail | Every state-changing action is written to `audit_logs` (actor, action, entity, timestamp, details) | `avrs/db`, `avrs/review` |

## 2. Prompt-injection posture

VCF INFO fields, annotations and free-text sample metadata are **data**, not
instructions. The review prompt states this explicitly, and the provider
layer sends only a structured JSON evidence block. The model never receives:

* shell access,
* arbitrary filesystem read/write,
* credentials,
* the raw FASTQ/BAM.

Recommended extension: keep it that way when adding new providers or
annotation sources.

## 3. Threat model (summary)

| Threat | Mitigation | Residual risk |
|---|---|---|
| Data exfiltration via the model | remote transfer opt-in, local providers preferred, explicit `visual_review_performed` flag | A misconfigured `allow_remote_screenshot_transfer: true` will send images |
| API key leakage | env-var only, never written to config/DB/logs | Keys in shell history or process lists |
| Untrusted annotation → prompt injection | annotations framed as untrusted data; schema validation of the response | A determined injection could still bias the text output; the deterministic evidence layer is unaffected |
| Cross-sample contamination of screenshots | per-sample/per-variant snapshot directories + metadata binding (sample, variant, build, checksums) | Human mis-interpretation of a manually copied file |
| Concurrent IGV control | `IGVInstanceLock` (file lock) serialises access to the single IGV instance | Lock file left behind after a crash → stale-lock reclaim after `max_age_s` |
| Tampering with results | AI and human conclusions stored separately, human decisions are append-only with a supersede chain | An operator with DB write access can still alter rows — protect the DB with OS permissions |
| Sensitive data in filenames/logs | operator-controlled sample IDs, log level configurable, no patient identifiers generated | Operator may still choose identifiable names |

## 4. Deployment responsibilities (not implemented in the MVP)

The following **must** be provided by the deploying organisation:

1. **Authentication & authorisation** — the workbench currently has no login;
   bind it to localhost or put it behind an authenticating reverse proxy.
2. **Transport security** — TLS for anything crossing a network boundary.
3. **Encryption at rest** — filesystem/DB encryption for the output directory.
4. **Key management** — a secrets manager for API keys, not shell exports.
5. **Access logging & retention policy** — how long reviews and screenshots
   are kept, and who may read them.
6. **Backup & recovery** — the SQLite DB plus the snapshot directory are the
   audit record; back them up together.
7. **Network segmentation** — the IGV port listener (`60151`) should not be
   exposed outside the analysis workstation.
8. **Incident response** — a documented procedure for data exposure.

## 5. IGV-specific notes

* IGV Desktop's command listener has **no authentication**. Keep it bound to
  `127.0.0.1` and reachable only from the analysis host.
* IGV caches loaded data in its own session; use `avrs igv check` to confirm
  the instance is the one you expect, and `new` (reset) between unrelated
  samples.
* Snapshots contain visual information about the sample's reads — treat them
  with the same protection as the BAM itself.

## 6. Model-specific notes

* **Local models (Ollama)** keep screenshots on the host; verify the Ollama
  server is not exposed to a wider network.
* **Remote APIs** are a data transfer to a third party. Record the vendor,
  the data contract and the retention policy in your SOP before enabling
  `allow_remote_screenshot_transfer`.
* Record `model_name`, `model_version` and `prompt_version` for every review
  (AVRS does this automatically) so a result can be reproduced or challenged.
* Prompts and model outputs can be wrong. AVRS stores the raw JSON result
  next to the recommendation for exactly this reason.

## 7. Reporting a security issue

Please report suspected vulnerabilities privately to the repository owner
rather than in a public issue. Include a description, reproduction steps and
the affected version (`avrs version`).
