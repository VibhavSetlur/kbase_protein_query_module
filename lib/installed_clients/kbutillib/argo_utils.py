"""KBase SDK utilities for working with KBase SDK environments and services."""

import logging
import os
import random
import re
import time
from typing import Any, Optional

import httpx

from .shared_env_utils import SharedEnvUtils

SYSTEM_MSG = (
    "You are an expert curator of gene functional annotations (SwissProt or RAST style). "
    "When you compare two annotations, classify their relationship according to "
    "ontology practice (EC, UniProt, KEGG, phage manuals). "
    "Respond in the format: Label — short justification (evidence or rule number)."
)

_ALIAS = {
    "Identical": "Exact",
    "Same": "Exact",
    "Equivalent": "Synonym",
    "Similar": "Synonym",
    "Partial": "Related",
}
_DASH = re.compile(r"\s*[—–-]\s*")

# Valid labels for annotation comparison
LABEL_SET = {"Exact", "Synonym", "Related", "Unknown"}

# Create logger
logger = logging.getLogger(__name__)

__all__ = ["ArgoUtils", "llm_label"]

# ---------------------------------------------------------------------------
# Constants & helpers
# ---------------------------------------------------------------------------

# Models that exist on both prod & dev – we prefer prod but can fall back
DUAL_ENV_MODELS = {"gpt4o", "gpt4olatest"}

# o-series models (reasoning family) need more time to respond
O_SERIES_TIMEOUT = 120.0  # sec – max per request
DEFAULT_TIMEOUT = 30.0

# HTTP codes indicating the request is still being processed
_PROCESSING = {102, 202}

# Polling settings (used only when we get 102/202)
POLL_EVERY = 3.0  # sec


class ArgoUtils(SharedEnvUtils):
    """Utilities for interfacing with Argo to run LLM queries."""

    def __init__(
        self,
        model: str = "gpto3mini",
        env: Optional[str] = None,
        stream: Optional[bool] = None,
        user: Optional[str] = None,
        api_key: Optional[str] = None,
        proxy_port: Optional[int] = 1080,
        timeout: Optional[float] = None,
        retries: int = 5,
        **kwargs: Any,
    ) -> None:
        """Initialize Argo utilities."""
        super().__init__(**kwargs)
        # ------------------------------------------------------------------
        # 1. Decide environment (prod vs dev)
        #    • caller can override via *env*
        #    • o-series default to dev
        #    • dual-env models: try prod first, fall back to dev on first 5xx
        # ------------------------------------------------------------------
        self._extra = {}
        self.model = model
        self.env = env or ("dev" if model.startswith("gpto") else "prod")

        # build helper to re-construct URLs when env flips later
        def _base_url(e: str) -> str:
            return (
                f"https://apps{'-' + e if e != 'prod' else ''}"
                ".inside.anl.gov/argoapi/api/v1/resource/"
            )

        self._base_url_fn = _base_url  # keep for later reuse

        base = _base_url(self.env)

        # Decide streaming endpoint
        stream_capable = {"gpto3mini"}  # expand when other o-series gain SSE
        self._stream = stream if stream is not None else (model in stream_capable)
        self.url = base + ("streamchat/" if self._stream else "chat/")

        # ------------------------------------------------------------------
        # 2. Auth & identity
        # ------------------------------------------------------------------
        self.user = user or os.getenv("ARGO_USER") or os.getlogin()
        self.retries = retries

        # optional kwargs (e.g. temperature for vote mode)
        self.temperature = kwargs.get("temperature")

        # ------------------------------------------------------------------
        # 3. Compute HTTP timeout (must be BEFORE first use of self.timeout)
        # ------------------------------------------------------------------
        # Compute timeout: caller override wins; else pick model-based default
        self.timeout = (
            timeout
            if (timeout is not None and timeout > 0)
            else (O_SERIES_TIMEOUT if model.startswith("gpto") else DEFAULT_TIMEOUT)
        )

        # enable verbose logging when ARGO_DEBUG=1 or debug kwarg
        if os.getenv("ARGO_DEBUG") or kwargs.get("debug"):
            logger.setLevel(logging.DEBUG)
            self.log_debug("Debug mode active for ArgoGatewayClient")

        self.log_info(
            f"ArgoGatewayClient initialised | model={self.model} env={self.env} "
            f"timeout={self.timeout:.1f}s url={self.url}"
        )

        # Instantiate HTTP client using socks5 proxy if needed
        proxies = None
        if proxy_port:
            proxies = {
                "http://": f"socks5://127.0.0.1:{proxy_port}",
                "https://": f"socks5://127.0.0.1:{proxy_port}",
            }
        self.cli = httpx.Client(
            proxies=proxies, timeout=self.timeout, follow_redirects=True
        )

        # headers (api key optional)
        self.headers = {"Content-Type": "application/json"}
        if api_key or os.getenv("ARGO_API_KEY"):
            self.headers["x-api-key"] = api_key or os.getenv("ARGO_API_KEY")

        # ------------------------------------------------------------------
        # 4. Dual-env prod→dev quick check (single ping)
        # ------------------------------------------------------------------
        if env is None and model in DUAL_ENV_MODELS and self.env == "prod":
            try:
                if not self.ping():
                    # switch to dev & rebuild URL
                    self.env = "dev"
                    base = _base_url("dev")
                    self.url = base + ("streamchat/" if self._stream else "chat/")
            except Exception:
                # network error → assume prod dead, flip to dev
                self.env = "dev"
                base = _base_url("dev")
                self.url = base + ("streamchat/" if self._stream else "chat/")

    @staticmethod
    def _extract_job_url(r: httpx.Response, base: str):  # helper
        """Return absolute URL that can be polled for job status/result."""
        if "Location" in r.headers:  # standard HTTP header
            loc = r.headers["Location"]
            return loc if loc.startswith("http") else base + loc.lstrip("/")
        try:
            j = r.json()
        except Exception:
            return None
        if "job_url" in j:
            return j["job_url"]
        if job_id := j.get("job_id"):
            return base + f"status/{job_id}"
        return None

    # ------------------------------------------------------------------
    def _payload(self, prompt: str, system: str) -> dict:
        if self.model.startswith("gpto") or self.model.startswith(
            "o"
        ):  # o-series models
            # Allow caller override; else default to a small number (32) to avoid
            # massive completions that sometimes trigger gateway bugs.
            payload = {
                "user": self.user,
                "model": self.model,
                "prompt": [prompt],
            }
            # NEW: only send if explicitly configured
            if "max_completion_tokens" in self._extra:
                payload["max_completion_tokens"] = int(
                    self._extra["max_completion_tokens"]
                )

            # o3 supports system; o1/o1-mini ignore it – harmless to include.
            if system:
                payload["system"] = system
            return payload
        return {  # GPT-style
            "user": self.user,
            "model": self.model,
            "messages": [
                {"role": "system", "content": system},
                {"role": "user", "content": prompt},
            ],
            "temperature": 0.0,
        }

    # ------------------------------------------------------------------
    def chat(self, prompt: str, *, system: str = "") -> str:
        """Send a chat request to the Argo LLM service.

        Args:
            prompt: The user prompt/question to send
            system: Optional system message for context

        Returns:
            The LLM response text
        """
        payload = self._payload(prompt, system)

        # Allow one automatic flip between /chat/ and /streamchat/ on blank reply
        endpoint_switched = False
        sentinel_injected = False

        for att in range(self.retries + 1):
            self.log_debug(f"Attempt {att + 1} POST → {self.url}")
            try:
                r = self.cli.post(self.url, json=payload, headers=self.headers)
            except httpx.TimeoutException:
                reason = "timeout"
                self.log_warning(f"Timeout on attempt {att + 1}")
                r = None
            else:
                reason = ""

            if r is not None:
                self.log_debug(f"Status {r.status_code} | body preview: {r.text[:120]}")
                # handling for server responses
                if r.status_code in _PROCESSING:
                    self.log_debug(
                        f"Processing accepted (HTTP {r.status_code}); enter poll loop"
                    )
                    txt = self._poll_for_result(r)
                    if txt:
                        self.log_debug(
                            f"Poll loop succeeded; returning text ({len(txt)} chars)"
                        )
                        return txt
                    reason = "processing timeout"
                elif r.status_code >= 500:
                    self.log_warning(f"5xx ({r.status_code}) on attempt {att + 1}")
                    # on prod failure & dual-env model, auto-retry once on dev
                    if (
                        self.model in DUAL_ENV_MODELS
                        and self.env == "prod"
                        and att == 0
                    ):
                        self.env = "dev"
                        self.url = self._base_url_fn("dev") + (
                            "streamchat/" if self._stream else "chat/"
                        )
                        continue  # retry immediately on dev
                    reason = f"{r.status_code}"
                else:
                    r.raise_for_status()
                    txt = self._extract_txt(r)
                    if txt:
                        self.log_debug(f"Received final text ({len(txt)} chars)")
                        return txt
                    reason = "blank reply"
                    self.log_warning(f"Blank reply body: {r.text[:200]}")

                    # (1) Flip endpoint once
                    if not endpoint_switched:
                        self._stream = not self._stream
                        base = self._base_url_fn(self.env)
                        self.url = base + ("streamchat/" if self._stream else "chat/")
                        payload = self._payload(prompt, system)
                        endpoint_switched = True
                        self.log_info(
                            f"Endpoint switched due to blank reply → {self.url}"
                        )
                        continue  # retry immediately without back-off

                    # (2) Inject sentinel once, force /chat/
                    if not sentinel_injected:
                        prompt = f"Label: {prompt}"
                        self._stream = False
                        base = self._base_url_fn(self.env)
                        self.url = base + "chat/"
                        payload = self._payload(prompt, system)
                        sentinel_injected = True
                        self.log_info(
                            'Sentinel injected → retry with /chat/ and "Label:" prefix'
                        )
                        continue  # retry again

            if att < self.retries:
                delay = 1.5 * 2**att + random.random()
                print(
                    f"[retry {att + 1}/{self.retries}] {reason}; sleeping {delay:.1f}s"
                )
                self.log_debug(f"Retrying after {delay:.1f}s due to {reason}")
                time.sleep(delay)
        logger.error("All attempts exhausted; final failure")
        raise RuntimeError("exhausted retries")

    # ------------------------------------------------------------------
    def ping(self) -> bool:
        """Test connectivity to the Argo service.

        Returns:
            True if service is responding, False otherwise
        """
        try:
            msg = self.chat("Say: Ready to work!")
            return "ready" in msg.lower()
        except Exception:
            return False

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def _extract_txt(self, r: httpx.Response) -> str:
        """Return cleaned text from any successful 200 response."""
        try:
            j = r.json()
        except Exception:
            # If not JSON but body present, treat whole text as the answer
            raw = r.text.strip()
            return raw

        if "choices" in j:
            return j["choices"][0]["message"]["content"].strip()
        for k in ("response", "content", "text"):
            if isinstance(j.get(k), str):
                return j[k].strip()
        # As a final fallback, return the entire trimmed body
        return r.text.strip()

    def _poll_for_result(self, first: httpx.Response):
        """Poll the gateway when we get 102/202 until we receive the final 200."""
        base = self._base_url_fn(self.env)
        poll_url = self._extract_job_url(first, base) or first.url
        self.log_debug(f"Start polling at {poll_url}")
        waited = 0.0
        while waited < self.timeout:
            time.sleep(POLL_EVERY)
            waited += POLL_EVERY
            try:
                r = self.cli.get(poll_url, headers=self.headers)
            except httpx.TimeoutException:
                self.log_debug(f"Timeout while polling (waited {waited:.1f}s)")
                continue
            if r.status_code in _PROCESSING:
                self.log_debug(
                    f"Still processing ({r.status_code}) after {waited:.1f}s"
                )
                continue  # still running
            if r.status_code == 200:
                txt = self._extract_txt(r)
                if txt:
                    self.log_debug(f"Polling succeeded in {waited:.1f}s")
                    return txt
                return None  # blank – treat as failure
            # any other status ⇒ break & let outer loop retry/backoff
            self.log_warning(
                f"Unexpected status {r.status_code} during polling after {waited:.1f}s"
            )
            break
        self.log_warning(f"Polling timed out after {waited:.1f}s")
        return None


# ─────────────────────────────────────────────────────────────────────
def _parse(raw: str) -> tuple[str, str]:
    if not raw:
        return "Unknown", ""
    # Remove leading/trailing whitespace and trailing period
    raw = raw.strip().rstrip(".")

    # Debug raw response
    logging.debug(f"Raw LLM response: {raw}")

    # Split on dash separators
    parts = _DASH.split(raw, maxsplit=1)

    # Clean and extract raw label (strip quotes/punctuation)
    label_part = parts[0].replace("*", "").strip()
    label_raw = label_part.split()[0].strip("\"'")
    # Normalize and apply alias
    label_norm = label_raw.capitalize()
    label_mapped = _ALIAS.get(label_norm, label_norm)
    # Case-insensitive match against known labels
    valid_map = {lbl.lower(): lbl for lbl in LABEL_SET}
    if label_mapped.lower() in valid_map:
        label = valid_map[label_mapped.lower()]
    else:
        logging.warning(f"Unknown label: {label_mapped} from response: {raw}")
        return "Unknown", raw

    # Extract evidence text if present
    note = parts[1].strip() if len(parts) > 1 else ""

    return label, note


def llm_label(raw_response: str) -> tuple[str, str]:
    """Parse LLM response to extract annotation comparison label and justification.

    Args:
        raw_response: Raw text response from LLM

    Returns:
        Tuple of (label, justification) where label is one of the valid annotation labels
    """
    return _parse(raw_response)
