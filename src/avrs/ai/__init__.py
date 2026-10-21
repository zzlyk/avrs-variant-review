"""AI review: schema, versioned prompts, providers, fusion engine."""
from .schema import AIReviewResult, AIRecommendation, ReviewStatus  # noqa: F401
from .providers import (  # noqa: F401
    AIProviderError, MockAIProvider, OpenAICompatibleProvider, OllamaProvider,
    make_provider, encode_image,
)
from .engine import ReviewEngine  # noqa: F401
from .prompts import PROMPT_VERSION, build_review_prompt  # noqa: F401
