from pathlib import Path
from pydantic_settings import BaseSettings, SettingsConfigDict

_ENV_FILE = Path(__file__).parent / ".env"


class Settings(BaseSettings):
    model_config = SettingsConfigDict(env_file=str(_ENV_FILE), extra="ignore")

    llm_api_key: str = ""
    llm_base_url: str = "https://api.openai.com/v1"
    llm_model: str = "gpt-4o-mini"
    r_scripts_dir: str = "../r_scripts"
    jobs_dir: str = "./jobs"
    log_level: str = "INFO"
    rscript_path: str = "Rscript"  # override in .env if not on PATH


settings = Settings()
