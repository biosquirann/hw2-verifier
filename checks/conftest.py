import os, sys, pathlib
from hypothesis import settings
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))
# default: детерминированный прогон (derandomize) 
# explore: случайные входы, seed задатся
settings.register_profile("default", max_examples=200, derandomize=True, deadline=None)
settings.register_profile("explore", max_examples=300, deadline=None)
settings.load_profile(os.environ.get("HYPOTHESIS_PROFILE", "default"))
