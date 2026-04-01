from setuptools import setup, find_packages
from pathlib import Path

this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text(encoding="utf-8")

setup(
    name="prsedm",
    version="1.3.0",
    description=(
        "PRSedm (Polygenic Risk Score Extension for Diabetes Mellitus) "
        "is a package for local and remote generation of "
        "Polygenic Risk Scores (PRS) for Diabetes Mellitus (DM)"
    ),
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Seth A. Sharp",
    author_email="ssharp@broadinstitute.org",
    packages=find_packages(include=["prsedm", "prsedm.*"]),
    entry_points={
        "console_scripts": [
            "prsedm = prsedm.extensions.score_dm:main",
        ],
    },
    install_requires=[
        "joblib>=1.3",
        "numpy<2",
        "pandas>=2.2",
        "pysam>=0.22",
    ],
    python_requires=">=3.9",
)