"""lpvid 源码可编辑安装元数据。"""

import json
from pathlib import Path

from setuptools import find_packages, setup

ROOT = Path(__file__).resolve().parent
with (ROOT / "config" / "project.json").open(encoding="utf-8") as handle:
    project = json.load(handle)

setup(
    name=project["distribution"],
    version=project["version"],
    description=project["description"],
    long_description=(ROOT / "README.md").read_text(encoding="utf-8"),
    long_description_content_type="text/markdown",
    python_requires=">=3.10",
    packages=find_packages(include=["cdc_pipeline", "cdc_pipeline.*"]),
    include_package_data=True,
    install_requires=[
        "PyYAML==6.0.3",
        "matplotlib==3.10.9",
        "ete3==3.1.3",
        "numpy==2.2.6",
        "pandas==2.3.3",
        "biopython==1.88",
        "tqdm==4.70.1",
    ],
    extras_require={
        "dev": ["pytest==9.2.0"],
    },
    entry_points={
        "console_scripts": [
            f"{project['command']}=cdc_pipeline.cli:main",
        ],
    },
)
