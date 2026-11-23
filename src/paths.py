#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""交付包内部路径解析：全部相对包根目录，不依赖本机任何安装位置。"""
import os
import sys
from pathlib import Path

BUNDLE = Path(__file__).resolve().parent.parent
SRC = BUNDLE / "src"
MODELS = BUNDLE / "models"
CONFIG = BUNDLE / "config"
KB_FILE = BUNDLE / "molecular_build_rules.yaml"
CODONBERT_DIR = SRC / "codonbert"
# 所有子进程复用当前解释器；如需指定，设置 BIO_SEQDESIGN_PYTHON
PY_OPT = Path(os.environ.get("BIO_SEQDESIGN_PYTHON") or sys.executable)
