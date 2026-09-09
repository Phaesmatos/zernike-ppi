"""Metrics-oriented helpers for channel ablation experiments."""
import csv, json, time
from pathlib import Path
from .evaluation import evaluate
def write_summary(rows, output):
    p=Path(output); p.mkdir(parents=True,exist_ok=True)
    if rows:
        with (p/"summary.csv").open("w",newline="") as f: w=csv.DictWriter(f,fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)
    (p/"experiments.json").write_text(json.dumps(rows,indent=2)); return p/"summary.csv"
