---
name: biomedical-research-team
description: Orchestrates a biomedical analysis-and-writing subagent team for publication-grade manuscripts. Use when building evidence-anchored workflows from project results and drafting CMH-style papers with role-separated analysis, methods, writing, and QA.
---

# Biomedical Research Team Cursor Skill (Reusable Template)

## 1) Purpose

This skill defines a reusable, leader-driven multi-agent workflow for:
- biomedical data analysis,
- mechanism-focused interpretation,
- manuscript writing,
- and submission-quality checks.

It is designed to be project-agnostic by replacing only the project root name.

---

## 2) Project Variable System

Use these variables in every path and instruction:

- `{{PROJECT_ROOT}}`: absolute root path of current project
- `{{PROJECT_NAME}}`: root directory name of current project

Example mapping:
- ITLAS project:
  - `{{PROJECT_ROOT}} = G:\내 드라이브\ITLAS`
  - `{{PROJECT_NAME}} = ITLAS`
- Another project:
  - `{{PROJECT_ROOT}} = D:\Research\LiverAtlasX`
  - `{{PROJECT_NAME}} = LiverAtlasX`

Only the root name/path changes; folder structure and workflow remain the same.

---

## 3) Fixed Folder Blueprint

Use this structure in every project:

```text
{{PROJECT_ROOT}}/
├── results/                        # evidence source (analysis inputs)
└── paper-CR/                       # manuscript workspace
    ├── reference/                  # references, literature notes, checklists
    └── images/                     # figures, graphical abstract, visual assets
```

Storage policy:
- Main manuscript and plan files -> `{{PROJECT_ROOT}}/paper-CR`
- Reference/literature files -> `{{PROJECT_ROOT}}/paper-CR/reference`
- Figure/image files -> `{{PROJECT_ROOT}}/paper-CR/images`

---

## 4) Evidence Scope Lock

Hard rule:
- Evidence extraction and quantitative support must come from:
  - `{{PROJECT_ROOT}}/results/**`

No mechanistic conclusion is allowed without traceable evidence artifacts.
Uncertain interpretations must be explicitly marked as hypothesis.

---

## 5) Team Architecture

Leader (main agent) orchestrates all stages and quality gates.

Subagent roles:
1. `analysis-stat`
   - develops analysis algorithms
   - runs independent statistical validation
   - outputs effect sizes/confidence and evidence lock
2. `methods-writer`
   - writes reproducible Methods from validated analysis-stat outputs
3. `writer-introduction`
   - writes concise question-focused Introduction
4. `writer-results-discussion`
   - writes Results/Discussion from validated outputs only
5. `reference-manager`
   - builds citation block, manages literature notes
6. `figure-table`
   - creates figures/tables/captions with format checks
7. `proofreader`
   - language/consistency polishing without changing evidence meaning
8. `fig-tab-matcher`
   - text vs figure/table numbering/order/legend consistency check
9. `reviewer`
   - internal peer-style review (logic, novelty, reproducibility, compliance risk)

---

## 6) Execution Order (Leader Workflow)

1. Build evidence index from `results`.
2. Run `analysis-stat` first and lock validated outputs.
3. Run `methods-writer` from locked outputs.
4. Run `writer-introduction` + `writer-results-discussion`.
5. Run `reference-manager` + `figure-table`.
6. Integrate manuscript in `paper-CR` root.
7. Quality gates:
   - `proofreader` -> `fig-tab-matcher` -> `reviewer`
8. Repeat until critical issues are cleared.

---

## 7) CMH-Oriented Output Constraints (Default)

Target template: CMH Original Article style.

Default checks:
- structured abstract (Background/Aims, Methods, Results, Conclusions)
- abstract <= 250 words, keywords 3-5
- highlights (3-4 sentences, <=100 words)
- graphical abstract prepared
- <= 6,000 words total
- <= 8 figures+tables combined
- <= 50 references
- numeric superscript citation order
- AI-use disclosure in Acknowledgement when applicable

If target journal changes, keep workflow intact and swap only journal-specific format constraints.

---

## 8) Reusability Rule (Most Important)

When porting this skill to a new project:
- Keep all role definitions and workflow unchanged.
- Keep `results -> paper-CR/reference/images` structure unchanged.
- Replace only `{{PROJECT_ROOT}}` / `{{PROJECT_NAME}}` values.

This guarantees cross-project reuse without redesigning the team.

---

## 9) Quick Start Prompt (Reusable)

Use this starter instruction in any project:

`Use the biomedical research team workflow with {{PROJECT_ROOT}} as root. Restrict evidence extraction to {{PROJECT_ROOT}}/results, write outputs to {{PROJECT_ROOT}}/paper-CR (reference/images policy enforced), run analysis-stat independently before methods writing, and produce a CMH-style manuscript draft with full quality gates.`

