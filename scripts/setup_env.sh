#!/usr/bin/env bash
# setup_env.sh — first-day environment check for new contributors.
#
# Validates LPC tool paths, checks group permissions, creates a Python venv,
# installs project dependencies. Idempotent — safe to re-run.
#
# Usage: bash scripts/setup_env.sh

set -uo pipefail

PROJECT_ROOT="/project/hall/analysis/hearing-loss-genomics"
cd "$PROJECT_ROOT" || { echo "ERROR: cannot cd to $PROJECT_ROOT"; exit 1; }

GREEN='\033[0;32m'; YELLOW='\033[0;33m'; RED='\033[0;31m'; NC='\033[0m'
OK="${GREEN}✓${NC}"; WARN="${YELLOW}⚠${NC}"; FAIL="${RED}✗${NC}"

FAILURES=0
ok()   { echo -e "  $OK $1"; }
warn() { echo -e "  $WARN $1"; }
fail() { echo -e "  $FAIL $1"; FAILURES=$((FAILURES+1)); }

echo "=================================================================="
echo " Hearing Loss Genomics — environment setup + validation"
echo "=================================================================="
echo "Project root: $PROJECT_ROOT"
echo "Hostname:     $(hostname)"
echo "User:         $(whoami)"
echo

# ------------------------------------------------------------------
# 1. Group memberships
# ------------------------------------------------------------------
echo "[1/6] Checking LPC group memberships..."
if groups | tr ' ' '\n' | grep -qx hall; then
    ok "Member of 'hall' group (project files)"
else
    fail "Not in 'hall' group — request access from Molly or LPC IT"
fi
if groups | tr ' ' '\n' | grep -qx ritchie; then
    ok "Member of 'ritchie' group (raw PMBB data)"
else
    fail "Not in 'ritchie' group — request access (needed for /static/PMBB/...)"
fi
echo

# ------------------------------------------------------------------
# 2. Critical filesystem paths
# ------------------------------------------------------------------
echo "[2/6] Checking critical filesystem paths..."
PATHS=(
    "/static/PMBB/PMBB-Release-2020-2.0:Raw PMBB v2 (read access expected)"
    "$PROJECT_ROOT/data/PMBB_Exome:Daniel's preserved v2 intermediates"
    "$PROJECT_ROOT/data/pmbb_v2:Symlink to raw PMBB v2"
    "$PROJECT_ROOT/analysis/daniel/runbook_hui2023.txt:Daniel's decompressed cookbook"
)
for entry in "${PATHS[@]}"; do
    p="${entry%%:*}"
    desc="${entry##*:}"
    if [[ -e "$p" ]]; then
        ok "$p exists ($desc)"
    else
        fail "$p MISSING ($desc)"
    fi
done
echo

# ------------------------------------------------------------------
# 3. LPC tools (versions matter)
# ------------------------------------------------------------------
echo "[3/6] Checking LPC tool availability..."

PLINK19="/appl/plink-1.90Beta6.18/plink"
if [[ -x "$PLINK19" ]]; then
    ok "plink 1.9 at $PLINK19"
else
    fail "plink 1.9 NOT FOUND at $PLINK19 — may need to find new path"
fi

R44="/appl/R-4.4/bin/Rscript"
if [[ -x "$R44" ]]; then
    ok "R-4.4 at $R44 (default R-3.6.3 is broken — use this explicitly)"
else
    fail "R-4.4 NOT FOUND at $R44"
fi

LOKI="/project/ritchie/datasets/loki/loki-20230816.db"
if [[ -f "$LOKI" ]]; then
    sz=$(du -h "$LOKI" | cut -f1)
    ok "LOKI database at $LOKI ($sz)"
else
    fail "LOKI database NOT FOUND at $LOKI — Phase 5+ will fail"
fi

BIOBIN_SHIM="$PROJECT_ROOT/analysis/daniel/configs/lib-shims/liblzma.so.0"
if [[ -L "$BIOBIN_SHIM" ]]; then
    target=$(readlink "$BIOBIN_SHIM")
    ok "biobin liblzma shim at $BIOBIN_SHIM → $target"
else
    warn "biobin liblzma shim missing — Phase 5+ will need it for RHEL 9"
    warn "  Create with: ln -s /usr/lib64/liblzma.so.5 $BIOBIN_SHIM"
fi

# biobin itself - try a couple common locations
BIOBIN_PATH=""
for cand in "/project/ritchie/env/modules/rlsoftware/latest/bin/biobin" "$(command -v biobin 2>/dev/null)"; do
    if [[ -n "$cand" && -x "$cand" ]]; then
        BIOBIN_PATH="$cand"
        break
    fi
done
if [[ -n "$BIOBIN_PATH" ]]; then
    ok "biobin found at $BIOBIN_PATH"
else
    warn "biobin not found in usual locations — check with Ritchie Lab module"
fi

LSF_BSUB=$(command -v bsub 2>/dev/null)
if [[ -n "$LSF_BSUB" ]]; then
    ok "LSF bsub at $LSF_BSUB (job submission available)"
else
    warn "LSF bsub NOT FOUND — local-only mode; heavy phases need LSF"
fi
echo

# ------------------------------------------------------------------
# 4. Python venv
# ------------------------------------------------------------------
echo "[4/6] Checking Python venv..."

VENV_DIR="$PROJECT_ROOT/venv"
PY312="/appl/python-3.12/bin/python3"

if [[ -d "$VENV_DIR" ]]; then
    OWNER=$(stat -c '%U' "$VENV_DIR")
    if [[ "$OWNER" == "$(whoami)" ]]; then
        ok "venv exists and belongs to you"
    else
        warn "venv exists but is owned by '$OWNER' — you likely need your own"
        warn "  Move/rename it: mv venv venv.${OWNER} && create a new one (below)"
    fi
else
    echo "  venv/ does not exist. Creating..."
    if [[ -x "$PY312" ]]; then
        "$PY312" -m venv "$VENV_DIR" && ok "venv created with $PY312"
    else
        warn "Python 3.12 at $PY312 not found — falling back to system python3"
        python3 -m venv "$VENV_DIR" && ok "venv created with system python3"
    fi
fi

if [[ -x "$VENV_DIR/bin/python3" ]]; then
    VENV_PY=$("$VENV_DIR/bin/python3" --version 2>&1)
    ok "venv Python: $VENV_PY"
fi
echo

# ------------------------------------------------------------------
# 5. Install Python dependencies
# ------------------------------------------------------------------
echo "[5/6] Installing Python dependencies from requirements.txt..."

if [[ -f "$VENV_DIR/bin/activate" ]]; then
    # shellcheck disable=SC1091
    source "$VENV_DIR/bin/activate"
    if [[ -f "$PROJECT_ROOT/requirements.txt" ]]; then
        pip install --quiet --upgrade pip
        if pip install --quiet -r "$PROJECT_ROOT/requirements.txt"; then
            n=$(pip list 2>/dev/null | wc -l)
            ok "Installed/verified dependencies ($n packages)"
        else
            fail "pip install failed — check network or pip version"
        fi
    else
        warn "requirements.txt missing"
    fi
    deactivate
else
    fail "venv activation script missing — venv setup failed"
fi
echo

# ------------------------------------------------------------------
# 6. Sanity check — light-mode Phase 1 reproduces a known intermediate
# ------------------------------------------------------------------
echo "[6/6] Sanity check — reading a known Daniel intermediate..."

SANITY_FILE="$PROJECT_ROOT/data/PMBB_Exome/all_genes_including_ShadisList.txt.gz"
if [[ -f "$SANITY_FILE" ]]; then
    n_genes=$(zcat "$SANITY_FILE" 2>/dev/null | wc -l)
    if [[ "$n_genes" -gt 0 ]]; then
        ok "Read Daniel's HL gene list ($n_genes genes)"
    else
        fail "HL gene list empty or unreadable"
    fi
else
    fail "HL gene list NOT FOUND at $SANITY_FILE"
fi
echo

# ------------------------------------------------------------------
# Summary
# ------------------------------------------------------------------
echo "=================================================================="
if [[ "$FAILURES" -eq 0 ]]; then
    echo -e "${GREEN}All critical checks passed.${NC}"
    echo
    echo "Next steps:"
    echo "  1. source venv/bin/activate"
    echo "  2. Read HANDOFF.md (you should have already)"
    echo "  3. Read REPRODUCTION_GUIDE.md to re-run any phase"
    echo "  4. Read STATUS_SNAPSHOT.md to see where we left off"
else
    echo -e "${RED}$FAILURES critical check(s) failed.${NC}"
    echo "Resolve before running any phase. Re-run this script after fixes."
    exit 1
fi
echo "=================================================================="
