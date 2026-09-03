#!/bin/bash
S=/private/tmp/claude-501/-Users-Kitty-git-depthSV/1db72424-c4a1-45f4-b81e-5b466f3ec8a4/scratchpad/review_shell
cd "$S" || exit 1
du -sh exwork
for x in exA exB exC exD exE exF exG; do rm -rf "$x"; cp -a exwork "$x"; done
# fake sbatch shim for X4
mkdir -p "$S/bin"
cat > "$S/bin/sbatch" <<'EOF'
#!/bin/bash
# Fake sbatch: record argv and the EX_/DSV_/SLURM_ environment a real job would inherit; emit a federated --parsable id.
n=$(ls "$SBATCH_SHIM_DIR" 2>/dev/null | wc -l | tr -d ' ')
rec="$SBATCH_SHIM_DIR/call.$n"
printf '%s\n' "$@" > "$rec.args"
env | grep -E '^(EX_|DSV_|SLURM_|SBATCH_)' | sort > "$rec.env"
echo "$((1000 + n));cluster1"
EOF
chmod +x "$S/bin/sbatch"
echo "setup done"; ls -d ex?
