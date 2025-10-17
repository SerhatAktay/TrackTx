#!/bin/bash
# Diagnose conda setup on Linux server

echo "=== Conda Installation ==="
which conda || echo "conda not in PATH"
which mamba || echo "mamba not in PATH"
which micromamba || echo "micromamba not in PATH"

if command -v conda >/dev/null 2>&1; then
  echo -e "\n=== Conda Info ==="
  conda info
  
  echo -e "\n=== Conda Environments ==="
  conda env list
fi

echo -e "\n=== Conda Activation Paths ==="
for path in \
  "/usr/bin/activate" \
  "/home/$USER/miniconda3/bin/activate" \
  "/home/$USER/anaconda3/bin/activate" \
  "$HOME/miniconda3/bin/activate" \
  "$HOME/anaconda3/bin/activate"; do
  if [ -f "$path" ]; then
    echo "✓ Found: $path"
  else
    echo "✗ Missing: $path"
  fi
done

echo -e "\n=== Storage Type Check ==="
df -T . 2>/dev/null || df .
mount | grep "$(df . | tail -1 | awk '{print $1}')" || echo "Could not determine mount type"

echo -e "\n=== Write Test ==="
test_file=".conda_write_test_$$"
if touch "$test_file" 2>/dev/null; then
  echo "✓ Write permission OK"
  rm -f "$test_file"
else
  echo "✗ Write permission FAILED"
fi

echo -e "\n=== Environment Variables ==="
echo "NXF_CONDA_ENV_PATH: ${NXF_CONDA_ENV_PATH:-not set}"
echo "CONDA_EXE: ${CONDA_EXE:-not set}"
echo "CONDA_PREFIX: ${CONDA_PREFIX:-not set}"

echo -e "\n=== Recommendations ==="
echo "1. Create prebuilt environment:"
echo "   conda env create -f envs/tracktx.yaml -p /tmp/$USER/tracktx-env"
echo ""
echo "2. Export environment path:"
echo "   export NXF_CONDA_ENV_PATH=\"/tmp/$USER/tracktx-env\""
echo ""
echo "3. Run pipeline:"
echo "   ./run_pipeline.sh -profile conda_server"
