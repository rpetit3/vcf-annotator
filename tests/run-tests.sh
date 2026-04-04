#!/bin/bash
set -euo pipefail

PASS=0
FAIL=0

for name in example example-multi example-version; do
    echo "Testing ${name}..."
    python bin/vcf-annotator.py "tests/data/${name}.vcf" "tests/data/${name}.gb" --output "test-${name}.vcf"
    EXPECTED=$(sha256sum "tests/data/${name}.annotated.vcf" | awk '{print $1}')
    ACTUAL=$(sha256sum "test-${name}.vcf" | awk '{print $1}')
    if [ "$EXPECTED" = "$ACTUAL" ]; then
        echo "  PASS (${EXPECTED})"
        PASS=$((PASS + 1))
    else
        echo "  FAIL"
        echo "    Expected: ${EXPECTED}"
        echo "    Actual:   ${ACTUAL}"
        FAIL=$((FAIL + 1))
    fi
    rm -f "test-${name}.vcf"
done

echo ""
echo "${PASS} passed, ${FAIL} failed"
[ "$FAIL" -eq 0 ]
