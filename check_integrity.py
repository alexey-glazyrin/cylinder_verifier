"""Check the SHA-256 hashes of immutable release files."""
from pathlib import Path
import hashlib
import json
import sys

root=Path(__file__).resolve().parent
manifest=json.loads((root/'MANIFEST.sha256.json').read_text(encoding='utf-8'))
failures=[]
for name, expected in manifest['files'].items():
    path=root/name
    if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest()!=expected:
        failures.append(name)
if failures:
    print('FAIL: missing or changed files:', *failures, sep='\n')
    sys.exit(1)
print(f'PASS: {len(manifest["files"])} immutable release files match their SHA-256 hashes.')
