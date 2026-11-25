#!/usr/bin/env python3
"""Copy schema artifacts emitted by the C++ authority; no Python policy logic."""
import argparse
import json
from pathlib import Path
import subprocess
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--cli',type=Path,required=True);p.add_argument('--output',type=Path,default=Path('schema'))
a=p.parse_args()
schemas=json.loads(subprocess.check_output([str(a.cli),'--schemas'],text=True))
a.output.mkdir(parents=True,exist_ok=True)
for name,schema in schemas.items():
    (a.output/f'{name}.schema.json').write_text(json.dumps(schema,indent=2,allow_nan=False)+'\n')
print(f'exported {len(schemas)} native schemas')
