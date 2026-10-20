#!/usr/bin/env python3
"""Stage 301 protocol amendment, registered BEFORE any timed run.

Two changes to the screens relative to the registered protocol, both of which
only ENLARGE POS/NEG (i.e. strengthen the necessary conditions) and are proved
in the stage301_common docstring:

  (1) the tolerance band is sharpened from tau = sum_u min(eps/l, 1) to the
      one-sided bands tau_up / tau_dn built from min(eps/l, (l -+ dh)/(2l)),
      dh = h(u) - h(v).  The old band is the special case dh = -l.
  (2) the covering requirement uses POS \\ Sigma rather than POS \\ N[Sigma],
      because for v outside Sigma every Sigma-neighbour has phi > eps and can
      only raise m(v).  This matches the frozen |positive| <= 3k - 2 count.

Also registered: the closed-form certified lower bound
  Theta_k >= max( (4k-7)-th largest critical eps of NEG,
                  (3k-1)-th largest critical eps of POS ),
which needs no support enumeration at all.
"""
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import stage301_common as s301   # noqa: E402

OUT = s301.OUT
PATH = OUT / 'stage301_protocol_amendment.json'


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    assert not PATH.exists(), 'Never overwrite a registered amendment'
    base = json.loads((OUT / 'stage301_protocol.json').read_text())
    payload = {
        'stage': 301,
        'amends': 'stage301_protocol.json',
        'created_utc': datetime.now(timezone.utc).isoformat(),
        'base_code_sha256_at_registration': base['code_sha256']['stage301_common.py'],
        'amended_code_sha256': {n: sha(ROOT / 'scripts' / n)
                                for n in ['stage301_common.py', 'stage301_validate.py',
                                          'stage301_protocol_amendment.py']},
        'changes': ['sharpened one-sided tolerance bands tau_up / tau_dn',
                    'covering requirement POS \\ Sigma instead of POS \\ N[Sigma]',
                    'closed-form refutation threshold Theta_k >= eps_R'],
        'direction': 'both changes enlarge POS/NEG, i.e. strengthen necessary '
                     'conditions; neither can admit a support the weaker screen '
                     'excluded, and both are proved before use',
        'timed_runs_before_amendment': 'none',
    }
    PATH.write_text(json.dumps(payload, indent=1) + '\n')
    print(json.dumps(payload, indent=1)[:900])


if __name__ == '__main__':
    main()
