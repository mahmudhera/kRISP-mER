"""
Rule Set 2 (Doench 2016) on-target scoring via Azimuth 2.0.

Azimuth 2.0 is vendored at vendor/Azimuth-2.0, so no separate
download is needed. It is BSD-3-Clause, copyright 2015 Microsoft Research;
see vendor/Azimuth-2.0/LICENSE.txt. The copy came from CRISPOR, which
carries a Python-3 re-pickled model alongside the original.

Why not krispmer/on_target_scoring: that code needs V3_model_nopos.pickle,
which is a Python-2 protocol-2 pickle referencing PresortBestSplitter, a
Cython class removed from scikit-learn after 0.17. It cannot be loaded on
sklearn 1.1.

Verified 6 Nov 2026: the CRISPOR doctest 30-mer
ACGCAGCGCCTCCCTCCACTCGGAAGGACT scores 0.5513089, i.e. 55 on CRISPOR's
0-100 scale, matching CRISPOR's expected fusi value.

Azimuth's util.py uses Python-2 implicit relative imports, so BOTH the
package directory and its parent must be on sys.path.

Cite: Doench et al., Nature Biotechnology 34:184-191, 2016.
"""

import os
import sys
import logging

# The vendored copy, found relative to this file, comes first. The
# CRISPOR locations stay as fallbacks for installs that predate vendoring.
_VENDORED = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    'vendor', 'Azimuth-2.0')

_AZIMUTH_DIRS = [
    _VENDORED,
    os.path.expanduser('~/crisporWebsite/bin/Azimuth-2.0'),
    os.path.expanduser('~/Crispor/crisporWebsite/bin/Azimuth-2.0'),
]
_AZIMUTH_ENV = os.environ.get('AZIMUTH_DIR')
if _AZIMUTH_ENV:
    _AZIMUTH_DIRS.insert(0, _AZIMUTH_ENV)

_predict = None
_unavailable_reason = None


def _init():
    """Import azimuth once. Returns the predict function, or None."""
    global _predict, _unavailable_reason
    if _predict is not None or _unavailable_reason is not None:
        return _predict

    for base in _AZIMUTH_DIRS:
        inner = os.path.join(base, 'azimuth')
        if not os.path.isdir(inner):
            continue
        for p in (base, inner):
            if p not in sys.path:
                sys.path.append(p)
        try:
            import azimuth.model_comparison as mc
            _predict = mc.predict
            logging.info('Rule Set 2: using Azimuth-2.0 at %s', base)
            return _predict
        except Exception as e:
            _unavailable_reason = '%s: %s' % (type(e).__name__, e)

    if _unavailable_reason is None:
        _unavailable_reason = 'Azimuth-2.0 not found in %s' % _AZIMUTH_DIRS
    logging.warning('Rule Set 2 unavailable (%s). Scores will be -1.',
                    _unavailable_reason)
    return None


def rs2_score(seq30, scale_100=True):
    """Rule Set 2 score for one 30-mer (4 bp + 20 bp guide + 3 bp PAM + 3 bp).

    Returns -1 when the sequence is unusable or Azimuth is unavailable,
    which is the convention CRISPOR uses.
    """
    predict = _init()
    if predict is None or seq30 is None:
        return -1
    seq30 = seq30.upper()
    if len(seq30) != 30 or 'N' in seq30:
        return -1
    try:
        import numpy
        s = float(predict(numpy.array([seq30]), None, None, pam_audit=False)[0])
    except Exception as e:
        logging.warning('Rule Set 2 failed on %s: %s', seq30, e)
        return -1
    return int(round(100 * s)) if scale_100 else s


def rs2_scores(seq30_list, scale_100=True):
    """Score many 30-mers in one call. Same -1 convention per entry."""
    predict = _init()
    if predict is None:
        return [-1] * len(seq30_list)

    import numpy
    out = [-1] * len(seq30_list)
    good_idx, good_seq = [], []
    for i, s in enumerate(seq30_list):
        if s and len(s) == 30 and 'N' not in s.upper():
            good_idx.append(i)
            good_seq.append(s.upper())
    if not good_seq:
        return out
    try:
        preds = predict(numpy.array(good_seq), None, None, pam_audit=False)
    except Exception as e:
        logging.warning('Rule Set 2 batch failed: %s', e)
        return out
    for i, p in zip(good_idx, preds):
        out[i] = int(round(100 * float(p))) if scale_100 else float(p)
    return out


if __name__ == '__main__':
    ref = 'ACGCAGCGCCTCCCTCCACTCGGAAGGACT'
    got = rs2_score(ref)
    print('reference 30-mer:', ref)
    print('score           :', got)
    print('expected        : 55  (CRISPOR doctest)')
    print('MATCH' if got == 55 else 'MISMATCH')
