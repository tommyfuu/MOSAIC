"""Clinical key terms for the PubTator literature check, read from the `pubtator_key_terms` section of each
dataset's config YAML, and the "strongly relevant" rule.

A dataset needs at least five key terms. With exactly five, a taxon is strongly relevant if all five terms return
>=1 publication or at least three return >=10. With more than five, only the second criterion is used, because
requiring every term to return a publication becomes arbitrary as the list grows.

Usage from the shell: python3 key_terms.py <dataset> [--config <yml>]   (prints one term per line)
"""
import glob
import os
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MIN_TERMS = 5


def _read_yaml(path):
    try:
        import yaml
        with open(path) as f:
            return yaml.safe_load(f) or {}
    except ImportError:  # the openpyxl env has no PyYAML: read only the two keys used here
        out, key = {}, None
        for line in open(path):
            s = line.split('#', 1)[0].rstrip() if not line.lstrip().startswith('- ') else line.rstrip()
            if not s.strip():
                continue
            if s.startswith('dataset_name:'):
                out['dataset_name'] = s.split(':', 1)[1].strip().strip('\'"')
                key = None
            elif s.startswith('pubtator_key_terms:'):
                out['pubtator_key_terms'], key = [], 'pubtator_key_terms'
            elif key and s.lstrip().startswith('- '):
                out[key].append(s.lstrip()[2:].strip().strip('\'"'))
            elif not s.startswith(' '):
                key = None
        return out


def key_terms(dataset, config=None):
    """Key terms for `dataset`: from `config` if given, otherwise from the repository-root config YAML whose
    dataset_name is `dataset` and which has a pubtator_key_terms section (exactly one must match)."""
    paths = [config] if config else sorted(glob.glob(os.path.join(REPO_ROOT, 'config*.yml')))
    hits = []
    for p in paths:
        d = _read_yaml(p)
        if d.get('dataset_name') == dataset and d.get('pubtator_key_terms'):
            hits.append((p, [str(t) for t in d['pubtator_key_terms']]))
    if len(hits) != 1:
        found = ', '.join(os.path.basename(p) for p, _ in hits) or 'none'
        raise SystemExit(f'{dataset}: expected one config with pubtator_key_terms, found {found}')
    path, terms = hits[0]
    if len(terms) < MIN_TERMS:
        raise SystemExit(f'{dataset}: {os.path.basename(path)} lists {len(terms)} key terms; at least {MIN_TERMS} are needed')
    return terms


def slug(term):
    """File-name form of a term, as used by run_pubtator_keyterms.sh: {dataset}__{slug}.csv."""
    return ''.join(ch for ch in term.replace(' ', '_') if ch.isalnum() or ch in '_-')


def strongly_relevant(n_terms_ge1, n_terms_ge10, k):
    """Default rule for k key terms (k >= 5); works on numbers or pandas Series."""
    if k == MIN_TERMS:
        return (n_terms_ge1 == k) | (n_terms_ge10 >= 3)
    return n_terms_ge10 >= 3


def rule_text(k):
    if k == MIN_TERMS:
        return f'all {k} terms >=1 publication, or >=3 terms >=10'
    return '>=3 terms >=10 publications'


if __name__ == '__main__':
    args = sys.argv[1:]
    cfg = None
    if '--config' in args:
        i = args.index('--config')
        cfg = args[i + 1]
        del args[i:i + 2]
    if len(args) != 1:
        raise SystemExit(__doc__)
    print('\n'.join(key_terms(args[0], cfg)))
