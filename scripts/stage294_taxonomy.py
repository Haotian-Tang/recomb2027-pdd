#!/usr/bin/env python3
"""Stage 294: provenance and match report for the Web of Life taxonomy.

The stage-287 bank stores WoL-style genome accessions as leaf labels.  The Web
of Life project publishes a per-genome ranked NCBI taxonomy table; this script
registers exactly which file was fetched, from where, with which checksums, and
how completely it matches the accessions the noisy panel actually uses.  It
writes one new manifest and reads nothing read-write.
"""
import collections
import hashlib
import json
import lzma
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import stage294_common as s294                      # noqa: E402

RAW = s294.TAXDIR / 'raw'
MANIFEST = s294.TAXDIR / 'stage294_taxonomy_manifest.json'
REPO = 'biocore/wol'
FILES = {
    'taxonomy_ncbi_ranks.tsv.xz': 'data/taxonomy/ncbi/ranks.tsv.xz',
    'taxonomy_ncbi_lineages.txt.xz': 'data/taxonomy/ncbi/lineages.txt.xz',
    'taxonomy_gtdb_lineages.txt.xz': 'data/taxonomy/gtdb/lineages.txt.xz',
    'genomes_metadata.tsv.xz': 'data/genomes/metadata.tsv.xz',
}


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def blob_sha(remote):
    """The git blob id GitHub reports for the path, recorded as an independent id."""
    url = f'https://api.github.com/repos/{REPO}/contents/{remote}'
    out = subprocess.run(['curl', '-sS', '-m', '30', url], capture_output=True, text=True).stdout
    try:
        return json.loads(out)['sha']
    except Exception:
        return None


def rows(path):
    with lzma.open(path, 'rt') as handle:
        return sum(1 for _ in handle)


def main():
    assert not MANIFEST.exists(), 'Never overwrite a registered manifest'
    files = {}
    for local, remote in FILES.items():
        path = RAW / local
        files[local] = dict(
            url=f'https://raw.githubusercontent.com/{REPO}/master/{remote}',
            repo_path=remote, bytes=path.stat().st_size, sha256=sha256(path),
            github_blob_sha=blob_sha(remote), rows=rows(path))
    ranks = s294.load_ranks()
    labels, observations, matched = set(), 0, 0
    rank_breadth = collections.defaultdict(list)
    blank = collections.Counter()
    for pair in s294.pair_panel():
        case = s294.load_noisy_case(pair['noisy'])
        leaves = sorted(case['observed'])
        observations += len(leaves)
        for leaf in leaves:
            accession = case['leaf_labels'][leaf]
            labels.add(accession)
            matched += accession in ranks
        for rank in ('kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'):
            groups = collections.Counter(ranks[case['leaf_labels'][leaf]][rank].strip()
                                         or s294.UNASSIGNED for leaf in leaves)
            rank_breadth[rank].append(len(groups))
            blank[rank] += sum(v for k, v in groups.items() if k == s294.UNASSIGNED)
    breadth = {}
    for rank, counts in rank_breadth.items():
        counts = sorted(counts)
        breadth[rank] = dict(min=counts[0], median=counts[len(counts) // 2], max=counts[-1],
                             cases_with_one_group=sum(1 for c in counts if c == 1),
                             cases_with_at_least_two=sum(1 for c in counts if c >= 2))
    manifest = dict(
        stage=294, registered_utc=datetime.now(timezone.utc).isoformat(),
        source=dict(project='Web of Life (WoL)', repository=f'https://github.com/{REPO}',
                    reference='Zhu et al. 2019, Nat Commun 10:5477',
                    repository_head_commit='ff92d2251c0ea72299a2ead1aecfd748a9acd698',
                    repository_head_commit_date='2023-02-16T20:28:37Z',
                    fetched_utc=datetime.now(timezone.utc).isoformat()),
        files=files,
        taxonomy_used=dict(file='taxonomy_ncbi_ranks.tsv.xz',
                           columns=['genome', 'kingdom', 'phylum', 'class', 'order', 'family',
                                    'genus', 'species'],
                           genomes=len(ranks)),
        match=dict(distinct_accessions_in_noisy_panel=len(labels),
                   matched_in_ranks_table=len(labels & set(ranks)),
                   observations=observations, matched_observations=matched,
                   match_fraction=matched / observations),
        panel_taxonomic_breadth=breadth,
        unassigned_rank_values_by_rank=dict(blank),
        operable_ranks=dict(
            rule='a rank is operable for a case when its panel spans at least two groups there',
            operable_cases_by_rank={rank: breadth[rank]['cases_with_at_least_two']
                                    for rank in breadth},
            conclusion=('by that rule species is operable in all 296 cases, genus in 192, family '
                        'and order in 134 each, phylum and class in 28 each, and kingdom in none. '
                        'Stage 294 drops whole groups at the FINEST operable rank -- genus where '
                        'the panel spans more than one genus, species elsewhere -- so every case '
                        'receives a taxonomic draw. The coarser ranks are operable in a minority '
                        'of cases (134 of 296 at family and order, 28 at phylum) and are left '
                        'unused by design, not by impossibility: a regime defined only where the '
                        'rank is operable would change the case set between levels and regimes '
                        'and break the draw-for-draw pairing with stage 288. Panel breadth is the '
                        'limiting fact either way -- 162 of 296 panels span a single family and '
                        '268 of 296 a single phylum.'),
            used_by_stage294=list(s294.TAXON_RANKS)))
    MANIFEST.write_text(json.dumps(manifest, indent=1) + '\n')
    print(json.dumps(dict(match=manifest['match'], breadth=breadth), indent=1))


if __name__ == '__main__':
    main()
