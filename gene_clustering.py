#!/usr/bin/env python3
"""Scaffold-stratified Monte Carlo tests of gene proximity and clustering.

Run --help for input options. See README.md for the statistical definitions,
input assumptions, complete column dictionary, and reproducibility.
"""
import argparse
import csv
import hashlib
import json
import platform
import time
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from urllib.parse import unquote

import numpy as np
import openpyxl

VERSION = '2.0.0'
SPECIES = {'Ag': 'a_gambiae', 'Am': 'a_mellifera', 'Bt': 'b_terrestris',
           'Dm': 'd_melanogaster', 'Fv': 'f_varia', 'Mb': 'm_bicolor',
           'Mq': 'm_quadrifasciata', 'Mr': 'm_rotundata', 'Ta': 't_angustula'}
# Input schema for the study workbook format.
WORKBOOK_SCHEMA = {
    'GST': ('S2-GSTs', {'family': 'familia'}),
    'P450': ('S3-P450', {'subfamily': 'Sub_familia', 'family': 'Familia', 'class': 'Classe'}),
    'CCE': ('S5-CCEs', {'clade': 'clado', 'class': 'Classe'}),
    'UGT': ('S6-UGTs', {'family': 'Familias'}),
    'ABC': ('S7-ABCs', {'subfamily': 'Sub_familia'}),
}
# Fixed level codes used to derive reproducible group-specific random streams.
SEED_LEVEL_CODES = {'family': 'familia', 'subfamily': 'subfamilia', 'clade': 'clado', 'class': 'classe'}
KEYS = ['species', 'superfamily', 'level', 'group']
DISTANCE_METRICS = ['mean_bp', 'q1_bp', 'median_bp', 'q3_bp', 'p90_bp', 'p95_bp',
                    'minimum_bp', 'maximum_bp', 'zero_fraction', 'within_threshold_fraction']
METRICS = DISTANCE_METRICS + ['clustered_fraction']
BASE_FIELDS = KEYS + ['n_genes', 'n_scaffolds', 'n_without_partner', 'n_evaluable',
    'counts_per_scaffold', 'threshold_bp', 'iterations', 'status',
    'p_proximity', 'q_proximity_bh', 'p_clustering', 'q_clustering_bh',
    'observed_n_clusters', 'observed_n_clustered_genes']
SUMMARY_FIELDS = (BASE_FIELDS + ['observed_' + m for m in METRICS] +
    ['null_' + m + '_' + suffix for m in METRICS for suffix in ['mean', 'p025', 'p50', 'p975']] +
    ['observed_to_null_median_ratio', 'null_median_variable', 'null_clustering_variable'])
DETAIL_FIELDS = KEYS + ['gene_id', 'scaffold', 'start', 'end', 'nearest_gene_id',
    'nearest_distance_bp', 'overlap', 'status', 'cluster_id', 'previous_group_gene_gap_bp']
CLUSTER_FIELDS = KEYS + ['cluster_id', 'scaffold', 'n_genes', 'start', 'end', 'span_bp',
    'maximum_link_gap_bp', 'gene_ids']
AUDIT_FIELDS = ['species', 'superfamily', 'input_gene', 'scaffold', 'start', 'end',
                'status', 'gene_id', 'name_candidates', 'coordinate_candidates']

@dataclass
class Gene:
    id: str
    scaffold: str
    start: int
    end: int
    aliases: set = field(default_factory=set)
    multipart: bool = False


def text(value):
    value = '' if value is None else str(value).strip()
    return '' if value.lower() in {'', '-', 'na', 'nan', 'none', 'n/a'} else value


def integer(value):
    x = float(value)
    if not np.isfinite(x) or x < 1 or x != int(x):
        raise ValueError(f'Invalid positive integer coordinate: {value!r}')
    return int(x)


def attrs(value):
    return {k: unquote(v) for k, v in
            (part.split('=', 1) for part in value.split(';') if '=' in part)}


def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def write_csv(path, rows, fields):
    with path.open('w', encoding='utf-8-sig', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def read_gff(path):
    """Read RefSeq-style gene records; keep separate overlapping loci."""
    genes, organelles, multipart_ids = {}, set(), set()
    biotypes, repeated = Counter(), 0
    with path.open(encoding='utf-8') as f:
        for line_no, line in enumerate(f, 1):
            if line.startswith('##FASTA'):
                break
            if line.startswith('#') or not line.strip():
                continue
            c = line.rstrip('\r\n').split('\t')
            if len(c) != 9:
                raise ValueError(f'{path.name}:{line_no}: expected 9 GFF columns.')
            if c[2] not in {'gene', 'region'}:
                continue
            a = attrs(c[8])
            if c[2] == 'region':
                genome = a.get('genome', '').lower()
                if any(s in genome for s in ('mitochond', 'chloroplast', 'plastid')):
                    organelles.add(c[0])
                continue
            biotypes[a.get('gene_biotype', 'MISSING_BIOTYPE')] += 1
            if a.get('gene_biotype') != 'protein_coding' or a.get('pseudo', '').lower() == 'true':
                continue
            gid = a.get('ID', '')
            if not gid:
                raise ValueError(f'{path.name}:{line_no}: gene without ID.')
            start, end = integer(c[3]), integer(c[4])
            if end < start:
                raise ValueError(f'{path.name}:{line_no}: end < start.')
            aliases = {gid, gid.removeprefix('gene-')}
            aliases.update(a[k] for k in ('Name', 'gene', 'locus_tag') if a.get(k))
            marked = bool(a.get('part')) and 'trans-splicing' in a.get('exception', '').lower()
            g = Gene(gid, c[0], start, end, aliases, marked)
            if gid in genes:
                old = genes[gid]
                if (old.scaffold, old.start, old.end) == (g.scaffold, start, end):
                    repeated += 1
                elif old.multipart and marked and old.scaffold == g.scaffold:
                    old.start, old.end = min(old.start, start), max(old.end, end)
                    multipart_ids.add(gid)
                else:
                    raise ValueError(f'{path.name}: repeated ID with conflicting coordinates: {gid}')
                old.aliases.update(aliases)
            else:
                genes[gid] = g
    excluded = sum(g.scaffold in organelles for g in genes.values())
    genes = {k: g for k, g in genes.items() if g.scaffold not in organelles}
    if not genes:
        raise ValueError(f'{path.name}: no eligible genes; check gene_biotype attributes.')
    by_scaffold = defaultdict(list)
    for g in genes.values():
        by_scaffold[g.scaffold].append(g)
    for group in by_scaffold.values():
        group.sort(key=lambda g: (g.start, g.end, g.id))
    audit = dict(file=str(path), sha256=sha256(path), gene_line_biotypes=dict(biotypes),
                 eligible_genes=len(genes), organelle_genes_excluded=excluded,
                 scaffolds=len(by_scaffold), detected_organelle_scaffolds=sorted(organelles),
                 merged_multipart_ids=sorted(multipart_ids), duplicate_lines=repeated)
    return genes, by_scaffold, audit


def workbook_records(path, args):
    """Convert the study workbook into normalized input records."""
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    try:
        for sf, (sheet, levels) in WORKBOOK_SCHEMA.items():
            if args.superfamily and sf != args.superfamily:
                continue
            if args.level:
                levels = {k: v for k, v in levels.items() if k == args.level}
            if not levels:
                continue
            if sheet not in wb.sheetnames:
                raise ValueError(f'Missing workbook sheet: {sheet}')
            rows = wb[sheet].iter_rows(values_only=True)
            headers = [text(x) for x in next(rows)]
            required = {'Species_Code', 'Gene_Name', 'Scaffold', 'Start', 'End', *levels.values()}
            if not required.issubset(headers):
                raise ValueError(f'{sheet}: missing columns: {required - set(headers)}')
            seen = set()
            for row_no, values in enumerate(rows, 2):
                if all(v is None for v in values):
                    continue
                r = dict(zip(headers, values))
                code, name = text(r['Species_Code']), text(r['Gene_Name'])
                if code not in SPECIES:
                    raise ValueError(f'{sheet}:{row_no}: unknown species code: {code!r}')
                if args.species and SPECIES[code] != args.species:
                    continue
                if not name or (code, name) in seen:
                    raise ValueError(f'{sheet}:{row_no}: empty or duplicate Gene_Name: {name}')
                seen.add((code, name))
                for level, column in levels.items():
                    yield dict(species=SPECIES[code], superfamily=sf, level=level,
                        group=text(r[column]), gene=name, scaffold=text(r['Scaffold']),
                        start=r['Start'], end=r['End'], alias=text(r.get('Alias_GFF')),
                        input_location=f'{sheet}:{row_no}')
    finally:
        wb.close()


def csv_records(path):
    with path.open(encoding='utf-8-sig', newline='') as f:
        reader = csv.DictReader(f)
        required = set(KEYS + ['gene', 'scaffold', 'start', 'end'])
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f'CSV missing columns: {required - set(reader.fieldnames or [])}')
        for row_no, row in enumerate(reader, 2):
            if None in row:
                raise ValueError(f'CSV line {row_no}: excess fields; check delimiters/quoting.')
            yield dict(row, input_location=f'CSV:{row_no}')


def read_input(args):
    groups, source, missing = defaultdict(list), {}, []
    membership = {}
    records = workbook_records(args.workbook, args) if args.workbook else csv_records(args.genes)
    for raw in records:
        r = {k: text(v) for k, v in raw.items()}
        for column in ['species', 'superfamily', 'level', 'gene']:
            if not r.get(column):
                raise ValueError(f'{r["input_location"]}: empty {column}.')
        if r['level'] not in SEED_LEVEL_CODES:
            raise ValueError(f'Unsupported level: {r["level"]}; use {list(SEED_LEVEL_CODES)}')
        if any(getattr(args, k) and r[k] != getattr(args, k)
               for k in ['species', 'superfamily', 'level']):
            continue
        if not r.get('group'):
            missing.append({k: r[k] for k in ['species', 'superfamily', 'level', 'gene', 'input_location']})
            continue
        if args.group and r['group'] != args.group:
            continue
        key = tuple(r[k] for k in KEYS)
        ref = (r['species'], r['superfamily'], r['gene'])
        member_key = (ref, r['level'])
        if member_key in membership:
            raise ValueError(f'Duplicate/conflicting group assignment: {member_key}')
        membership[member_key] = r['group']
        start, end = integer(r.get('start')), integer(r.get('end'))
        if end < start or not r.get('scaffold'):
            raise ValueError(f'{r["input_location"]}: invalid scaffold/coordinates.')
        coords = (r['scaffold'], start, end)
        aliases = {r['gene'], r.get('alias', '')} - {''}
        if ref in source and source[ref][:3] != coords:
            raise ValueError(f'Conflicting coordinates across levels: {ref}')
        if ref in source:
            aliases |= source[ref][3]
        source[ref] = (*coords, aliases)
        groups[key].append(ref)
    return groups, source, missing


def match_genes(source, universes):
    matched, audit = {}, []
    for species in sorted(universes):
        alias_index, coord_index = defaultdict(set), defaultdict(set)
        for gid, g in universes[species].items():
            coord_index[(g.scaffold, g.start, g.end)].add(gid)
            for alias in g.aliases:
                alias_index[alias].add(gid)
        for ref, (scaffold, start, end, aliases) in source.items():
            if ref[0] != species:
                continue
            named = set().union(*(alias_index[a] for a in aliases))
            coords = coord_index[(scaffold, start, end)]
            candidates = named & coords if named else coords
            status = ('matched_name_and_coordinates' if named else 'matched_coordinates_only') if len(candidates) == 1 else 'error'
            gid = next(iter(candidates)) if len(candidates) == 1 else ''
            audit.append(dict(species=species, superfamily=ref[1], input_gene=ref[2],
                scaffold=scaffold, start=start, end=end, status=status, gene_id=gid,
                name_candidates=';'.join(sorted(named)), coordinate_candidates=';'.join(sorted(coords))))
            if gid:
                matched[ref] = universes[species][gid]
    return matched, audit


def nearest(coordinates):
    n = len(coordinates)
    if n < 2:
        return np.empty(0, dtype=np.int64), np.empty(0, dtype=np.int64)
    starts, ends = coordinates[:, 0], coordinates[:, 1]
    distances = np.maximum(0, np.maximum(starts[:, None], starts[None, :])
                           - np.minimum(ends[:, None], ends[None, :]) - 1)
    np.fill_diagonal(distances, np.iinfo(np.int64).max)
    neighbors = distances.argmin(axis=1)
    return distances[np.arange(n), neighbors], neighbors


def cluster_blocks(coordinates, threshold):
    """Immediate-predecessor runs, NOT a maximum-end interval sweep."""
    if len(coordinates) < 2:
        return []
    order = np.lexsort((coordinates[:, 1], coordinates[:, 0]))
    coords = coordinates[order]
    gaps = coords[1:, 0] - coords[:-1, 1]
    breaks = np.flatnonzero(gaps > threshold) + 1
    return [part for part in np.split(order, breaks) if len(part) >= 2]


def cluster_counts(coordinates, threshold):
    blocks = cluster_blocks(coordinates, threshold)
    return len(blocks), sum(map(len, blocks))


def stats(distances, threshold):
    q = np.percentile(distances, [25, 50, 75, 90, 95])
    values = [float(np.mean(distances)), *map(float, q), int(min(distances)),
              int(max(distances)), float(np.mean(distances == 0)), float(np.mean(distances <= threshold))]
    return dict(zip(DISTANCE_METRICS, values))


def adjust_bh(rows, p_field, q_field):
    tested = sorted((r for r in rows if r['status'] == 'tested'), key=lambda r: r[p_field])
    q = 1.0
    for i in range(len(tested) - 1, -1, -1):
        q = min(q, tested[i][p_field] * len(tested) / (i + 1))
        tested[i][q_field] = q


def group_rng(seed, key):
    species, sf, level, group = key
    seed_key = [seed, species, sf, SEED_LEVEL_CODES[level], group]
    digest = hashlib.sha256(json.dumps(seed_key).encode()).digest()
    return np.random.default_rng(int.from_bytes(digest[:16], 'big'))


def save_xlsx(path, rows):
    """Numeric workbook for full results, not the potentially huge simulations."""
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = 'Results'
    ws.append(SUMMARY_FIELDS)
    for r in rows:
        ws.append([None if r.get(k, '') == '' else r[k] for k in SUMMARY_FIELDS])
    ws.freeze_panes = 'E2'
    ws.auto_filter.ref = ws.dimensions
    for cell in ws[1]:
        cell.font = openpyxl.styles.Font(bold=True)
    for i, name in enumerate(SUMMARY_FIELDS, 1):
        ws.column_dimensions[openpyxl.utils.get_column_letter(i)].width = min(48, max(16, len(name) + 2))
    wb.save(path)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument('--workbook', type=Path, help='Study workbook; use the sheet and column names specified in README.md.')
    source.add_argument('--genes', type=Path, help='Normalized gene CSV; see examples/genes.csv.')
    parser.add_argument('--gff-dir', type=Path, required=True, help='Directory of SPECIES_proteome_genomic.gff files.')
    parser.add_argument('--output', type=Path, help='New output directory; must not exist.')
    parser.add_argument('--species')
    parser.add_argument('--superfamily')
    parser.add_argument('--level', choices=list(SEED_LEVEL_CODES))
    parser.add_argument('--group')
    parser.add_argument('--iterations', type=int, default=9999)
    parser.add_argument('--seed', type=int, default=20261103)
    parser.add_argument('--threshold', type=int, default=35000, help='bp; default 35000, used for both threshold statistics.')
    parser.add_argument('--no-simulations', action='store_true', help='Do not save per-draw CSV (tests unchanged).')
    parser.add_argument('--xlsx', action='store_true', help='Also export a typed Excel summary.')
    args = parser.parse_args(argv)
    if args.iterations < 1 or args.threshold < 0:
        parser.error('Iterations must be >=1 and threshold must be >=0.')
    args.gff_dir = args.gff_dir.resolve()
    input_path = (args.workbook or args.genes).resolve()
    groups, source, missing = read_input(args)
    if not groups:
        raise ValueError('No groups selected; check inputs and filters.')
    out = (args.output or Path('results') / ('run_' + datetime.now().strftime('%Y%m%d_%H%M%S_%f'))).resolve()
    out.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    metadata = dict(script_version=VERSION, parameters={k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()},
        input_file=str(input_path), input_sha256=sha256(input_path), script_sha256=sha256(Path(__file__)),
        python=platform.python_version(), numpy=np.__version__, openpyxl=openpyxl.__version__,
        completed=False, gffs={}, rng_namespace='fixed classification-level seed codes; SHA256; NumPy default_rng',
        clustering_rule='consecutive start minus previous end, no -1; no maximum-end sweep',
        nearest_rule='max(0, max(starts)-min(ends)-1)',
        bh_scope='all tested groups in this execution, separately for each statistic')
    def save_metadata():
        (out / 'parameters.json').write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding='utf-8')
    save_metadata()
    universes, pools = {}, {}
    for species in sorted({key[0] for key in groups}):
        if Path(species).name != species or '\\' in species:
            raise ValueError(f'Invalid species file key: {species!r}')
        path = args.gff_dir / f'{species}_proteome_genomic.gff'
        genes, by_scaffold, audit = read_gff(path)
        universes[species] = genes
        pools[species] = {s: np.array([(g.start, g.end) for g in gs], dtype=np.int64) for s, gs in by_scaffold.items()}
        metadata['gffs'][species] = audit
        print(f'{species}: {len(genes)} eligible loci.', flush=True)
    save_metadata()
    matched, audit = match_genes(source, universes)
    write_csv(out / 'gene_matching_audit.csv', audit, AUDIT_FIELDS)
    write_csv(out / 'missing_classifications.csv', missing, ['species', 'superfamily', 'level', 'gene', 'input_location'])
    if any(r['status'] == 'error' for r in audit):
        raise ValueError(f'Unresolved gene matches; inspect {out}/gene_matching_audit.csv. No simulations started.')
    summary, details, catalog = [], [], []
    sim_file = None
    try:
        if not args.no_simulations:
            sim_file = (out / 'simulation_statistics.csv').open('w', encoding='utf-8-sig', newline='')
            writer = csv.DictWriter(sim_file, fieldnames=KEYS + ['iteration'] + METRICS + ['n_clusters', 'n_clustered_genes'])
            writer.writeheader()
        for number, key in enumerate(sorted(groups), 1):
            identity = dict(zip(KEYS, key))
            observed = [matched[ref] for ref in groups[key]]
            if len({g.id for g in observed}) != len(observed):
                raise ValueError(f'{key}: different input genes resolve to the same locus.')
            by_scaffold = defaultdict(list)
            for g in observed:
                by_scaffold[g.scaffold].append(g)
            blocks, observed_distances = [], []
            singletons, cluster_number, clustered_count = 0, 0, 0
            for scaffold, genes in sorted(by_scaffold.items()):
                genes.sort(key=lambda g: (g.start, g.end, g.id))
                coords = np.array([(g.start, g.end) for g in genes], dtype=np.int64)
                distances, neighbors = nearest(coords)
                memberships, link_gaps = {}, {}
                for block in cluster_blocks(coords, args.threshold):
                    cluster_number += 1
                    clustered_count += len(block)
                    cid = f'cluster_{cluster_number:03d}'
                    members = [genes[int(j)] for j in block]
                    gaps = [members[i].start - members[i-1].end for i in range(1, len(members))]
                    for j in block:
                        memberships[int(j)] = cid
                    for j, gap in zip(block[1:], gaps):
                        link_gaps[int(j)] = gap
                    catalog.append(dict(**identity, cluster_id=cid, scaffold=scaffold, n_genes=len(block),
                        start=min(g.start for g in members), end=max(g.end for g in members),
                        span_bp=max(g.end for g in members)-min(g.start for g in members)+1,
                        maximum_link_gap_bp=max(gaps), gene_ids=';'.join(g.id for g in members)))
                if len(genes) > 1:
                    blocks.append((pools[key[0]][scaffold], len(genes)))
                    observed_distances.extend(distances)
                else:
                    singletons += 1
                for j, g in enumerate(genes):
                    h = genes[neighbors[j]] if len(genes) > 1 else None
                    details.append(dict(**identity, gene_id=g.id, scaffold=scaffold, start=g.start, end=g.end,
                        nearest_gene_id=h.id if h else '', nearest_distance_bp=int(distances[j]) if h else '',
                        overlap=bool(max(g.start, h.start) <= min(g.end, h.end)) if h else '',
                        status='evaluable' if h else 'no_partner_on_scaffold', cluster_id=memberships.get(j, ''),
                        previous_group_gene_gap_bp=link_gaps.get(j, '')))
            row = dict(**identity, n_genes=len(observed), n_scaffolds=len(by_scaffold), n_without_partner=singletons,
                n_evaluable=len(observed_distances), counts_per_scaffold=json.dumps({s: len(gs) for s, gs in sorted(by_scaffold.items())}),
                threshold_bp=args.threshold, iterations=0, status='no_pairs_on_scaffold',
                p_proximity='', q_proximity_bh='', p_clustering='', q_clustering_bh='',
                observed_n_clusters=cluster_number, observed_n_clustered_genes=clustered_count,
                observed_clustered_fraction=clustered_count / len(observed))
            if observed_distances:
                obs = stats(np.array(observed_distances), args.threshold)
                obs['clustered_fraction'] = row['observed_clustered_fraction']
                row.update({f'observed_{k}': v for k, v in obs.items()})
                rng = group_rng(args.seed, key)
                simulations = np.empty((args.iterations, len(METRICS)))
                exceedances = 0
                last_report = time.monotonic()
                for iteration in range(args.iterations):
                    values, nc, ng = [], 0, 0
                    for universe, n in blocks:
                        indices = rng.choice(len(universe), size=n, replace=False)
                        selected = universe[indices]
                        d, _ = nearest(selected)
                        values.extend(d)
                        c, g = cluster_counts(selected, args.threshold)
                        nc, ng = nc + c, ng + g
                    sim = stats(np.array(values), args.threshold)
                    sim['clustered_fraction'] = ng / len(observed)
                    simulations[iteration] = [sim[m] for m in METRICS]
                    exceedances += int(ng >= clustered_count)
                    if sim_file:
                        writer.writerow(dict(**identity, iteration=iteration+1, **sim, n_clusters=nc, n_clustered_genes=ng))
                    if time.monotonic() - last_report > 10:
                        print(f'{key}: {iteration+1}/{args.iterations} draws', flush=True)
                        last_report = time.monotonic()
                for j, metric in enumerate(METRICS):
                    col = simulations[:, j]
                    row[f'null_{metric}_mean'] = float(col.mean())
                    for label, pct in [('p025', 2.5), ('p50', 50), ('p975', 97.5)]:
                        row[f'null_{metric}_{label}'] = float(np.percentile(col, pct))
                medians = simulations[:, METRICS.index('median_bp')]
                null_median = float(np.median(medians))
                row.update(status='tested', iterations=args.iterations,
                    p_proximity=(1 + int(np.sum(medians <= obs['median_bp']))) / (args.iterations + 1),
                    p_clustering=(1 + exceedances) / (args.iterations + 1),
                    observed_to_null_median_ratio=obs['median_bp'] / null_median if null_median else '',
                    null_median_variable=bool(np.ptp(medians) > 0),
                    null_clustering_variable=bool(np.ptp(simulations[:, -1]) > 0))
            summary.append(row)
            print(f'[{number}/{len(groups)}] {key}: {row["status"]}; {clustered_count}/{len(observed)} clustered.', flush=True)
            if sim_file:
                sim_file.flush()
    finally:
        if sim_file:
            sim_file.close()
    adjust_bh(summary, 'p_proximity', 'q_proximity_bh')
    adjust_bh(summary, 'p_clustering', 'q_clustering_bh')
    write_csv(out / 'summary.csv', summary, SUMMARY_FIELDS)
    write_csv(out / 'observed_genes.csv', details, DETAIL_FIELDS)
    write_csv(out / 'observed_clusters.csv', catalog, CLUSTER_FIELDS)
    if args.xlsx:
        save_xlsx(out / 'summary.xlsx', summary)
    metadata.update(completed=True, seconds=time.monotonic()-started, n_groups=len(summary),
                    n_tests=sum(r['status'] == 'tested' for r in summary))
    save_metadata()
    print(f'Completed. Results: {out}', flush=True)


if __name__ == '__main__':
    try:
        main()
    except (ValueError, FileNotFoundError, FileExistsError, KeyError) as exc:
        raise SystemExit(f'ERROR: {exc}') from exc
