"""Numerical and end-to-end checks; run python -m unittest discover -s tests -v."""
import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import openpyxl

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import gene_clustering as g


class AnalysisTests(unittest.TestCase):
    def test_nearest_against_scalar_oracle(self):
        rng = np.random.default_rng(42)
        for n in [2, 3, 15]:
            starts = rng.integers(1, 1000, n)
            ends = starts + rng.integers(1, 700, n)
            coords = np.column_stack([starts, ends])
            d, _ = g.nearest(coords)
            expected = [min(max(0, max(int(starts[i]), int(starts[j])) -
                                      min(int(ends[i]), int(ends[j])) - 1)
                            for j in range(n) if i != j) for i in range(n)]
            np.testing.assert_array_equal(d, expected)

    def test_gap_boundary_and_nested_intervals(self):
        # Clustering counts start-end=35000 as linked, 35001 as split.
        self.assertEqual(g.cluster_counts(np.array([[1, 100], [35100, 35200]]), 35000), (1, 2))
        coords = np.array([[1, 100], [35101, 35200]])
        self.assertEqual(g.cluster_counts(coords, 35000), (0, 0))
        self.assertEqual(g.nearest(coords)[0].tolist(), [35000, 35000])
        # Documents immediate-predecessor behavior for nested intervals:
        nested = np.array([[1, 100000], [100, 200], [50000, 50100]])
        self.assertEqual(g.cluster_counts(nested, 35000), (1, 2))
        self.assertEqual(g.nearest(nested)[0].tolist(), [0, 0, 0])
        self.assertEqual(g.cluster_counts(np.array([[1, 2]]), 35000), (0, 0))

    def test_bh(self):
        rows = [dict(status='tested', p=p) for p in [.04, .001, .02, 1]]
        rows.append(dict(status='no_pairs_on_scaffold', p=''))
        g.adjust_bh(rows, 'p', 'q')
        np.testing.assert_allclose([r['q'] for r in rows[:-1]], [.04*4/3, .004, .04, 1])
        self.assertNotIn('q', rows[-1])

    def test_gff_filters_and_multipart(self):
        lines = [
            's\t.\tgene\t100\t200\t.\t+\t.\tID=gene-a;gene_biotype=protein_coding;part=1;exception=trans-splicing',
            's\t.\tgene\t300\t400\t.\t-\t.\tID=gene-a;gene_biotype=protein_coding;part=2;exception=trans-splicing',
            's\t.\tgene\t500\t600\t.\t+\t.\tID=gene-p;gene_biotype=protein_coding;pseudo=true',
            's\t.\tgene\t700\t800\t.\t+\t.\tID=gene-n;gene_biotype=lncRNA',
            'mt\t.\tgene\t1\t90\t.\t+\t.\tID=gene-m;gene_biotype=protein_coding',
            'mt\t.\tregion\t1\t1000\t.\t+\t.\tID=mt;genome=mitochondrion']
        with tempfile.TemporaryDirectory() as t:
            path = Path(t)/'test.gff'
            path.write_text('\n'.join(lines)+'\n')
            genes, _, audit = g.read_gff(path)
            self.assertEqual(set(genes), {'gene-a'})
            self.assertEqual((genes['gene-a'].start, genes['gene-a'].end), (100, 400))
            self.assertEqual(audit['organelle_genes_excluded'], 1)
            path.write_text('\n'.join(lines + [lines[0].replace('100\t200', '900\t950').replace(';part=1;exception=trans-splicing','')]))
            with self.assertRaises(ValueError):
                g.read_gff(path)

    def test_cli_outputs_and_filtered_reproducibility(self):
        with tempfile.TemporaryDirectory() as t:
            base = [sys.executable, str(ROOT/'gene_clustering.py'), '--genes', str(ROOT/'examples/genes.csv'),
                    '--gff-dir', str(ROOT/'examples'), '--iterations', '99']
            full, filtered = Path(t)/'full', Path(t)/'filtered'
            subprocess.run(base+['--output', str(full), '--xlsx'], check=True, capture_output=True)
            subprocess.run(base+['--output', str(filtered), '--group', 'DEMO_CLOSE'], check=True, capture_output=True)
            def read(path):
                with path.open(encoding='utf-8-sig',newline='') as f:
                    return list(csv.DictReader(f))
            rows = read(full/'summary.csv')
            self.assertEqual(len(rows), 3)
            self.assertEqual(list(rows[0]), g.SUMMARY_FIELDS)
            sims = read(full/'simulation_statistics.csv')
            self.assertEqual(len(sims), 198)
            self.assertEqual([r for r in sims if r['group']=='DEMO_CLOSE'], read(filtered/'simulation_statistics.csv'))
            for r in rows:
                self.assertEqual(int(r['n_evaluable'])+int(r['n_without_partner']), int(r['n_genes']))
                if r['status'] != 'tested':
                    self.assertEqual(r['p_clustering'], '')
                    continue
                subset = [s for s in sims if s['group']==r['group']]
                lower = sum(float(s['median_bp']) <= float(r['observed_median_bp']) for s in subset)
                upper = sum(int(s['n_clustered_genes']) >= int(r['observed_n_clustered_genes']) for s in subset)
                self.assertAlmostEqual(float(r['p_proximity']), (lower+1)/100)
                self.assertAlmostEqual(float(r['p_clustering']), (upper+1)/100)
                fractions = [float(s['clustered_fraction']) for s in subset]
                self.assertAlmostEqual(float(r['null_clustered_fraction_p975']), np.percentile(fractions,97.5))
            spread = next(r for r in rows if r['group']=='DEMO_SPREAD')
            self.assertEqual(float(spread['p_clustering']), 1)
            close = next(r for r in rows if r['group']=='DEMO_CLOSE')
            self.assertEqual(float(close['observed_clustered_fraction']), .75)
            self.assertEqual(float(close['observed_within_threshold_fraction']), 1)
            metadata=json.loads((full/'parameters.json').read_text())
            self.assertTrue(metadata['completed'])
            wb=openpyxl.load_workbook(full/'summary.xlsx',read_only=True,data_only=True)
            self.assertIsInstance(wb['Results']['E2'].value,int)
            wb.close()
            retry=subprocess.run(base+['--output',str(full)],capture_output=True)
            self.assertNotEqual(retry.returncode,0)


if __name__ == '__main__':
    unittest.main()
