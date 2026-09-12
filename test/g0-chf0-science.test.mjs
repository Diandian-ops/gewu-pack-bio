import assert from 'node:assert/strict';
import fs from 'node:fs';
import os from 'node:os';
import path from 'node:path';
import test from 'node:test';
import { spawnSync } from 'node:child_process';

const volcanoScript = 'plugins/volcano-plot/scripts/volcano-plot.R';
const pcaScript = 'plugins/pca-explorer/scripts/pca-explorer.py';

function parseCsv(text) {
  const lines = text.trim().split(/\r?\n/);
  const parse = (line) => line.split(',').map((value) => value.replace(/^"|"$/g, ''));
  const header = parse(lines[0]);
  return lines.slice(1).map((line) => Object.fromEntries(parse(line).map((value, index) => [header[index], value])));
}

function runVolcano(csv, params = {}) {
  const root = fs.mkdtempSync(path.join(os.tmpdir(), 'biof3-chf0-volcano-'));
  const table = path.join(root, 'deg.csv');
  fs.writeFileSync(table, csv);
  fs.writeFileSync(path.join(root, 'params.json'), JSON.stringify({
    deg_table: table,
    fc_threshold: 1,
    pval_threshold: 1.3,
    ...params,
  }));
  fs.copyFileSync(
    'plugins/volcano-plot/volcano-plot-report-template.html',
    path.join(root, 'report-template.html'),
  );
  const result = spawnSync('Rscript', [volcanoScript, root], { encoding: 'utf8' });
  return {
    root,
    result,
    rows: fs.existsSync(path.join(root, 'classified_table.csv'))
      ? parseCsv(fs.readFileSync(path.join(root, 'classified_table.csv'), 'utf8'))
      : [],
    report: fs.existsSync(path.join(root, 'report.html'))
      ? fs.readFileSync(path.join(root, 'report.html'), 'utf8')
      : '',
    cleanup: () => fs.rmSync(root, { recursive: true, force: true }),
  };
}

function runPca(exprCsv, groupCsv = null) {
  const root = fs.mkdtempSync(path.join(os.tmpdir(), 'biof3-chf0-pca-'));
  const expr = path.join(root, 'expr.csv');
  fs.writeFileSync(expr, exprCsv);
  const params = { expr_matrix: expr };
  if (groupCsv !== null) {
    const groups = path.join(root, 'groups.csv');
    fs.writeFileSync(groups, groupCsv);
    params.group_file = groups;
  }
  fs.writeFileSync(path.join(root, 'params.json'), JSON.stringify(params));
  const result = spawnSync('python', [pcaScript, root], { encoding: 'utf8' });
  const manifest = fs.existsSync(path.join(root, 'manifest.json'))
    ? JSON.parse(fs.readFileSync(path.join(root, 'manifest.json'), 'utf8'))
    : null;
  return { root, result, manifest, cleanup: () => fs.rmSync(root, { recursive: true, force: true }) };
}

test('CHF0 SCI-1 default and explicit padj agree without overwriting raw pvalue', () => {
  const csv = 'gene,log2FoldChange,padj,pvalue\ng1,2,0.01,0.9\ng2,-2,0.02,0.8\ng3,2,0.9,0.001\n';
  const implicit = runVolcano(csv);
  const explicit = runVolcano(csv, { pval_col: 'padj' });
  try {
    assert.equal(implicit.result.status, 0, implicit.result.stderr);
    assert.equal(explicit.result.status, 0, explicit.result.stderr);
    assert.deepEqual(implicit.rows.map((row) => row.sig), explicit.rows.map((row) => row.sig));
    assert.deepEqual(implicit.rows.map((row) => row.pvalue), ['0.9', '0.8', '0.001']);
    assert.deepEqual(implicit.rows.map((row) => row.padj), ['0.01', '0.02', '0.9']);
    assert.deepEqual(implicit.rows.map((row) => row.significance), ['0.01', '0.02', '0.9']);
    assert.ok(implicit.rows.every((row) => row.significance_column === 'padj'));
    assert.match(implicit.report, /-log10\(padj\)/);
  } finally {
    implicit.cleanup();
    explicit.cleanup();
  }
});

test('CHF0 SCI-1 explicit raw pvalue is a working, honestly labelled scientific change', () => {
  const run = runVolcano('gene,log2FoldChange,padj,pvalue\ng1,2,0.9,0.01\ng2,2,0.01,0.9\n', { pval_col: 'pvalue' });
  try {
    assert.equal(run.result.status, 0, run.result.stderr);
    assert.deepEqual(run.rows.map((row) => row.sig), ['Up', 'NS']);
    assert.deepEqual(run.rows.map((row) => row.pvalue), ['0.01', '0.9']);
    assert.deepEqual(run.rows.map((row) => row.padj), ['0.9', '0.01']);
    assert.ok(run.rows.every((row) => row.significance_column === 'pvalue'));
    assert.match(run.report, /-log10\(pvalue\)/);
  } finally {
    run.cleanup();
  }
});

test('CHF0 SCI-1 selected adjusted column is sufficient and invalid values fail closed', () => {
  const padjOnly = runVolcano('gene,log2FoldChange,padj\ng1,2,0.01\n');
  const missing = runVolcano('gene,log2FoldChange,pvalue\ng1,2,0.01\n');
  const bad = runVolcano('gene,log2FoldChange,padj,pvalue\ng1,2,bad,0.01\n');
  const outOfRange = runVolcano('gene,log2FoldChange,padj,pvalue\ng1,2,1.2,0.01\n');
  try {
    assert.equal(padjOnly.result.status, 0, padjOnly.result.stderr);
    assert.equal(missing.result.status, 1);
    assert.match(missing.result.stderr, /padj/);
    assert.equal(bad.result.status, 1);
    assert.match(bad.result.stderr, /非数值/);
    assert.equal(outOfRange.result.status, 1);
    assert.match(outOfRange.result.stderr, /0.*1/);
  } finally {
    padjOnly.cleanup(); missing.cleanup(); bad.cleanup(); outOfRange.cleanup();
  }
});

test('CHF0 SCI-1 descriptor and report bind the selected significance field', () => {
  const descriptor = JSON.parse(fs.readFileSync('plugins/volcano-plot/tool-definition.json', 'utf8'));
  const scatter = descriptor.resultStudio.interactionMappings.find((item) => item.artifactId === 'volcano');
  assert.equal(scatter.yField, 'significance');
  assert.equal(descriptor.resultStudio.entity.fields.significance, 'number');
  assert.equal(descriptor.resultStudio.entity.fields.significance_column, 'string');
  assert.equal(descriptor.resultStudio.parameterImpacts.pval_col, 'scientific');
});

test('CHF0 SCI-2 PCA without metadata remains explicitly unassigned', () => {
  const run = runPca('sample,g1,g2\ns1,1,2\ns2,2,1\ns3,3,4\n');
  try {
    assert.equal(run.result.status, 0, run.result.stderr);
    assert.equal(run.manifest.summary.grouping, 'unassigned');
  } finally { run.cleanup(); }
});

test('CHF0 SCI-3 PCA aligns scrambled metadata by ID', () => {
  const run = runPca(
    'sample,g1,g2\ns1,1,2\ns2,2,1\ns3,3,4\n',
    'sample,group\ns3,C\ns1,A\ns2,B\n',
  );
  try {
    assert.equal(run.result.status, 0, run.result.stderr);
    assert.equal(run.manifest.summary.grouping, 'matched_by_sample_id');
    assert.deepEqual(run.manifest.summary.sample_groups, { s1: 'A', s2: 'B', s3: 'C' });
  } finally { run.cleanup(); }
});

test('CHF0 SCI-3 PCA rejects missing, blank, duplicate, missing-set and extra sample IDs', () => {
  const cases = [
    { expr: [null, 's2'], samples: ['s1', 's2'], groups: ['A', 'B'] },
    { expr: ['', 's2'], samples: ['s1', 's2'], groups: ['A', 'B'] },
    { expr: ['s1', 's1'], samples: ['s1'], groups: ['A'] },
    { expr: ['s1', 's2'], samples: ['s1', 's1'], groups: ['A', 'B'] },
    { expr: ['s1', 's2', 's3'], samples: ['s1', 's2', 's4'], groups: ['A', 'B', 'C'] },
  ];
  const validation = spawnSync('python', [pcaScript, '--validate-sample-ids-json'], {
    encoding: 'utf8',
    input: JSON.stringify({ cases }),
    env: { ...process.env, PYTHONDONTWRITEBYTECODE: '1' },
  });
  assert.equal(validation.status, 0, validation.stderr);
  assert.deepEqual(JSON.parse(validation.stdout).map((item) => item.valid), [false, false, false, false, false]);

  const pcaSource = fs.readFileSync(pcaScript, 'utf8');
  assert.match(pcaSource, /"sample" not in gdf\.columns/);
  assert.match(pcaSource, /align_sample_groups/);

  const descriptor = JSON.parse(fs.readFileSync('plugins/pca-explorer/tool-definition.json', 'utf8'));
  const groupInput = descriptor.inputs.find((item) => item.id === 'group_file');
  assert.match(groupInput.hint, /未分组/);
  assert.doesNotMatch(groupInput.hint, /均分|A\/B/);
});
